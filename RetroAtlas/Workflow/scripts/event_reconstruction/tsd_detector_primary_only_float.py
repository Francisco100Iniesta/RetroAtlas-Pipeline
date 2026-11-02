#!/usr/bin/env python3
"""
PRIMARY-only TSD detector.

This script preserves the PRIMARY/v2 TSD caller from tsd_detector_v5.py and
removes all expensive rescue/refinement logic and all poly(A/T) analysis.

PRIMARY TSD calling:
  1. read BLAST outfmt6 HSPs;
  2. build same-strand, non-overlapping flanking anchor pairs;
  3. infer TSD length from overlap on the reference axis;
  4. validate both contig copies against the SAME reference-overlap sequence;
  5. require consistency with the reported insertion size.

No Smith-Waterman.
No parasail.
No local refinement.
No terminal rescue.
No ALT FASTA.
No poly(A/T).

Single-event mode accepts --ref-fasta and --contig-fasta only for command-line
compatibility with older pipelines. They are ignored.

Batch mode can process independent events in parallel using only the Python
standard library (ProcessPoolExecutor).
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import re
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict, dataclass
from typing import Dict, List, Optional, Tuple


DNA = set("ACGT")
_COMP = str.maketrans("ACGTNacgtn-", "TGCANtgcan-")


def comp_base(base: str) -> str:
    return base.translate(_COMP)


@dataclass
class HSP:
    qseqid: str
    sseqid: str
    pident: float
    length: int
    mismatch: int
    gapopen: int
    qstart: int
    qend: int
    sstart: int
    send: int
    evalue: float
    bitscore: float


@dataclass
class Anchor:
    hsp: HSP
    contig_start: int
    contig_end: int
    ref_at_contig_start: int
    ref_at_contig_end: int
    strand: str
    ref_axis_start: int
    ref_axis_end: int

    @property
    def length(self) -> int:
        return self.contig_end - self.contig_start + 1


@dataclass
class PairwiseBlock:
    qstart: int
    qend: int
    sstart: int
    send: int
    qmap: Dict[int, str]
    smap: Dict[int, str]
    q_to_sbase: Dict[int, str]
    s_to_qbase: Dict[int, str]


@dataclass
class PairCandidate:
    left: Anchor
    right: Anchor
    contig_gap: int
    reference_overlap: int
    predicted_net_gain: int
    size_residual: int
    reference_tsd_sequence: Optional[str] = None
    tsd_left_copy: Optional[str] = None
    tsd_right_copy: Optional[str] = None
    left_ref_identity: Optional[float] = None
    right_ref_identity: Optional[float] = None
    pair_identity: Optional[float] = None
    left_callable_fraction: Optional[float] = None
    right_callable_fraction: Optional[float] = None
    pair_callable_fraction: Optional[float] = None
    reference_consistency_identity: Optional[float] = None

    @property
    def min_ref_identity(self) -> Optional[float]:
        vals = [
            x
            for x in (self.left_ref_identity, self.right_ref_identity)
            if x is not None
        ]
        return min(vals) if len(vals) == 2 else None

    @property
    def min_callable_fraction(self) -> Optional[float]:
        vals = [
            x
            for x in (self.left_callable_fraction, self.right_callable_fraction)
            if x is not None
        ]
        return min(vals) if len(vals) == 2 else None

    @property
    def anchor_quality(self) -> float:
        return (
            self.left.hsp.bitscore
            + self.right.hsp.bitscore
            + 0.25 * (self.left.length + self.right.length)
        )


@dataclass
class Result:
    event_id: str
    ins_size: float
    tsd_detected: bool
    tsd_length: Optional[int]
    tsd_sequence: Optional[str]
    classification: str
    detection_mode: str
    inserted_core_length: Optional[int]
    reference_overlap: Optional[int]
    predicted_net_gain: Optional[int]
    size_residual: Optional[int]
    tsd_left_copy: Optional[str]
    tsd_right_copy: Optional[str]
    tsd_identity_callable: Optional[float]
    tsd_callable_fraction: Optional[float]
    tsd_left_ref_identity: Optional[float]
    tsd_right_ref_identity: Optional[float]
    tsd_left_callable_fraction: Optional[float]
    tsd_right_callable_fraction: Optional[float]
    reference_consistency_identity: Optional[float]
    strand: Optional[str]
    left_contig_start: Optional[int]
    left_contig_end: Optional[int]
    right_contig_start: Optional[int]
    right_contig_end: Optional[int]
    left_ref_start: Optional[int]
    left_ref_end: Optional[int]
    right_ref_start: Optional[int]
    right_ref_end: Optional[int]
    n_hsps: int
    n_anchor_pairs: int
    n_valid_tsd_pairs: int
    ambiguous_pair: bool
    diagnostic: str


def read_blast_outfmt6(path: str) -> List[HSP]:
    hsps: List[HSP] = []

    with open(path, "r", encoding="utf-8") as fh:
        for line_no, line in enumerate(fh, 1):
            line = line.strip()

            if not line or line.startswith("#"):
                continue

            fields = line.split("\t")
            if len(fields) < 12:
                fields = line.split()

            if len(fields) < 12:
                raise ValueError(
                    f"{path}:{line_no}: expected >=12 BLAST columns"
                )

            hsps.append(
                HSP(
                    qseqid=fields[0],
                    sseqid=fields[1],
                    pident=float(fields[2]),
                    length=int(fields[3]),
                    mismatch=int(fields[4]),
                    gapopen=int(fields[5]),
                    qstart=int(fields[6]),
                    qend=int(fields[7]),
                    sstart=int(fields[8]),
                    send=int(fields[9]),
                    evalue=float(fields[10]),
                    bitscore=float(fields[11]),
                )
            )

    return hsps


def detect_contig_side(
    hsps: List[HSP],
    requested: str,
    contig_id: Optional[str],
) -> str:
    if requested in {"query", "subject"}:
        return requested

    if not hsps:
        raise ValueError("No HSPs available to infer contig side")

    qids = {h.qseqid for h in hsps}
    sids = {h.sseqid for h in hsps}

    if contig_id:
        in_q = contig_id in qids
        in_s = contig_id in sids

        if in_q ^ in_s:
            return "query" if in_q else "subject"

        raise ValueError(
            f"Cannot uniquely locate contig-id {contig_id!r}"
        )

    q_like = [
        x
        for x in qids
        if x.lower() in {"contig", "seq"} or "contig" in x.lower()
    ]
    s_like = [
        x
        for x in sids
        if x.lower() in {"contig", "seq"} or "contig" in x.lower()
    ]

    if bool(q_like) ^ bool(s_like):
        return "query" if q_like else "subject"

    raise ValueError(
        "Could not infer contig side; use --contig-side query|subject"
    )


def normalize_anchor(h: HSP, contig_side: str) -> Anchor:
    if contig_side == "query":
        ep1, ep2 = (h.qstart, h.sstart), (h.qend, h.send)
    else:
        ep1, ep2 = (h.sstart, h.qstart), (h.send, h.qend)

    (c1, r1), (c2, r2) = sorted(
        (ep1, ep2),
        key=lambda x: x[0],
    )

    strand = "+" if r2 >= r1 else "-"
    mult = 1 if strand == "+" else -1

    return Anchor(
        hsp=h,
        contig_start=c1,
        contig_end=c2,
        ref_at_contig_start=r1,
        ref_at_contig_end=r2,
        strand=strand,
        ref_axis_start=mult * r1,
        ref_axis_end=mult * r2,
    )


def anchor_ref_bounds(a: Anchor) -> Tuple[int, int]:
    return (
        min(a.ref_at_contig_start, a.ref_at_contig_end),
        max(a.ref_at_contig_start, a.ref_at_contig_end),
    )


def parse_pairwise_blocks(path: str) -> List[PairwiseBlock]:
    score_re = re.compile(r"^\s*Score\s*=")
    seq_re = re.compile(
        r"^(Query|Sbjct)\s+(\d+)\s+([A-Za-z-]+)\s+(\d+)\s*$"
    )

    with open(
        path,
        "r",
        encoding="utf-8",
        errors="replace",
    ) as fh:
        lines = list(fh)

    blocks: List[PairwiseBlock] = []
    cur = None
    i = 0

    while i < len(lines):
        line = lines[i].rstrip("\n")

        if score_re.match(line):
            if cur and cur["qmap"]:
                blocks.append(PairwiseBlock(**cur))

            cur = {
                "qstart": None,
                "qend": None,
                "sstart": None,
                "send": None,
                "qmap": {},
                "smap": {},
                "q_to_sbase": {},
                "s_to_qbase": {},
            }

            i += 1
            continue

        mq = seq_re.match(line)

        if mq and mq.group(1) == "Query" and cur is not None:
            qs = int(mq.group(2))
            qseq = mq.group(3).upper()
            qe = int(mq.group(4))

            ms = None
            j = i + 1

            while j < min(i + 4, len(lines)):
                m = seq_re.match(lines[j].rstrip("\n"))

                if m and m.group(1) == "Sbjct":
                    ms = m
                    break

                j += 1

            if ms is None:
                i += 1
                continue

            ss = int(ms.group(2))
            sseq = ms.group(3).upper()
            se = int(ms.group(4))

            if cur["qstart"] is None:
                cur["qstart"] = qs
                cur["sstart"] = ss

            cur["qend"] = qe
            cur["send"] = se

            qstep = 1 if qe >= qs else -1
            sstep = 1 if se >= ss else -1
            qp = qs
            sp = ss

            for qb, sb in zip(qseq, sseq):
                qpos = None
                spos = None

                if qb != "-":
                    qpos = qp
                    cur["qmap"][qp] = qb
                    qp += qstep

                if sb != "-":
                    spos = sp
                    cur["smap"][sp] = sb
                    sp += sstep

                if qpos is not None:
                    cur["q_to_sbase"][qpos] = (
                        sb if sb != "-" else "-"
                    )

                if spos is not None:
                    cur["s_to_qbase"][spos] = (
                        qb if qb != "-" else "-"
                    )

            i = j + 1
            continue

        i += 1

    if cur and cur["qmap"]:
        blocks.append(PairwiseBlock(**cur))

    return blocks


def hsp_block_key(
    h: HSP,
) -> Tuple[int, int, int, int]:
    return (
        h.qstart,
        h.qend,
        h.sstart,
        h.send,
    )


def _canonical_base(
    displayed: str,
    reference_forward: bool,
) -> str:
    if displayed == "-":
        return "-"

    return (
        displayed
        if reference_forward
        else comp_base(displayed)
    )


def overlap_metrics(
    block: PairwiseBlock,
    ref_side: str,
    ref_lo: int,
    ref_hi: int,
) -> Tuple[Optional[float], float, str, str]:
    matches = 0
    callable_n = 0
    ref_out: List[str] = []
    copy_out: List[str] = []

    if ref_side == "query":
        reference_forward = block.qend >= block.qstart

        for pos in range(ref_lo, ref_hi + 1):
            rb = block.qmap.get(pos, "N").upper()
            cb = block.q_to_sbase.get(pos, "N").upper()

            ref_out.append(
                _canonical_base(rb, reference_forward)
            )
            copy_out.append(
                _canonical_base(cb, reference_forward)
            )

            if rb in DNA and cb in DNA:
                callable_n += 1
                matches += int(rb == cb)

    else:
        reference_forward = block.send >= block.sstart

        for pos in range(ref_lo, ref_hi + 1):
            rb = block.smap.get(pos, "N").upper()
            cb = block.s_to_qbase.get(pos, "N").upper()

            ref_out.append(
                _canonical_base(rb, reference_forward)
            )
            copy_out.append(
                _canonical_base(cb, reference_forward)
            )

            if rb in DNA and cb in DNA:
                callable_n += 1
                matches += int(rb == cb)

    total = ref_hi - ref_lo + 1

    identity = (
        matches / callable_n
        if callable_n
        else None
    )

    callable_fraction = (
        callable_n / total
        if total > 0
        else 0.0
    )

    return (
        identity,
        callable_fraction,
        "".join(ref_out),
        "".join(copy_out),
    )


def compare_callable(
    a: str,
    b: str,
) -> Tuple[Optional[float], float]:
    n = min(len(a), len(b))

    if n == 0:
        return None, 0.0

    pairs = [
        (x, y)
        for x, y in zip(
            a[:n].upper(),
            b[:n].upper(),
        )
        if x in DNA and y in DNA
    ]

    identity = (
        sum(x == y for x, y in pairs) / len(pairs)
        if pairs
        else None
    )

    callable_fraction = len(pairs) / n

    return identity, callable_fraction


def generate_pairs(
    anchors: List[Anchor],
    min_anchor: int,
    min_anchor_identity: float,
) -> List[PairCandidate]:
    """
    PRIMARY/v2 pairing only.

    Candidate anchors must:
      * pass minimum anchor length;
      * pass minimum BLAST identity;
      * have the same strand;
      * be non-overlapping on the contig.
    """
    good = [
        a
        for a in anchors
        if (
            a.length >= min_anchor
            and a.hsp.pident >= min_anchor_identity
        )
    ]

    good.sort(
        key=lambda a: (
            a.contig_start,
            a.contig_end,
        )
    )

    out: List[PairCandidate] = []

    for i, left in enumerate(good):
        for right in good[i + 1:]:
            if left.strand != right.strand:
                continue

            gap = (
                right.contig_start
                - left.contig_end
                - 1
            )

            if gap < 0:
                continue

            overlap = (
                left.ref_axis_end
                - right.ref_axis_start
                + 1
            )

            out.append(
                PairCandidate(
                    left=left,
                    right=right,
                    contig_gap=gap,
                    reference_overlap=overlap,
                    predicted_net_gain=gap + overlap,
                    size_residual=0,
                )
            )

    return out


def enrich_candidate(
    p: PairCandidate,
    ins_size: float,
    block_map: Dict[
        Tuple[int, int, int, int],
        PairwiseBlock,
    ],
    contig_side: str,
) -> None:
    p.size_residual = (
        p.predicted_net_gain
        - ins_size
    )

    if p.reference_overlap <= 0:
        return

    llo, lhi = anchor_ref_bounds(p.left)
    rlo, rhi = anchor_ref_bounds(p.right)

    ref_lo = max(llo, rlo)
    ref_hi = min(lhi, rhi)

    if ref_hi < ref_lo:
        return

    bl = block_map.get(
        hsp_block_key(p.left.hsp)
    )
    br = block_map.get(
        hsp_block_key(p.right.hsp)
    )

    if bl is None or br is None:
        return

    ref_side = (
        "subject"
        if contig_side == "query"
        else "query"
    )

    li, lc, lref, lcopy = overlap_metrics(
        bl,
        ref_side,
        ref_lo,
        ref_hi,
    )

    ri, rc, rref, rcopy = overlap_metrics(
        br,
        ref_side,
        ref_lo,
        ref_hi,
    )

    pair_id, pair_call = compare_callable(
        lcopy,
        rcopy,
    )

    ref_id, _ = compare_callable(
        lref,
        rref,
    )

    p.left_ref_identity = li
    p.right_ref_identity = ri
    p.left_callable_fraction = lc
    p.right_callable_fraction = rc
    p.pair_identity = pair_id
    p.pair_callable_fraction = pair_call
    p.reference_consistency_identity = ref_id
    p.reference_tsd_sequence = lref
    p.tsd_left_copy = lcopy
    p.tsd_right_copy = rcopy


def candidate_valid(
    p: PairCandidate,
    max_tsd: int,
    min_tsd_identity: float,
    min_callable_fraction: float,
    size_tolerance: int,
) -> bool:
    return bool(
        1 <= p.reference_overlap <= max_tsd
        and abs(p.size_residual) <= size_tolerance
        and p.min_ref_identity is not None
        and p.min_ref_identity >= min_tsd_identity
        and p.min_callable_fraction is not None
        and p.min_callable_fraction >= min_callable_fraction
        and p.reference_consistency_identity is not None
        and p.reference_consistency_identity >= 0.99
    )


def rank_key(
    p: PairCandidate,
):
    return (
        abs(p.size_residual),
        -(
            p.min_ref_identity
            if p.min_ref_identity is not None
            else -1
        ),
        -(
            p.pair_identity
            if p.pair_identity is not None
            else -1
        ),
        -p.anchor_quality,
        p.reference_overlap,
    )


def pair_ambiguous(
    a: PairCandidate,
    b: PairCandidate,
    identity_margin: float,
) -> bool:
    if abs(a.size_residual) != abs(b.size_residual):
        return False

    ai = (
        a.min_ref_identity
        if a.min_ref_identity is not None
        else -1
    )
    bi = (
        b.min_ref_identity
        if b.min_ref_identity is not None
        else -1
    )

    if abs(ai - bi) > identity_margin:
        return False

    denom = max(
        abs(a.anchor_quality),
        1.0,
    )

    return (
        abs(
            a.anchor_quality
            - b.anchor_quality
        )
        / denom
        <= 0.05
    )


def choose_diagnostic_pair(
    pairs: List[PairCandidate],
    max_tsd: int,
) -> PairCandidate:
    def key(p: PairCandidate):
        len_penalty = (
            0
            if 1 <= p.reference_overlap <= max_tsd
            else 1
        )

        seq_identity = (
            p.min_ref_identity
            if p.min_ref_identity is not None
            else -1
        )

        return (
            len_penalty,
            abs(p.size_residual),
            -seq_identity,
            -p.anchor_quality,
        )

    return sorted(
        pairs,
        key=key,
    )[0]


def empty_result(
    event_id: str,
    ins_size: float,
    n_hsps: int,
    classification: str,
    diagnostic: str,
) -> Result:
    return Result(
        event_id=event_id,
        ins_size=ins_size,
        tsd_detected=False,
        tsd_length=None,
        tsd_sequence=None,
        classification=classification,
        detection_mode="NONE",
        inserted_core_length=None,
        reference_overlap=None,
        predicted_net_gain=None,
        size_residual=None,
        tsd_left_copy=None,
        tsd_right_copy=None,
        tsd_identity_callable=None,
        tsd_callable_fraction=None,
        tsd_left_ref_identity=None,
        tsd_right_ref_identity=None,
        tsd_left_callable_fraction=None,
        tsd_right_callable_fraction=None,
        reference_consistency_identity=None,
        strand=None,
        left_contig_start=None,
        left_contig_end=None,
        right_contig_start=None,
        right_contig_end=None,
        left_ref_start=None,
        left_ref_end=None,
        right_ref_start=None,
        right_ref_end=None,
        n_hsps=n_hsps,
        n_anchor_pairs=0,
        n_valid_tsd_pairs=0,
        ambiguous_pair=False,
        diagnostic=diagnostic,
    )


def analyze_primary(
    blast_path: str,
    alignment_path: str,
    ins_size: float,
    event_id: str,
    contig_side_req: str,
    contig_id: Optional[str],
    min_anchor: int,
    min_anchor_identity: float,
    max_tsd: int,
    min_tsd_identity: float,
    min_callable_fraction: float,
    size_tolerance: int,
    ambiguity_identity_margin: float,
) -> Result:
    hsps = read_blast_outfmt6(
        blast_path
    )

    if not hsps:
        return empty_result(
            event_id,
            ins_size,
            0,
            "INSUFFICIENT_ALIGNMENT",
            "no_hsps",
        )

    contig_side = detect_contig_side(
        hsps,
        contig_side_req,
        contig_id,
    )

    anchors = [
        normalize_anchor(
            h,
            contig_side,
        )
        for h in hsps
    ]

    blocks = parse_pairwise_blocks(
        alignment_path
    )

    block_map = {
        (
            b.qstart,
            b.qend,
            b.sstart,
            b.send,
        ): b
        for b in blocks
    }

    pairs = generate_pairs(
        anchors,
        min_anchor,
        min_anchor_identity,
    )

    for p in pairs:
        enrich_candidate(
            p,
            ins_size,
            block_map,
            contig_side,
        )

    if not pairs:
        result = empty_result(
            event_id,
            ins_size,
            len(hsps),
            "INSUFFICIENT_ALIGNMENT",
            "no_nonoverlapping_flanking_pair_after_anchor_filters",
        )
        result.n_anchor_pairs = 0
        return result

    valid = [
        p
        for p in pairs
        if candidate_valid(
            p,
            max_tsd,
            min_tsd_identity,
            min_callable_fraction,
            size_tolerance,
        )
    ]

    valid.sort(
        key=rank_key
    )

    ambiguous = False

    if valid:
        best = valid[0]

        if len(valid) > 1:
            ambiguous = pair_ambiguous(
                valid[0],
                valid[1],
                ambiguity_identity_margin,
            )

        classification = (
            "AMBIGUOUS_ANCHORS"
            if ambiguous
            else "PASS_TSD"
        )

        detected = not ambiguous

    else:
        best = choose_diagnostic_pair(
            pairs,
            max_tsd,
        )

        detected = False

        if best.reference_overlap < 0:
            classification = (
                "INS_WITH_TARGET_SITE_DELETION"
            )

        elif best.reference_overlap == 0:
            classification = "NO_TSD"

        elif best.reference_overlap > max_tsd:
            classification = "OVERLAP_TOO_LONG"

        elif abs(best.size_residual) > size_tolerance:
            classification = (
                "POSSIBLE_TSD_SIZE_MISMATCH"
            )

        elif (
            best.min_callable_fraction is None
            or best.min_callable_fraction
            < min_callable_fraction
        ):
            classification = "AMBIGUOUS_SEQUENCE"

        elif (
            best.min_ref_identity is None
            or best.min_ref_identity
            < min_tsd_identity
        ):
            classification = (
                "POSSIBLE_TSD_SEQUENCE_MISMATCH"
            )

        else:
            classification = "POSSIBLE_TSD"

    llo, lhi = anchor_ref_bounds(
        best.left
    )
    rlo, rhi = anchor_ref_bounds(
        best.right
    )

    diag = [
        f"geometry_overlap={best.reference_overlap}",
        f"contig_gap={best.contig_gap}",
        f"predicted_gain={best.predicted_net_gain}",
        f"size_residual={best.size_residual}",
        f"valid_pairs={len(valid)}",
    ]

    if best.left_ref_identity is not None:
        diag.append(
            f"left_ref_identity={best.left_ref_identity:.3f}"
        )

    if best.right_ref_identity is not None:
        diag.append(
            f"right_ref_identity={best.right_ref_identity:.3f}"
        )

    if best.pair_identity is not None:
        diag.append(
            f"copy_identity={best.pair_identity:.3f}"
        )

    if best.min_callable_fraction is not None:
        diag.append(
            "min_callable_fraction="
            f"{best.min_callable_fraction:.3f}"
        )

    if ambiguous:
        diag.append(
            "multiple_near_equivalent_valid_pairs"
        )

    return Result(
        event_id=event_id,
        ins_size=ins_size,
        tsd_detected=detected,
        tsd_length=(
            best.reference_overlap
            if detected
            else None
        ),
        tsd_sequence=(
            best.reference_tsd_sequence
            if detected
            else None
        ),
        classification=classification,
        detection_mode=(
            "PRIMARY"
            if detected
            else "NONE"
        ),
        inserted_core_length=best.contig_gap,
        reference_overlap=best.reference_overlap,
        predicted_net_gain=best.predicted_net_gain,
        size_residual=best.size_residual,
        tsd_left_copy=best.tsd_left_copy,
        tsd_right_copy=best.tsd_right_copy,
        tsd_identity_callable=best.pair_identity,
        tsd_callable_fraction=best.pair_callable_fraction,
        tsd_left_ref_identity=best.left_ref_identity,
        tsd_right_ref_identity=best.right_ref_identity,
        tsd_left_callable_fraction=best.left_callable_fraction,
        tsd_right_callable_fraction=best.right_callable_fraction,
        reference_consistency_identity=best.reference_consistency_identity,
        strand=best.left.strand,
        left_contig_start=best.left.contig_start,
        left_contig_end=best.left.contig_end,
        right_contig_start=best.right.contig_start,
        right_contig_end=best.right.contig_end,
        left_ref_start=llo,
        left_ref_end=lhi,
        right_ref_start=rlo,
        right_ref_end=rhi,
        n_hsps=len(hsps),
        n_anchor_pairs=len(pairs),
        n_valid_tsd_pairs=len(valid),
        ambiguous_pair=ambiguous,
        diagnostic=";".join(diag),
    )


def fmt_value(value):
    if value is None:
        return "NA"

    if isinstance(value, bool):
        return (
            "TRUE"
            if value
            else "FALSE"
        )

    if isinstance(value, float):
        if math.isnan(value):
            return "NA"

        return f"{value:.6g}"

    return str(value)


def write_tsv(
    result: Result,
    output: str,
) -> None:
    parent = os.path.dirname(
        os.path.abspath(output)
    )

    if parent:
        os.makedirs(
            parent,
            exist_ok=True,
        )

    row = asdict(result)

    with open(
        output,
        "w",
        encoding="utf-8",
        newline="",
    ) as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=list(row.keys()),
            delimiter="\t",
        )

        writer.writeheader()

        writer.writerow(
            {
                k: fmt_value(v)
                for k, v in row.items()
            }
        )


def _tsd_cfg_from_args(
    a: argparse.Namespace,
) -> dict:
    return {
        "min_anchor": a.min_anchor,
        "min_anchor_identity": a.min_anchor_identity,
        "max_tsd": a.max_tsd,
        "min_tsd_identity": a.min_tsd_identity,
        "min_callable_fraction": a.min_callable_fraction,
        "size_tolerance": a.size_tolerance,
        "ambiguity_identity_margin": a.ambiguity_identity_margin,
    }


def process_event(
    *,
    blast: str,
    alignment: str,
    ins_size: float,
    event_id: str,
    output: str,
    contig_side: str,
    contig_id: Optional[str],
    tsd_cfg: dict,
) -> Result:
    result = analyze_primary(
        blast_path=blast,
        alignment_path=alignment,
        ins_size=ins_size,
        event_id=event_id,
        contig_side_req=contig_side,
        contig_id=contig_id,
        **tsd_cfg,
    )

    write_tsv(
        result,
        output,
    )

    return result


def _worker_process_event(
    payload: dict,
) -> Tuple[str, str, str, int]:
    try:
        result = process_event(
            **payload
        )

        return (
            result.event_id,
            result.classification,
            payload["output"],
            0,
        )

    except Exception as exc:
        return (
            str(
                payload.get(
                    "event_id",
                    "event",
                )
            ),
            str(exc),
            str(
                payload.get(
                    "output",
                    "",
                )
            ),
            2,
        )


def read_batch_manifest(
    path: str,
    default_contig_side: str,
    default_contig_id: Optional[str],
    tsd_cfg: dict,
) -> List[dict]:
    required = {
        "blast",
        "alignment",
        "ins_size",
        "output",
    }

    with open(
        path,
        "r",
        encoding="utf-8",
        newline="",
    ) as fh:
        reader = csv.DictReader(
            fh,
            delimiter="\t",
        )

        if reader.fieldnames is None:
            raise ValueError(
                f"{path}: empty manifest"
            )

        missing = (
            required
            - set(reader.fieldnames)
        )

        if missing:
            raise ValueError(
                f"{path}: missing required columns: "
                + ", ".join(
                    sorted(missing)
                )
            )

        payloads: List[dict] = []

        for line_no, row in enumerate(
            reader,
            2,
        ):
            if not any(
                (v or "").strip()
                for v in row.values()
            ):
                continue

            try:
                ins_size = float(
                    row["ins_size"]
                )
            except Exception as exc:
                raise ValueError(
                    f"{path}:{line_no}: "
                    f"invalid ins_size="
                    f"{row.get('ins_size')!r}"
                ) from exc

            event_id = (
                (
                    row.get("event_id")
                    or ""
                ).strip()
                or f"event_{line_no - 1}"
            )

            contig_side = (
                (
                    row.get("contig_side")
                    or ""
                ).strip()
                or default_contig_side
            )

            contig_id = (
                (
                    row.get("contig_id")
                    or ""
                ).strip()
                or default_contig_id
            )

            if contig_side not in {
                "auto",
                "query",
                "subject",
            }:
                raise ValueError(
                    f"{path}:{line_no}: "
                    "contig_side must be "
                    "auto|query|subject"
                )

            payloads.append(
                {
                    "blast": row["blast"],
                    "alignment": row["alignment"],
                    "ins_size": ins_size,
                    "event_id": event_id,
                    "output": row["output"],
                    "contig_side": contig_side,
                    "contig_id": contig_id,
                    "tsd_cfg": tsd_cfg,
                }
            )

    return payloads


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=(
            "PRIMARY-only TSD detection "
            "from BLAST geometry"
        )
    )

    # Single-event mode.
    p.add_argument("--blast")
    p.add_argument("--alignment")
    p.add_argument(
        "--ins-size",
        type=float,
    )
    p.add_argument(
        "--event-id",
        default="event",
    )
    p.add_argument("--output")
    p.add_argument(
        "--contig-side",
        choices=[
            "auto",
            "query",
            "subject",
        ],
        default="auto",
    )
    p.add_argument("--contig-id")

    # Compatibility only: accepted but never opened/read.
    p.add_argument(
        "--ref-fasta",
        help=(
            "Ignored. Kept only for "
            "compatibility with old commands."
        ),
    )
    p.add_argument(
        "--contig-fasta",
        help=(
            "Ignored. Kept only for "
            "compatibility with old commands."
        ),
    )

    # PRIMARY/v2 thresholds.
    p.add_argument(
        "--min-anchor",
        type=int,
        default=100,
    )
    p.add_argument(
        "--min-anchor-identity",
        type=float,
        default=90.0,
    )
    p.add_argument(
        "--max-tsd",
        type=int,
        default=50,
    )
    p.add_argument(
        "--min-tsd-identity",
        type=float,
        default=0.85,
    )
    p.add_argument(
        "--min-callable-fraction",
        type=float,
        default=0.80,
    )
    p.add_argument(
        "--size-tolerance",
        type=int,
        default=5,
    )
    p.add_argument(
        "--ambiguity-identity-margin",
        type=float,
        default=0.02,
    )

    # Multi-event / multi-core mode.
    p.add_argument(
        "--batch-manifest",
        help=(
            "TSV columns: "
            "blast, alignment, ins_size, output; "
            "optional event_id, contig_side, contig_id"
        ),
    )
    p.add_argument(
        "--jobs",
        type=int,
        default=1,
        help=(
            "Parallel event processes in batch mode; "
            "0 = all available CPUs"
        ),
    )

    return p


def run_single(
    a: argparse.Namespace,
) -> int:
    required = {
        "--blast": a.blast,
        "--alignment": a.alignment,
        "--ins-size": a.ins_size,
        "--output": a.output,
    }

    missing = [
        flag
        for flag, value
        in required.items()
        if value is None
    ]

    if missing:
        raise ValueError(
            "single-event mode requires: "
            + ", ".join(missing)
        )

    result = process_event(
        blast=a.blast,
        alignment=a.alignment,
        ins_size=a.ins_size,
        event_id=a.event_id,
        output=a.output,
        contig_side=a.contig_side,
        contig_id=a.contig_id,
        tsd_cfg=_tsd_cfg_from_args(a),
    )

    print(
        f"{result.event_id}: "
        f"{result.classification}; "
        f"mode={result.detection_mode}; "
        f"TSD={result.tsd_sequence or 'NA'}; "
        f"TSD_len="
        f"{result.tsd_length if result.tsd_length is not None else 'NA'}; "
        f"output={a.output}"
    )

    return 0


def run_batch(
    a: argparse.Namespace,
) -> int:
    tsd_cfg = _tsd_cfg_from_args(a)

    payloads = read_batch_manifest(
        path=a.batch_manifest,
        default_contig_side=a.contig_side,
        default_contig_id=a.contig_id,
        tsd_cfg=tsd_cfg,
    )

    if not payloads:
        print(
            "No events found in batch manifest."
        )
        return 0

    if a.jobs < 0:
        raise ValueError(
            "--jobs must be >= 0"
        )

    if a.jobs == 0:
        jobs = os.cpu_count() or 1
    else:
        jobs = a.jobs

    jobs = max(
        1,
        min(
            jobs,
            len(payloads),
        ),
    )

    failures = 0

    if jobs == 1:
        results = [
            _worker_process_event(
                payload
            )
            for payload in payloads
        ]

    else:
        results = []

        with ProcessPoolExecutor(
            max_workers=jobs
        ) as executor:
            future_map = {
                executor.submit(
                    _worker_process_event,
                    payload,
                ): payload
                for payload in payloads
            }

            for future in as_completed(
                future_map
            ):
                results.append(
                    future.result()
                )

    for (
        event_id,
        classification,
        output,
        status,
    ) in sorted(results):
        if status != 0:
            failures += 1
            print(
                f"{event_id}: "
                f"ERROR: {classification}",
                file=sys.stderr,
            )
        else:
            print(
                f"{event_id}: "
                f"TSD={classification}; "
                f"output={output}"
            )

    print(
        f"Batch complete: "
        f"events={len(results)}; "
        f"jobs={jobs}; "
        f"failures={failures}"
    )

    return (
        2
        if failures
        else 0
    )


def main(
    argv: Optional[List[str]] = None,
) -> int:
    a = build_parser().parse_args(
        argv
    )

    try:
        if a.batch_manifest:
            return run_batch(a)

        return run_single(a)

    except Exception as exc:
        print(
            f"ERROR: {exc}",
            file=sys.stderr,
        )
        return 2


if __name__ == "__main__":
    raise SystemExit(
        main()
    )
