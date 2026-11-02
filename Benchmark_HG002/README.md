## Benchmarking and biological characterization — RetroAtlas v9

RetroAtlas v9 was evaluated on the HG002 (NA24385) human genome (GRCh38) using published retrotransposon insertion and structural variation callsets from **xTea** [1] and the **Human Genome Structural Variation Consortium Phase 3 (HGSVC3)** [2].

### Concordance with external reference callsets

| Reference callset | Reference insertions | Matched by RetroAtlas v9 | Positional recovery | Family-concordant recovery |
|:---|---:|---:|---:|---:|
| xTea — ALL insertions [1] | 1,642 | 1,592 | **96.95%** | **96.59%** (1,586/1,642) |
| xTea — TPRT subset [1] | 1,531 | 1,498 | **97.84%** | **97.58%** (1,494/1,531) |
| HGSVC3 — HG002 insertions [2] | 1,460 | 1,426 | **97.67%** | **97.53%** (1,424/1,460) |

**Matching protocol:** Calls were compared on harmonized GRCh38 coordinates using chromosome-consistent, one-to-one positional matching within **±150 bp**. *Positional recovery* is the percentage of reference insertions with a matching RetroAtlas call. *Family-concordant recovery* additionally requires agreement in the assigned retroelement family. The xTea ALL, xTea TPRT, and HGSVC3 callsets overlap and must **not** be treated as independent or added together.

> **Interpretation:** These values quantify overlap with external callsets, not experimentally established sensitivity, specificity, or positive predictive value. The reference catalogs are not exhaustive, and a call absent from them is not automatically a false positive.

### Additional retroelement-associated insertion characterization

| Metric | RetroAtlas v9 (HG002) |
|:---|---:|
| **Total reported insertion loci** | **2,769** |
| Alu-associated insertions | 2,040 |
| LINE-1-associated insertions | 305 |
| SVA-associated insertions | 414 |
| HERVK-associated insertions | 6 |
| Multi-family (`multi_active`) calls | 4 |
| Putative multicopy events (`multi_same`) | 40 |
| Putative LINE-1 inversions (`L1_inversion`) | 35 |
| Embedded retroelement architectures (`embedded`) | 487 |
| Loci with both poly(A/T) and target-site duplication evidence (`BIO=TRUE`) | 1,283 |

**Classification note:** Family counts refer to the principal assigned class. Architecture counts are event labels and may span more than one family; therefore, family and architecture rows are **not additive**. Architectural labels describe computationally inferred structures and should not be interpreted as experimentally confirmed mechanisms.

### Beyond positional detection

RetroAtlas integrates insertion detection with downstream, sequence-aware annotation and preserves evidence for later reinterpretation:

- **Retroelement family assignment:** identification of Alu, LINE-1, SVA, and HERVK-associated sequence, including partial and fragmented RepeatMasker annotations.
- **Multifragment annotation integration:** aggregation of nonredundant repeat evidence across multiple alignments within an inserted sequence.
- **Architecture inference:** annotation of putative single-component, multicopy, inverted, and embedded events where the evidence supports such interpretations.
- **Molecular signature reporting:** poly(A/T) and target-site duplication (TSD) evidence, recorded separately from family identification and architecture.
- **Genotype and sequence preservation:** retention of insertion genotype and available inserted allele sequences, enabling later reannotation and refinement.
- **Conservative handling of uncertainty:** plausible retroelement-associated insertions can be retained even when a specific architecture or retrotransposition mechanism cannot be established.

The presence of retroelement-derived sequence does **not**, by itself, establish that a newly detected insertion arose through an active retrotransposition mechanism. Likewise, absence of a detectable poly(A) tail or TSD does not automatically exclude a complex or truncated retroelement-associated event.

### Scope and limitations

These numbers refer to **RetroAtlas v9 on HG002** and represent a benchmark of the current experimental classifier, not a claim of superiority over xTea or HGSVC3. The pipeline's complex-architecture assignments remain under review, particularly for embedded and multi-family events. Independent validation and explicit false-positive assessment are required before making publication-level claims about precision or overall performance.

### References

1. **Chu, C., Borges-Monroy, R., Viswanadham, V. V., et al. (2021).** Comprehensive identification of transposable element insertions using multiple sequencing technologies. *Nature Communications*, **12**, 3836. https://doi.org/10.1038/s41467-021-24041-8
2. **Logsdon, G. A., Ebert, P., Audano, P. A., et al. (2025).** Complex genetic variation in nearly complete human genomes. *Nature*, **644**, 430–441. https://doi.org/10.1038/s41586-025-09140-6

