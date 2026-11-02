#!/usr/bin/env Rscript

# ==============================================================================
# TE event classifier for long-read reconstructed insertions (v9 family geometry + short rescue + copy safety)
# ==============================================================================
#
# PURPOSE
# -------
# Convert a RepeatMasker .out file into three deliberately separate descriptors:
#
#   class         : mobile class driving the interpretation
#   event         : insertion architecture
#   completeness  : completeness of the principal reconstructed unit/model
#   repeat_elements: explicit RepeatMasker model(s) for compact/resolved events
#
# Target mobile groups for the population-MEI classifier:
#   Alu, L1, SVA, HERV-K(HML-2)
#
# Ancient/interspersed repeats (MIR, L2, ERVL/MaLR, DNA transposons, etc.) are
# retained as sequence context but do not by themselves create a mobile-event call.
#
# IMPORTANT
# ---------
# - Poly(A), TSD and source-locus/transduction logic are intentionally NOT included.
# - RepeatMasker "*" rows are NOT globally discarded.
# - query + RepeatMasker ID is treated as an initial linkage unit, not as a final
#   biological-copy definition.
# - "single" requires a coherent monotonic single-copy chain through the RM model.
#
# FAMILY-AWARE EVENT-DEFINING EVIDENCE (v7)
# ------------------------------------------------
# Public run_te_classifier() inputs are unchanged. Internal post-processing adds:
#   Alu partial       : >=50 bp coherent reconstructed query sequence
#   SVA partial       : >=100 bp coherent reconstructed query sequence
#   SVA composite     : >=2 effective SVA fragments, each >=100 bp, jointly Q80
#   L1 partial        : >=200 bp coherent L1, no 3' anchoring veto;
#                       3' distance retained as a diagnostic.
#   L1 inversion      : two substantial opposite-orientation model-compatible
#                       arms; 3' anchoring is a diagnostic, NOT a gate.
#   Family rescue     : Q80 non-overlapping union of effective RM intervals,
#                       including up to two Alu blocks, L1 inversion-compatible
#                       orientation and SVA fragments; unresolved topology is
#                       explicitly represented, never called near_full.
#   HERV-K(HML-2)     : architecture-driven; isolated fragments do not become
#                       generic single_partial calls; <50 bp has no class power
# Sub-threshold evidence is retained in RepeatMasker/segment QC but cannot by
# itself define single, mosaic, multi-same, or embedded event identity.
#
# MAIN OUTPUTS
# ------------
# <prefix>_event_classification.csv
# <prefix>_units.csv
# <prefix>_query_metrics.csv
# <prefix>_segments.csv          (optional QC)
#
# DEPENDENCY
# ----------
# install.packages("data.table")
#
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
})


# ==============================================================================
# 0. CONFIGURATION
# ==============================================================================

te_config <- function(
    query_threshold = 80,
    model_threshold = 80,
    group_dominance_threshold = 0.80,
    explicit_not_applicable = TRUE,
    max_repeat_elements = 10L
) {
  list(
    query_threshold = query_threshold,
    model_threshold = model_threshold,
    group_dominance_threshold = group_dominance_threshold,
    explicit_not_applicable = explicit_not_applicable,
    max_repeat_elements = as.integer(max_repeat_elements),

    # Family-aware class-power rules. These are deliberately internal so the
    # Retroposon_aggregate launch contract remains unchanged.
    alu_min_partial_bp = 50L,
    sva_min_partial_bp = 100L,
    sva_composite_min_fragment_bp = 100L,
    sva_composite_min_fragments = 2L,
    l1_3prime_tolerance_bp = 100L, # diagnostic only, not an exclusion criterion
    l1_strong_3prime_tolerance_bp = 25L,
    l1_min_partial_bp = 200L,
    l1_inversion_min_arm_bp = 100L,
    # Preserve previously resolved short-arm inversions if 3' evidence exists.
    l1_inversion_restore_short_anchor = TRUE,
    hml2_min_fragment_bp = 50L,

    # Independent family evidence gate: does NOT infer an architecture.
    # Applied only after all legacy architecture rules fail.
    family_pending_enabled = TRUE,
    family_pending_min_query_bp = 200L,
    family_pending_query_pct = 80,
    family_pending_dominance = 0.80,
    family_pending_min_primary_fragment_bp = 100L,
    family_pending_max_star_only_fraction = 0.15,
    family_pending_groups = c("Alu", "L1", "SVA"),
    # v9: preserve the former short L1 calls without making the L1 3' anchor
  # a global exclusion criterion. Applies ONLY to a short partial unit.
  l1_short_min_coherent_bp = 50L,
  l1_short_restore_3prime = TRUE,
  # Alu 100-199 bp: high-coverage family-preserving rescue (not near_full).
  family_pending_alu_short_min_query_bp = 100L,
  family_pending_alu_short_max_query_bp = 199L,
  family_pending_alu_short_min_primary_fragment_bp = 60L,
  # A possible second Alu copy is not proof of a colinear single element.
  alu_copy_reset_min_piece_bp = 60L,
  alu_copy_reset_min_model_overlap_bp = 30L,
  alu_copy_reset_max_query_overlap_fraction = 0.20,
  alu_multicopy_min_unique_bp = 80L,
  # v8: No artificial requirement for two different L1 model names.
    # A single substantial internal alignment may support the family; 3' is
    # diagnostic only.  Deprecated switch retained for CLI compatibility.
    family_pending_l1_require_multifragment = FALSE,
    family_pending_alu_max_query_bp = 450L,
    family_pending_alu_max_blocks = 2L,
    family_pending_alu_max_gap_bp = 40L,
    family_pending_sva_max_blocks = 5L,
    family_pending_sva_max_gap_bp = 600L,
    family_pending_l1_max_gap_bp = 1000L,
    family_pending_l1_max_orientation_runs = 2L, # audit threshold, not hard veto
    family_pending_l1_min_primary_fraction = 0.35,
    family_pending_l1_min_piece_bp = 100L,
    family_pending_l1_max_pair_overlap = 0.20,
    family_pending_l1_min_star_score = 500,

    hml2_ltr_names = c("LTR5_Hs", "LTR5A", "LTR5B"),
    hml2_internal_name = "HERVK-int",

    interspersed_groups = c(
      "Alu", "L1", "SVA", "MIR", "tRNA", "5S",
      "LTR", "DNA", "Other_SINE", "Other_LINE", "Other_retro"
    ),

    mobile_groups = c("Alu", "L1", "SVA", "HERVK_HML2")
  )
}


# ==============================================================================
# 1. SMALL HELPERS
# ==============================================================================
# 1. SMALL HELPERS
# ==============================================================================

rm_num <- function(x) {
  x <- gsub("[()]", "", x)
  suppressWarnings(as.numeric(x))
}


`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L) y else x
}

named_value <- function(x, key, default = 0) {
  if (is.null(x) || length(x) == 0L || !(key %in% names(x))) {
    return(default)
  }
  unname(x[key])
}


id_key <- function(query, id) {
  paste(query, id, sep = "|||")
}


group_key <- function(query, group) {
  paste(query, group, sep = "|||")
}


add_named_count <- function(x, key, value) {
  if (length(key) == 0L || is.na(key) || !nzchar(key)) {
    return(x)
  }

  if (!(key %in% names(x))) {
    x[key] <- 0
  }

  x[key] <- x[key] + value
  x
}


reduce_intervals <- function(start, end = NULL, merge_adjacent = TRUE) {
  if (is.null(end)) {
    if (is.null(start) || nrow(start) == 0L) {
      return(data.table(start = numeric(), end = numeric()))
    }

    x <- as.data.table(start)
    if (!all(c("start", "end") %in% names(x))) {
      stop("Interval table must contain columns 'start' and 'end'.")
    }

    a <- pmin(x$start, x$end)
    b <- pmax(x$start, x$end)
  } else {
    if (length(start) == 0L) {
      return(data.table(start = numeric(), end = numeric()))
    }

    a <- pmin(start, end)
    b <- pmax(start, end)
  }

  ok <- is.finite(a) & is.finite(b)
  a <- a[ok]
  b <- b[ok]

  if (length(a) == 0L) {
    return(data.table(start = numeric(), end = numeric()))
  }

  x <- data.table(start = a, end = b)
  setorder(x, start, end)

  adjacency <- if (merge_adjacent) 1 else 0

  out_start <- numeric()
  out_end <- numeric()

  cur_s <- x$start[1]
  cur_e <- x$end[1]

  if (nrow(x) > 1L) {
    for (i in 2:nrow(x)) {
      if (x$start[i] <= cur_e + adjacency) {
        cur_e <- max(cur_e, x$end[i])
      } else {
        out_start <- c(out_start, cur_s)
        out_end <- c(out_end, cur_e)
        cur_s <- x$start[i]
        cur_e <- x$end[i]
      }
    }
  }

  out_start <- c(out_start, cur_s)
  out_end <- c(out_end, cur_e)

  data.table(start = out_start, end = out_end)
}


interval_width <- function(iv) {
  if (is.null(iv) || nrow(iv) == 0L) {
    return(0)
  }

  sum(iv$end - iv$start + 1)
}


interval_intersection_width <- function(A, B) {
  A <- reduce_intervals(A)
  B <- reduce_intervals(B)

  if (nrow(A) == 0L || nrow(B) == 0L) {
    return(0)
  }

  i <- 1L
  j <- 1L
  total <- 0

  while (i <= nrow(A) && j <= nrow(B)) {
    lo <- max(A$start[i], B$start[j])
    hi <- min(A$end[i], B$end[j])

    if (lo <= hi) {
      total <- total + hi - lo + 1
    }

    if (A$end[i] < B$end[j]) {
      i <- i + 1L
    } else {
      j <- j + 1L
    }
  }

  total
}


append_interval_env <- function(env, key, start, end) {
  if (exists(key, envir = env, inherits = FALSE)) {
    old <- get(key, envir = env, inherits = FALSE)
    assign(
      key,
      rbind(old, data.table(start = start, end = end)),
      envir = env
    )
  } else {
    assign(
      key,
      data.table(start = start, end = end),
      envir = env
    )
  }

  invisible(NULL)
}


get_interval_env <- function(env, key) {
  if (!exists(key, envir = env, inherits = FALSE)) {
    return(data.table(start = numeric(), end = numeric()))
  }

  get(key, envir = env, inherits = FALSE)
}


reduce_interval_env <- function(env) {
  keys <- ls(env, all.names = TRUE)

  for (k in keys) {
    iv <- get(k, envir = env, inherits = FALSE)
    assign(k, reduce_intervals(iv), envir = env)
  }

  invisible(env)
}


interval_string <- function(iv) {
  if (is.null(iv) || nrow(iv) == 0L) {
    return("")
  }

  paste0(iv$start, "-", iv$end, collapse = ";")
}


# ==============================================================================
# 2. REPEATMAPPER GROUPING
# ==============================================================================

repeat_group_one <- function(class_family) {
  repeat_class <- strsplit(class_family, "/", fixed = TRUE)[[1]][1]

  if (identical(class_family, "Retroposon/SVA")) return("SVA")
  if (identical(class_family, "SINE/Alu"))       return("Alu")
  if (identical(class_family, "SINE/MIR"))       return("MIR")
  if (identical(class_family, "SINE/tRNA"))      return("tRNA")
  if (identical(class_family, "SINE/5S"))        return("5S")
  if (identical(class_family, "LINE/L1"))        return("L1")

  if (identical(repeat_class, "LTR"))            return("LTR")
  if (identical(repeat_class, "DNA"))            return("DNA")
  if (identical(repeat_class, "SINE"))           return("Other_SINE")
  if (identical(repeat_class, "LINE"))           return("Other_LINE")
  if (identical(repeat_class, "Retroposon"))     return("Other_retro")

  if (repeat_class %in% c("Simple_repeat", "Low_complexity")) {
    return("Simple_repeat")
  }

  "Other"
}


mobile_group_one <- function(element_group, class_family, repeat_name, cfg) {
  if (element_group %in% c("Alu", "L1", "SVA")) {
    return(element_group)
  }

  if (
    identical(class_family, "LTR/ERVK") &&
    (
      repeat_name %in% cfg$hml2_ltr_names ||
      identical(repeat_name, cfg$hml2_internal_name)
    )
  ) {
    return("HERVK_HML2")
  }

  NA_character_
}


# ==============================================================================
# 3. READ REPEATMASKER .OUT
# ==============================================================================

read_repeatmasker_out <- function(path, skip = 3L, cfg = te_config()) {
  if (!file.exists(path)) {
    stop("RepeatMasker .out does not exist: ", path)
  }

  lines <- readLines(path, warn = FALSE)

  if (length(lines) <= skip) {
    stop("RepeatMasker file appears empty or malformed.")
  }

  lines <- lines[-seq_len(skip)]
  lines <- trimws(lines)
  lines <- lines[nzchar(lines)]

  tokens <- strsplit(lines, "\\s+")

  # Keep standard .out rows only.
  tokens <- tokens[lengths(tokens) >= 15L]

  if (length(tokens) == 0L) {
    stop("No RepeatMasker rows could be parsed.")
  }

  padded <- lapply(tokens, function(x) {
    if (length(x) < 16L) {
      x <- c(x, rep(NA_character_, 16L - length(x)))
    }
    x[1:16]
  })

  mat <- do.call(rbind, padded)
  rm <- as.data.table(mat)

  setnames(
    rm,
    c(
      "score", "div", "del", "ins",
      "query", "q_start_raw", "q_end_raw", "q_left_raw",
      "strand", "repeat_name", "class_family",
      "rep_pos1_raw", "rep_pos2_raw", "rep_pos3_raw",
      "id_raw", "secondary"
    )
  )

  rm[, row_id := .I]

  rm[, `:=`(
    score = rm_num(score),
    div = suppressWarnings(as.numeric(div)),
    del = suppressWarnings(as.numeric(del)),
    ins = suppressWarnings(as.numeric(ins)),
    q_start = rm_num(q_start_raw),
    q_end = rm_num(q_end_raw),
    q_left = rm_num(q_left_raw),
    rep_pos1 = rm_num(rep_pos1_raw),
    rep_pos2 = rm_num(rep_pos2_raw),
    rep_pos3 = rm_num(rep_pos3_raw),
    id = rm_num(id_raw)
  )]

  rm[, star := !is.na(secondary) & secondary == "*"]

  rm[, `:=`(
    query_start = pmin(q_start, q_end),
    query_end = pmax(q_start, q_end)
  )]

  rm[, query_length := query_end + q_left]

  rm[, `:=`(
    ref_coord_1 = fifelse(strand == "C", rep_pos3, rep_pos1),
    ref_coord_2 = rep_pos2,
    ref_left = fifelse(strand == "C", rep_pos1, rep_pos3)
  )]

  rm[, `:=`(
    ref_start = pmin(ref_coord_1, ref_coord_2),
    ref_end = pmax(ref_coord_1, ref_coord_2)
  )]

  rm[, ref_length := ref_end + ref_left]

  rm[, element_group := vapply(
    class_family,
    repeat_group_one,
    character(1)
  )]

  rm[, mobile_group := mapply(
    mobile_group_one,
    element_group,
    class_family,
    repeat_name,
    MoreArgs = list(cfg = cfg),
    USE.NAMES = FALSE
  )]

  # Effective bp are assigned later at atomic-segment level.
  rm[, effective_bp := 0]

  setorder(rm, query, query_start, query_end, row_id)

  # From this point row_id is deliberately identical to the physical row index.
  # Several chaining functions use it as a direct integer index.
  rm[, row_id := .I]

  rm
}


# ==============================================================================
# 4. ATOMIC SEGMENTATION OF REPEATMASKER EVIDENCE
# ==============================================================================

build_effective_evidence <- function(rm, cfg = te_config(), keep_segments = TRUE) {
  rm <- copy(rm)

  id_potential <- new.env(hash = TRUE, parent = emptyenv())
  id_unique <- new.env(hash = TRUE, parent = emptyenv())
  mobile_intervals <- new.env(hash = TRUE, parent = emptyenv())
  simple_intervals <- new.env(hash = TRUE, parent = emptyenv())
  query_state <- new.env(hash = TRUE, parent = emptyenv())

  seg_acc <- list()
  seg_n <- 0L

  queries <- unique(rm$query)
  # Keep physical row indices intact: chain paths address rm by row_id.
  query_rows <- split(seq_len(nrow(rm)), factor(rm$query, levels = queries))

  for (query_i in seq_along(queries)) {
    q <- queries[[query_i]]
    idx <- query_rows[[query_i]]
    rr <- rm[idx]

    bounds <- sort(unique(c(
      rr$query_start,
      rr$query_end + 1
    )))

    rm_bp <- 0
    interspersed_bp <- 0
    mobile_bp <- 0
    star_only_bp <- 0

    inter_unique <- numeric()
    inter_potential <- numeric()
    mobile_unique <- numeric()
    mobile_potential <- numeric()

    if (length(bounds) >= 2L) {
      for (b in seq_len(length(bounds) - 1L)) {
        s <- bounds[b]
        e <- bounds[b + 1L] - 1
        L <- e - s + 1

        cov_local <- which(
          rr$query_start <= s &
          rr$query_end >= e
        )

        if (length(cov_local) == 0L) {
          next
        }

        cov_idx <- idx[cov_local]

        primary_idx <- cov_idx[!rm$star[cov_idx]]

        # LOCAL rule:
        # non-star evidence wins where it exists;
        # starred evidence is retained where it is the only evidence.
        candidate_idx <- if (length(primary_idx) > 0L) {
          primary_idx
        } else {
          cov_idx
        }

        rm_bp <- rm_bp + L

        if (length(primary_idx) == 0L) {
          star_only_bp <- star_only_bp + L
        }

        candidate_ids <- unique(rm$id[candidate_idx])

        for (iid in candidate_ids) {
          append_interval_env(
            id_potential,
            id_key(q, iid),
            s,
            e
          )
        }

        if (length(candidate_ids) == 1L) {
          append_interval_env(
            id_unique,
            id_key(q, candidate_ids[1]),
            s,
            e
          )
        }

        # Every candidate annotation covering the effective atomic segment receives
        # its segment length as effective evidence, matching the dry-run logic.
        rm[candidate_idx, effective_bp := effective_bp + L]

        # Interspersed-repeat context
        inter_idx <- candidate_idx[
          rm$element_group[candidate_idx] %in% cfg$interspersed_groups
        ]

        if (length(inter_idx) > 0L) {
          interspersed_bp <- interspersed_bp + L

          groups <- unique(rm$element_group[inter_idx])

          for (g in groups) {
            inter_potential <- add_named_count(
              inter_potential,
              g,
              L
            )
          }

          if (length(groups) == 1L) {
            inter_unique <- add_named_count(
              inter_unique,
              groups[1],
              L
            )
          }
        }

        # Population-mobile groups
        mobile_idx <- candidate_idx[
          !is.na(rm$mobile_group[candidate_idx])
        ]

        if (length(mobile_idx) > 0L) {
          mobile_bp <- mobile_bp + L

          groups <- unique(rm$mobile_group[mobile_idx])

          for (g in groups) {
            mobile_potential <- add_named_count(
              mobile_potential,
              g,
              L
            )

            append_interval_env(
              mobile_intervals,
              group_key(q, g),
              s,
              e
            )
          }

          if (length(groups) == 1L) {
            mobile_unique <- add_named_count(
              mobile_unique,
              groups[1],
              L
            )
          }
        }

        if (
          any(
            rm$element_group[candidate_idx] == "Simple_repeat",
            na.rm = TRUE
          )
        ) {
          append_interval_env(
            simple_intervals,
            q,
            s,
            e
          )
        }

        if (keep_segments) {
          seg_n <- seg_n + 1L

          seg_acc[[seg_n]] <- data.table(
            query = q,
            start = s,
            end = e,
            bp = L,
            star_only = length(primary_idx) == 0L,
            candidate_ids = paste(candidate_ids, collapse = ","),
            candidate_repeat_names = paste(
              unique(rm$repeat_name[candidate_idx]),
              collapse = ","
            ),
            candidate_groups = paste(
              unique(rm$element_group[candidate_idx]),
              collapse = ","
            ),
            candidate_mobile_groups = paste(
              unique(na.omit(rm$mobile_group[candidate_idx])),
              collapse = ","
            )
          )
        }
      }
    }

    qlen <- max(rr$query_length, na.rm = TRUE)

    assign(
      q,
      list(
        query_length = qlen,
        repeatmasker_coverage_pct = 100 * rm_bp / qlen,
        interspersed_TE_coverage_pct = 100 * interspersed_bp / qlen,
        mobile_coverage_pct = 100 * mobile_bp / qlen,
        repeatmasker_bp = rm_bp,
        interspersed_bp = interspersed_bp,
        mobile_bp = mobile_bp,
        star_only_bp = star_only_bp,
        inter_unique = inter_unique,
        inter_potential = inter_potential,
        mobile_unique = mobile_unique,
        mobile_potential = mobile_potential
      ),
      envir = query_state
    )
  }

  reduce_interval_env(id_potential)
  reduce_interval_env(id_unique)
  reduce_interval_env(mobile_intervals)
  reduce_interval_env(simple_intervals)

  segments <- if (keep_segments && length(seg_acc) > 0L) {
    rbindlist(seg_acc, use.names = TRUE, fill = TRUE)
  } else {
    data.table()
  }

  list(
    rm = rm,
    query_state = query_state,
    id_potential = id_potential,
    id_unique = id_unique,
    mobile_intervals = mobile_intervals,
    simple_intervals = simple_intervals,
    segments = segments
  )
}


# ==============================================================================
# 5. SINGLE-COPY CHAIN
# ==============================================================================

single_copy_chain <- function(rm, row_indices) {
  if (length(row_indices) == 0L) {
    return(integer())
  }

  idx <- row_indices[
    rm$effective_bp[row_indices] > 0
  ]

  if (length(idx) == 0L) {
    return(integer())
  }

  rr <- rm[idx]

  # A generic simple unit must use one RM repeat model and one orientation.
  if (
    uniqueN(rr$repeat_name) != 1L ||
    uniqueN(rr$strand) != 1L
  ) {
    return(integer())
  }

  strand <- rr$strand[1]

  setorder(rr, query_start, query_end, row_id)

  if (strand == "+") {
    oriented_start <- rr$ref_start
    oriented_end <- rr$ref_end
  } else if (strand == "C") {
    oriented_start <- -rr$ref_end
    oriented_end <- -rr$ref_start
  } else {
    return(integer())
  }

  a <- data.table(
    row_id = rr$row_id,
    query_start = rr$query_start,
    query_end = rr$query_end,
    oriented_start = oriented_start,
    oriented_end = oriented_end,
    effective_bp = rr$effective_bp
  )

  n <- nrow(a)

  if (n == 0L) {
    return(integer())
  }

  dp <- a$effective_bp
  previous <- rep(NA_integer_, n)

  if (n > 1L) {
    for (i in 2:n) {
      for (j in 1:(i - 1L)) {
        # Coherent forward progress through query and model.
        if (
          a$query_end[j] <= a$query_start[i] &&
          a$oriented_end[j] <= a$oriented_start[i]
        ) {
          candidate_score <- dp[j] + a$effective_bp[i]

          if (candidate_score > dp[i]) {
            dp[i] <- candidate_score
            previous[i] <- j
          }
        }
      }
    }
  }

  k <- which.max(dp)
  path <- integer()

  while (!is.na(k)) {
    path <- c(a$row_id[k], path)
    k <- previous[k]
  }

  path
}


# ==============================================================================
# 6. RECONSTRUCT QUERY + ID UNITS
# ==============================================================================

reconstruct_units <- function(evidence, cfg = te_config()) {
  rm <- evidence$rm

  split_idx <- split(
    seq_len(nrow(rm)),
    id_key(rm$query, rm$id)
  )

  out <- vector("list", length(split_idx))
  n_out <- 0L

  for (idx in split_idx) {
    allr <- rm[idx]

    q <- allr$query[1]
    iid <- allr$id[1]

    path <- single_copy_chain(rm, idx)

    qpath <- if (length(path) > 0L) {
      reduce_intervals(
        rm$query_start[path],
        rm$query_end[path]
      )
    } else {
      data.table(start = numeric(), end = numeric())
    }

    unique_iv <- get_interval_env(
      evidence$id_unique,
      id_key(q, iid)
    )

    potential_iv <- get_interval_env(
      evidence$id_potential,
      id_key(q, iid)
    )

    unique_chain_bp <- if (length(path) > 0L) {
      interval_intersection_width(qpath, unique_iv)
    } else {
      0
    }

    potential_chain_bp <- if (length(path) > 0L) {
      interval_intersection_width(qpath, potential_iv)
    } else {
      0
    }

    ref_iv <- if (length(path) > 0L) {
      reduce_intervals(
        rm$ref_start[path],
        rm$ref_end[path]
      )
    } else {
      data.table(start = numeric(), end = numeric())
    }

    ref_length <- if (length(path) > 0L) {
      max(rm$ref_length[path], na.rm = TRUE)
    } else {
      max(allr$ref_length, na.rm = TRUE)
    }

    model_coverage_pct <- if (
      length(path) > 0L &&
      is.finite(ref_length) &&
      ref_length > 0
    ) {
      100 * interval_width(ref_iv) / ref_length
    } else {
      0
    }

    state <- get(q, envir = evidence$query_state, inherits = FALSE)
    qlen <- state$query_length

    mobile_groups <- unique(na.omit(allr$mobile_group))

    mobile <- if (length(mobile_groups) == 1L) {
      mobile_groups[1]
    } else {
      NA_character_
    }

    repeat_name <- if (uniqueN(allr$repeat_name) == 1L) {
      allr$repeat_name[1]
    } else {
      "mixed"
    }

    strand <- if (uniqueN(allr$strand) == 1L) {
      allr$strand[1]
    } else {
      "mixed"
    }

    n_out <- n_out + 1L

    out[[n_out]] <- data.table(
      query = q,
      id = iid,
      repeat_name = repeat_name,
      strand = strand,
      mobile_group = mobile,

      # Unique and potential query support intersected with the monotonic chain.
      query_coverage_unique_pct =
        100 * unique_chain_bp / qlen,

      query_coverage_potential_pct =
        100 * potential_chain_bp / qlen,

      model_coverage_pct = model_coverage_pct,

      # Absolute coherent evidence used by family-aware post-processing.
      query_chain_bp = interval_width(qpath),
      query_support_unique_bp = unique_chain_bp,
      query_support_potential_bp = potential_chain_bp,
      model_chain_bp = interval_width(ref_iv),
      ref_length = ref_length,
      model_3prime_distance_bp = if (
        length(path) > 0L &&
        is.finite(ref_length) &&
        ref_length > 0 &&
        nrow(ref_iv) > 0L
      ) {
        max(0, ref_length - max(ref_iv$end))
      } else {
        NA_real_
      },

      near_full =
        length(path) > 0L &&
        model_coverage_pct >= cfg$model_threshold,

      n_rm_rows = length(idx),
      n_chain_rows = length(path),

      query_chain = list(qpath),
      model_chain = list(ref_iv),
      chain_row_ids = list(path)
    )
  }

  rbindlist(out, use.names = TRUE, fill = TRUE)
}


# ==============================================================================
# 7. SPECIALIZED ARCHITECTURE DETECTORS
# ==============================================================================

mobile_interval_for <- function(evidence, query, mobile_group) {
  get_interval_env(
    evidence$mobile_intervals,
    group_key(query, mobile_group)
  )
}


l1_inversion_like <- function(
    query,
    evidence,
    units,
    cfg = te_config()
) {
  rm <- evidence$rm

  fail <- function(coverage_pct = 0) {
    list(
      call = FALSE,
      coverage_pct = coverage_pct,
      min_3prime_distance_bp = NA_real_,
      anchor_strength = NA_character_
    )
  }

  idx <- which(
    rm$query == query &
    rm$mobile_group == "L1" &
    rm$effective_bp > 0
  )
  if (length(idx) < 2L) return(fail())

  rr <- rm[idx]
  setorder(rr, query_start, query_end, row_id)

  # Both arms must refer to the same RepeatMasker L1 model.
  if (uniqueN(rr$repeat_name) != 1L) return(fail())

  run_id <- cumsum(c(TRUE, rr$strand[-1] != rr$strand[-nrow(rr)]))
  rr[, run_id := run_id]

  if (uniqueN(rr$run_id) != 2L || !setequal(unique(rr$strand), c("+", "C"))) {
    return(fail())
  }

  ref_spans <- vector("list", 2L)
  ref_lengths <- numeric(2L)
  has_short_arm <- FALSE

  for (r in 1:2) {
    run_rows <- rr[run_id == r]
    run_idx <- run_rows$row_id
    path <- single_copy_chain(rm, run_idx)
    if (length(path) == 0L) return(fail())

    path_query_iv <- reduce_intervals(rm$query_start[path], rm$query_end[path])
    if (interval_width(path_query_iv) < cfg$l1_inversion_min_arm_bp)
      has_short_arm <- TRUE
    all_run_iv <- reduce_intervals(rm$query_start[run_idx], rm$query_end[run_idx])
    run_width <- interval_width(all_run_iv)

    if (
      run_width <= 0 ||
      interval_intersection_width(path_query_iv, all_run_iv) / run_width < 0.80
    ) return(fail())

    ref_iv <- reduce_intervals(rm$ref_start[path], rm$ref_end[path])
    ref_spans[[r]] <- c(start = min(ref_iv$start), end = max(ref_iv$end))
    ref_lengths[r] <- max(rm$ref_length[path], na.rm = TRUE)
  }

  state <- get(query, envir = evidence$query_state, inherits = FALSE)
  l1_iv <- mobile_interval_for(evidence, query, "L1")
  coverage_pct <- 100 * interval_width(l1_iv) / state$query_length
  if (coverage_pct < cfg$query_threshold) return(fail(coverage_pct))

  a <- ref_spans[[1]]
  b <- ref_spans[[2]]
  overlap <- max(0, min(a["end"], b["end"]) - max(a["start"], b["start"]) + 1)
  smaller <- min(a["end"] - a["start"] + 1, b["end"] - b["start"] + 1)

  if (smaller > 0 && overlap / smaller >= cfg$group_dominance_threshold) {
    return(fail(coverage_pct))
  }

  # Twin-priming geometry remains strict.  3' anchoring is informative only.
  d3 <- pmax(0, ref_lengths - vapply(ref_spans, function(x) x["end"], numeric(1)))
  min_d3 <- min(d3, na.rm = TRUE)
  # For a short inversion arm, keep the old anchor-supported classification.
  # No-anchor inversions still require two substantial (>=100 bp) arms.
  if (has_short_arm && (!isTRUE(cfg$l1_inversion_restore_short_anchor) ||
      !is.finite(min_d3) || min_d3 > cfg$l1_3prime_tolerance_bp))
    return(fail(coverage_pct))
  list(
    call = TRUE,
    coverage_pct = coverage_pct,
    min_3prime_distance_bp = if (is.finite(min_d3)) min_d3 else NA_real_,
    anchor_strength = if (is.finite(min_d3) && min_d3 <= cfg$l1_strong_3prime_tolerance_bp) {
      "strong"
    } else if (is.finite(min_d3) && min_d3 <= cfg$l1_3prime_tolerance_bp) {
      "borderline"
    } else {
      "not_observed"
    }
  )
}


sva_composite_rescue <- function(
    query,
    evidence,
    cfg = te_config()
) {
  q_local <- query
  state <- get(q_local, envir = evidence$query_state, inherits = FALSE)
  rm <- evidence$rm

  # Only fragments with event-defining SVA power are allowed to build the rescue.
  rr <- rm[
    query == q_local &
    mobile_group == "SVA" &
    effective_bp >= cfg$sva_composite_min_fragment_bp
  ]

  if (nrow(rr) == 0L) {
    return(list(
      call = FALSE, coverage_pct = 0,
      n_qualifying_fragments = 0L, repeat_names = character()
    ))
  }

  iv <- reduce_intervals(rr$query_start, rr$query_end)
  coverage_pct <- 100 * interval_width(iv) / state$query_length
  setorder(rr, query_start, query_end, row_id)

  # SVA subfamily fragments may change names, but a single canonical SVA
  # orientation should be consistent; mixed orientation stays unresolved.
  list(
    call =
      coverage_pct >= cfg$query_threshold &&
      nrow(rr) >= cfg$sva_composite_min_fragments &&
      uniqueN(rr$strand) == 1L,
    coverage_pct = coverage_pct,
    n_qualifying_fragments = nrow(rr),
    repeat_names = collapse_adjacent_names(rr$repeat_name)
  )
}


hervk_architecture <- function(
    query,
    evidence,
    units,
    cfg = te_config()
) {
  rm <- evidence$rm
  q_local <- query

  idx <- which(
    rm$query == q_local &
    rm$mobile_group == "HERVK_HML2" &
    rm$effective_bp > 0
  )

  if (length(idx) == 0L) {
    return(list(
      event_detail = NA_character_, coverage_pct = 0,
      repeat_names = character(), architecture_power = FALSE
    ))
  }

  state <- get(q_local, envir = evidence$query_state, inherits = FALSE)
  uq <- units[units$query == q_local]

  # Solo-LTR: require both Q80 and M80, preserving the established rule.
  solo <- uq[
    repeat_name %in% cfg$hml2_ltr_names &
    query_coverage_unique_pct >= cfg$query_threshold &
    model_coverage_pct >= cfg$model_threshold
  ]

  if (nrow(solo) > 0L) {
    setorder(solo, -query_coverage_unique_pct, -model_coverage_pct)
    best <- solo[1]
    return(list(
      event_detail = "HERVK_solo_LTR_like",
      coverage_pct = best$query_coverage_unique_pct,
      repeat_names = best$repeat_name,
      architecture_power = TRUE
    ))
  }

  # Provirus-like: architecture-driven. Each participating RM fragment must have
  # >=50 effective bp so tiny LTR/internal microhits cannot manufacture a call.
  rr <- rm[idx]
  rr <- rr[effective_bp >= cfg$hml2_min_fragment_bp]

  if (nrow(rr) == 0L) {
    return(list(
      event_detail = NA_character_, coverage_pct = 0,
      repeat_names = character(), architecture_power = FALSE
    ))
  }

  ltr <- rr[repeat_name %in% cfg$hml2_ltr_names]
  internal <- rr[repeat_name == cfg$hml2_internal_name]
  has_ltr <- nrow(ltr) > 0L
  has_internal <- nrow(internal) > 0L
  one_strand <- uniqueN(rr$strand) == 1L
  arch_iv <- reduce_intervals(rr$query_start, rr$query_end)
  coverage_pct <- 100 * interval_width(arch_iv) / state$query_length

  if (has_ltr && has_internal && one_strand && coverage_pct >= cfg$query_threshold) {
    setorder(rr, query_start, query_end, row_id)
    return(list(
      event_detail = "HERVK_provirus_like",
      coverage_pct = coverage_pct,
      repeat_names = collapse_adjacent_names(rr$repeat_name),
      architecture_power = TRUE
    ))
  }

  list(
    event_detail = NA_character_, coverage_pct = coverage_pct,
    repeat_names = character(), architecture_power = FALSE
  )
}


# ==============================================================================
# 8. FINAL QUERY CLASSIFICATION
# ==============================================================================
# 8. FINAL QUERY CLASSIFICATION
# ==============================================================================

coarse_event <- function(event_detail) {
  if (
    event_detail %in% c(
      "single_near_full",
      "single_partial",
      "SVA_composite"
    )
  ) {
    return("single")
  }

  if (event_detail == "L1_5prime_inversion_like") {
    return("L1_inversion")
  }

  if (event_detail == "HERVK_solo_LTR_like") {
    return("HERVK_solo_LTR")
  }

  if (event_detail == "HERVK_provirus_like") {
    return("HERVK_provirus")
  }

  if (event_detail == "multi_same_class") {
    return("multi_same")
  }

  if (event_detail == "multi_class_mosaic") {
    return("mosaic")
  }

  if (event_detail %in% c(
    "embedded_single",
    "embedded_multi"
  )) {
    return("embedded")
  }

  "unresolved"
}


event_completeness <- function(
    event_detail,
    explicit_not_applicable = TRUE
) {
  if (event_detail == "single_near_full") {
    return("near_full")
  }

  if (event_detail == "single_partial") {
    return("partial")
  }

  if (event_detail == "SVA_composite") {
    return("composite_model")
  }

  if (explicit_not_applicable) {
    return("not_applicable")
  }

  NA_character_
}


# Specific RepeatMasker model(s) that define a compact event architecture.
# The value is exported as a single CSV field using "|" as a vector separator.
# Very large architectures are truncated rather than dumping dozens of pieces.
format_repeat_elements <- function(x, max_elements = 10L) {
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(x) & x != "mixed"]

  if (length(x) == 0L) {
    return(NA_character_)
  }

  max_elements <- suppressWarnings(as.integer(max_elements[1]))

  if (!is.finite(max_elements) || max_elements < 1L) {
    max_elements <- 10L
  }

  if (length(x) > max_elements) {
    return(paste0(
      paste(x[seq_len(max_elements)], collapse = "|"),
      "|...(+", length(x) - max_elements, ")"
    ))
  }

  paste(x, collapse = "|")
}


unit_query_start <- function(unit_row) {
  iv <- unit_row$query_chain[[1]]

  if (is.null(iv) || nrow(iv) == 0L) {
    return(Inf)
  }

  min(iv$start)
}


ordered_unit_repeat_names <- function(unit_table) {
  if (is.null(unit_table) || nrow(unit_table) == 0L) {
    return(character())
  }

  starts <- vapply(
    seq_len(nrow(unit_table)),
    function(i) unit_query_start(unit_table[i]),
    numeric(1)
  )

  unit_table$repeat_name[order(starts, unit_table$id)]
}


collapse_adjacent_names <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(x)]

  if (length(x) <= 1L) {
    return(x)
  }

  x[c(TRUE, x[-1] != x[-length(x)])]
}


unit_class_power_one <- function(unit_row, cfg = te_config()) {
  if (is.null(unit_row) || nrow(unit_row) != 1L) return(FALSE)

  mg <- as.character(unit_row$mobile_group[1])
  rn <- as.character(unit_row$repeat_name[1])
  bp <- suppressWarnings(as.numeric(unit_row$query_chain_bp[1]))
  d3 <- suppressWarnings(as.numeric(unit_row$model_3prime_distance_bp[1]))
  nf <- isTRUE(unit_row$near_full[1])

  if (is.na(mg) || !nzchar(mg)) return(FALSE)
  if (!is.finite(bp)) bp <- 0

  if (nf) {
    if (mg == "HERVK_HML2") return(rn %in% cfg$hml2_ltr_names)
    return(TRUE)
  }

  if (mg == "Alu") return(bp >= cfg$alu_min_partial_bp)
  if (mg == "SVA") return(bp >= cfg$sva_min_partial_bp)
  if (mg == "L1") {
    if (bp >= cfg$l1_min_partial_bp) return(TRUE)
    # Retain short L1 partials accepted before v8. A 3' model endpoint
    # is supporting evidence only in this low-information short stratum.
    return(isTRUE(cfg$l1_short_restore_3prime) &&
      bp >= cfg$l1_short_min_coherent_bp &&
      is.finite(d3) && d3 <= cfg$l1_3prime_tolerance_bp)
  }

  # HML-2 requires recognized solo-LTR/proviral architecture.
  FALSE
}


annotate_unit_class_power <- function(units, cfg = te_config()) {
  out <- copy(units)

  if (nrow(out) == 0L) {
    out[, class_power_pass := logical()]
    out[, class_power_reason := character()]
    out[, l1_3prime_anchor_strength := character()]
    return(out)
  }

  out[, class_power_pass := vapply(
    seq_len(.N),
    function(i) unit_class_power_one(out[i], cfg),
    logical(1)
  )]

  out[, l1_3prime_anchor_strength := fifelse(
    mobile_group == "L1" &
      is.finite(model_3prime_distance_bp) &
      model_3prime_distance_bp <= cfg$l1_strong_3prime_tolerance_bp,
    "strong",
    fifelse(
      mobile_group == "L1" &
        is.finite(model_3prime_distance_bp) &
        model_3prime_distance_bp <= cfg$l1_3prime_tolerance_bp,
      "borderline",
      NA_character_
    )
  )]

  out[, class_power_reason := vapply(seq_len(.N), function(i) {
    u <- out[i]
    mg <- as.character(u$mobile_group[1])
    if (is.na(mg) || !nzchar(mg)) return("non_mobile_context")
    if (isTRUE(u$class_power_pass[1])) {
      if (isTRUE(u$near_full[1])) return("near_full_model")
      if (mg == "Alu") return("Alu_partial_ge_50bp")
      if (mg == "SVA") return("SVA_partial_ge_100bp")
      if (mg == "L1") return(if (is.finite(u$query_chain_bp[1]) &&
             u$query_chain_bp[1] < cfg$l1_min_partial_bp)
             "L1_short_partial_restored_with_3prime_support"
             else "L1_coherent_partial_no_3prime_requirement")
      return("family_specific_power")
    }
    if (mg == "Alu") return("Alu_microfragment_lt_50bp")
    if (mg == "SVA") return("SVA_microfragment_lt_100bp")
    if (mg == "L1") return("L1_coherent_fragment_lt_minimum_bp")
    if (mg == "HERVK_HML2") return("HML2_fragment_without_architecture")
    "insufficient_class_power"
  }, character(1))]

  out
}


powered_group_geometry <- function(active_powered, state) {
  if (is.null(active_powered) || nrow(active_powered) == 0L) {
    return(list(
      groups = character(), group_bp = numeric(), total_bp = 0,
      coverage_pct = 0, n_unique_groups = 0L, max_fraction = 0
    ))
  }

  groups <- sort(unique(active_powered$mobile_group))
  group_bp <- setNames(numeric(length(groups)), groups)
  all_iv <- list()

  for (g in groups) {
    gg <- active_powered[mobile_group == g]
    iv <- if (nrow(gg) > 0L) {
      reduce_intervals(rbindlist(gg$query_chain, use.names = TRUE, fill = TRUE))
    } else {
      data.table(start = numeric(), end = numeric())
    }
    group_bp[g] <- interval_width(iv)
    if (nrow(iv) > 0L) all_iv[[length(all_iv) + 1L]] <- iv
  }

  union_all <- if (length(all_iv) > 0L) {
    reduce_intervals(rbindlist(all_iv, use.names = TRUE, fill = TRUE))
  } else {
    data.table(start = numeric(), end = numeric())
  }
  total_bp <- interval_width(union_all)
  unique_powered_groups <- intersect(
    names(state$mobile_unique)[state$mobile_unique > 0], groups
  )

  list(
    groups = groups, group_bp = group_bp, total_bp = total_bp,
    coverage_pct = if (state$query_length > 0) 100 * total_bp / state$query_length else 0,
    n_unique_groups = length(unique_powered_groups),
    max_fraction = if (total_bp > 0 && length(group_bp) > 0L) max(group_bp) / total_bp else 0
  )
}


# For unresolved L1, two distinct substantial parts are safer than accepting
# any internal piece with Q80.  This does not assert a correct L1 model chain.
l1_pending_has_distinct_pieces <- function(rr, cfg) {
  if (nrow(rr) < 2L) return(FALSE)
  eligible <- which(
    rr$effective_bp >= cfg$family_pending_l1_min_piece_bp &
    (!rr$star | rr$score >= cfg$family_pending_l1_min_star_score)
  )
  if (length(eligible) < 2L) return(FALSE)
  for (i in seq_len(length(eligible) - 1L)) {
    aa <- rr[eligible[[i]]]
    for (j in seq.int(i + 1L, length(eligible))) {
      bb <- rr[eligible[[j]]]
      distinct <- !identical(aa$repeat_name[[1L]], bb$repeat_name[[1L]]) ||
        !identical(aa$strand[[1L]], bb$strand[[1L]])
      if (!distinct) next
      overlap <- max(0, min(aa$query_end, bb$query_end) -
        max(aa$query_start, bb$query_start) + 1)
      shorter <- min(aa$query_end - aa$query_start + 1,
                     bb$query_end - bb$query_start + 1)
      if (shorter > 0 && overlap / shorter <= cfg$family_pending_l1_max_pair_overlap) {
        return(TRUE)
      }
    }
  }
  FALSE
}

# v8: make the family-level decision on physical, nonredundant query bases.
# RepeatMasker IDs are annotation/linkage hypotheses, NOT copy identifiers.
# A primary alignment outranks a starred overlapping alternative; the highest
# scoring annotation wins at each atomic query interval.  Orientation blocks
# are then defined by *non-overlapping* physical segments, not raw hit counts.
family_fragment_geometry <- function(rr, qlen) {
  empty <- list(coverage_bp = 0, coverage_pct = 0, blocks = 0L,
                orientation_runs = 0L, max_gap_bp = 0,
                block_bp = numeric(), strands = character(),
                models = character())
  rr <- rr[effective_bp > 0 & is.finite(query_start) & is.finite(query_end) &
           strand %in% c("+", "C")]
  if (!nrow(rr) || !is.finite(qlen) || qlen <= 0) return(empty)

  bounds <- sort(unique(c(rr$query_start, rr$query_end + 1)))
  atoms <- vector("list", max(length(bounds) - 1L, 0L))
  used <- 0L
  for (i in seq_len(length(bounds) - 1L)) {
    s <- bounds[i]
    e <- bounds[i + 1L] - 1
    candidates <- rr[query_start <= s & query_end >= e]
    if (!nrow(candidates)) next
    # A secondary annotation can explain a genuinely unmasked interval but
    # may not displace a primary annotation of the same bases.
    setorder(candidates, star, -score, -effective_bp, row_id)
    winner <- candidates[1L]
    used <- used + 1L
    atoms[[used]] <- data.table(start = s, end = e,
      strand = as.character(winner$strand[1L]),
      repeat_name = as.character(winner$repeat_name[1L]))
  }
  if (!used) return(empty)
  aa <- rbindlist(atoms[seq_len(used)])
  setorder(aa, start, end)
  coverage_bp <- sum(aa$end - aa$start + 1)
  gap <- c(0, pmax(0, aa$start[-1L] - aa$end[-nrow(aa)] - 1))
  switch <- c(FALSE, aa$strand[-1L] != aa$strand[-nrow(aa)])
  # Only a change of orientation OR an uncovered interval starts a new block.
  new_block <- c(TRUE, (gap[-1L] > 0) | switch[-1L])
  aa[, block_id := cumsum(new_block)]
  bb <- aa[, .(bp = sum(end - start + 1), strand = strand[1L]), by = block_id]
  # Adjacent blocks of the same orientation separated by short gaps remain
  # distinct physical fragments but one orientation run.
  list(
    coverage_bp = coverage_bp, coverage_pct = 100 * coverage_bp / qlen,
    blocks = as.integer(nrow(bb)),
    orientation_runs = as.integer(1L + sum(aa$strand[-1L] != aa$strand[-nrow(aa)])),
    max_gap_bp = if (length(gap)) max(gap) else 0,
    block_bp = bb$bp,
    strands = unique(aa$strand),
    models = unique(as.character(aa$repeat_name))
  )
}


# Physical Alu copy-ambiguity guard.  A repeated segment of the SAME Alu
# consensus at independent query positions may be two partial copies, not a
# two-piece representation of one copy. It is an AMBIGUITY flag, not an
# automatic multi_same call: repeated consensus territory alone does not
# establish the number of biological copies.
alu_possible_copy_reset <- function(q, evidence, cfg) {
  rr <- evidence$rm[query == q & mobile_group == "Alu" &
    effective_bp >= cfg$alu_copy_reset_min_piece_bp & star == FALSE]
  if (nrow(rr) < 2L) return(FALSE)
  setorder(rr, query_start, query_end, row_id)
  for (i in seq_len(nrow(rr)-1L)) {
    for (j in seq.int(i+1L, nrow(rr))) {
      aa <- rr[i]; bb <- rr[j]
      if (aa$repeat_name[1L] != bb$repeat_name[1L] ||
          aa$strand[1L] != bb$strand[1L]) next
      q_overlap <- max(0, min(aa$query_end[1L], bb$query_end[1L]) -
        max(aa$query_start[1L], bb$query_start[1L]) + 1)
      q_short <- min(aa$query_end[1L] - aa$query_start[1L] + 1,
                     bb$query_end[1L] - bb$query_start[1L] + 1)
      if (q_short <= 0 || q_overlap/q_short >
          cfg$alu_copy_reset_max_query_overlap_fraction) next
      ref_overlap <- max(0, min(aa$ref_end[1L], bb$ref_end[1L]) -
        max(aa$ref_start[1L], bb$ref_start[1L]) + 1)
      if (is.finite(ref_overlap) &&
          ref_overlap >= cfg$alu_copy_reset_min_model_overlap_bp)
        return(TRUE)
    }
  }
  FALSE
}

# Try to restore a genuinely colinear Alu split across RM *IDs* when the
# underlying consensus model and orientation agree.  Mixed subfamily models
# cannot be projected safely onto one consensus: they remain family-supported,
# architecture-unresolved, rather than an invented single_near_full.
alu_same_model_two_piece <- function(q, evidence, cfg) {
  fail <- list(call = FALSE, coverage_pct = 0, model_pct = NA_real_,
               names = character(), near_full = FALSE)
  rr <- evidence$rm[query == q & mobile_group == "Alu" & effective_bp > 0]
  if (nrow(rr) < 2L || alu_possible_copy_reset(q, evidence, cfg))
    return(fail)
  state <- get(q, envir = evidence$query_state, inherits = FALSE)
  qlen <- state$query_length
  if (qlen < cfg$family_pending_min_query_bp ||
      qlen > cfg$family_pending_alu_max_query_bp) return(fail)
  candidates <- list()
  n <- 0L
  for (name in unique(rr$repeat_name)) {
    for (orientation in c("+", "C")) {
      sub <- rr[repeat_name == name & strand == orientation]
      if (nrow(sub) < 2L) next
      path <- single_copy_chain(evidence$rm, sub$row_id)
      if (length(path) != 2L) next
      path_rr <- evidence$rm[path]
      qiv <- reduce_intervals(path_rr$query_start, path_rr$query_end)
      r_iv <- reduce_intervals(path_rr$ref_start, path_rr$ref_end)
      gap <- max(0, path_rr$query_start[2L] - path_rr$query_end[1L] - 1)
      coverage_pct <- 100 * interval_width(qiv) / qlen
      model_pct <- 100 * interval_width(r_iv) / max(path_rr$ref_length)
      if (gap > cfg$family_pending_alu_max_gap_bp ||
          coverage_pct < cfg$query_threshold ||
          min(path_rr$effective_bp) < cfg$alu_min_partial_bp ||
          sum(path_rr$star) > 0L) next
      n <- n + 1L
      candidates[[n]] <- list(call = TRUE, coverage_pct = coverage_pct,
           model_pct = model_pct, names = name,
           near_full = model_pct >= cfg$model_threshold)
    }
  }
  if (!length(candidates)) return(fail)
  candidates[[which.max(vapply(candidates, `[[`, numeric(1), "coverage_pct"))]]
}

# Separate identity evidence from architecture.  Never sum repeated bases twice,
# never interpret the number of RM matches as independent biological fragments,
# and do not use 3' proximity to veto L1.  >1 L1 orientation switch is
# reported as an ambiguous topology, NOT promoted to a twin-priming structure.
# Family-specific limits restrain Alu/SVA reconstruction and large gaps.
family_pending_gate <- function(q, evidence, cfg) {
  state <- get(q, envir = evidence$query_state, inherits = FALSE)
  qlen <- as.numeric(state$query_length)
  fallback <- list(
    pass = FALSE, family = NA_character_, family_bp = 0,
    family_pct = 0, dominance = 0, n_hits = 0L,
    primary_fragment_bp = 0, star_only_fraction = 0,
    l1_distinct_pieces = FALSE, blocks = 0L,
    orientation_runs = 0L, max_gap_bp = 0,
    models = character(), reason = "NO_MOBILE_FAMILY_EVIDENCE",
    alu_possible_multicopy = FALSE, alu_short_rescue = FALSE
  )
  if (!is.finite(qlen) || qlen <= 0L) return(fallback)

  fam_bp <- state$mobile_unique
  fam_bp <- fam_bp[names(fam_bp) %in% cfg$family_pending_groups]
  fam_bp <- fam_bp[is.finite(fam_bp) & fam_bp > 0]
  if (!length(fam_bp)) return(fallback)
  fam <- names(fam_bp)[which.max(fam_bp)]
  rr <- evidence$rm[query == q & mobile_group == fam & effective_bp > 0]
  profile <- family_fragment_geometry(rr, qlen)
  fbp <- profile$coverage_bp
  fpct <- profile$coverage_pct
  dominance <- as.numeric(fam_bp[[fam]]) / max(as.numeric(state$mobile_bp), 1)
  primary <- rr[star == FALSE]
  primary_bp <- if (nrow(primary)) max(primary$effective_bp, na.rm = TRUE) else 0
  star_frac <- as.numeric(state$star_only_bp) / qlen
  l1_pieces <- if (identical(fam, "L1")) l1_pending_has_distinct_pieces(rr, cfg) else FALSE
  alu_short <- identical(fam, "Alu") &&
    qlen >= cfg$family_pending_alu_short_min_query_bp &&
    qlen <= cfg$family_pending_alu_short_max_query_bp
  min_query <- if (alu_short) cfg$family_pending_alu_short_min_query_bp else
    cfg$family_pending_min_query_bp
  min_primary <- if (alu_short) cfg$family_pending_alu_short_min_primary_fragment_bp else
    cfg$family_pending_min_primary_fragment_bp
  possible_alu_multicopy <- if (identical(fam, "Alu"))
    alu_possible_copy_reset(q, evidence, cfg) else FALSE
  reason <- if (!isTRUE(cfg$family_pending_enabled)) {
    "FAMILY_PENDING_DISABLED"
  } else if (qlen < min_query) {
    "QUERY_TOO_SHORT"
  } else if (fpct < cfg$family_pending_query_pct) {
    "FAMILY_COVERAGE_BELOW_MIN"
  } else if (dominance < cfg$family_pending_dominance) {
    "AMBIGUOUS_MOBILE_FAMILY"
  } else if (primary_bp < min_primary) {
    "NO_SUBSTANTIAL_PRIMARY_ALIGNMENT"
  } else if (identical(fam, "L1") &&
             isTRUE(cfg$family_pending_l1_require_multifragment) && !l1_pieces) {
    "L1_LEGACY_REQUIRE_MULTIFRAGMENT"
  } else if (identical(fam, "L1") &&
             (100 * interval_width(reduce_intervals(primary$query_start, primary$query_end)) / qlen <
              100 * cfg$family_pending_l1_min_primary_fraction)) {
    "L1_INSUFFICIENT_PRIMARY_COVERAGE"
  } else if (!identical(fam, "L1") &&
             star_frac > cfg$family_pending_max_star_only_fraction) {
    "EXCESS_STAR_ONLY_EVIDENCE"
  } else if (identical(fam, "Alu") &&
             qlen > cfg$family_pending_alu_max_query_bp) {
    "ALU_QUERY_TOO_LONG_FOR_SINGLE"
  } else if (identical(fam, "Alu") &&
             profile$blocks > cfg$family_pending_alu_max_blocks) {
    "ALU_MORE_THAN_TWO_PHYSICAL_BLOCKS"
  } else if (identical(fam, "Alu") && profile$orientation_runs > 1L) {
    "ALU_MIXED_ORIENTATION"
  } else if (identical(fam, "Alu") &&
             profile$max_gap_bp > cfg$family_pending_alu_max_gap_bp) {
    "ALU_LARGE_INTERNAL_GAP"
  } else if (identical(fam, "SVA") &&
             profile$orientation_runs > 1L) {
    "SVA_MIXED_ORIENTATION"
  } else if (identical(fam, "SVA") &&
             profile$blocks > cfg$family_pending_sva_max_blocks) {
    "SVA_EXCESS_FRAGMENTATION"
  } else if (identical(fam, "SVA") &&
             profile$max_gap_bp > cfg$family_pending_sva_max_gap_bp) {
    "SVA_LARGE_INTERNAL_GAP"
  } else if (identical(fam, "L1") &&
             profile$max_gap_bp > cfg$family_pending_l1_max_gap_bp) {
    "L1_LARGE_INTERNAL_GAP"
  } else if (possible_alu_multicopy) {
    "FAMILY_SUPPORTED_ALU_POSSIBLE_MULTICOPY"
  } else if (alu_short) {
    "FAMILY_SUPPORTED_ALU_SHORT_ARCHITECTURE_PENDING"
  } else {
    "FAMILY_SUPPORTED_ARCHITECTURE_PENDING"
  }
  list(
    pass = reason %in% c("FAMILY_SUPPORTED_ARCHITECTURE_PENDING",
      "FAMILY_SUPPORTED_ALU_SHORT_ARCHITECTURE_PENDING",
      "FAMILY_SUPPORTED_ALU_POSSIBLE_MULTICOPY"),
    family = fam, family_bp = fbp, family_pct = fpct,
    dominance = dominance, n_hits = as.integer(nrow(rr)),
    primary_fragment_bp = as.numeric(primary_bp),
    star_only_fraction = star_frac,
    l1_distinct_pieces = l1_pieces,
    blocks = profile$blocks,
    orientation_runs = profile$orientation_runs,
    max_gap_bp = profile$max_gap_bp,
    models = profile$models, reason = reason,
    alu_possible_multicopy = possible_alu_multicopy,
    alu_short_rescue = alu_short
  )
}


classify_queries <- function(
    evidence,
    units,
    cfg = te_config()
) {
  queries <- ls(evidence$query_state, all.names = TRUE)
  out <- vector("list", length(queries))
  n_out <- 0L

  if (!("class_power_pass" %in% names(units))) {
    units <- annotate_unit_class_power(units, cfg)
  }
  # Include empty groups for queries without units and preserve within-query order.
  unit_rows <- split(seq_len(nrow(units)), factor(units$query, levels = queries))

  for (query_i in seq_along(queries)) {
    q <- queries[[query_i]]
    state <- get(q, envir = evidence$query_state, inherits = FALSE)
    qlen <- state$query_length
    uq <- units[unit_rows[[query_i]]]

    active <- uq[!is.na(mobile_group) & mobile_group %in% cfg$mobile_groups]
    active_powered <- active[class_power_pass == TRUE]
    near_full <- active_powered[near_full == TRUE & query_coverage_unique_pct > 0]
    dominant_raw <- active[query_coverage_unique_pct >= cfg$query_threshold]
    dominant <- dominant_raw[class_power_pass == TRUE]

    if (nrow(dominant_raw) > 0L) {
      setorder(dominant_raw, -query_coverage_unique_pct, -model_coverage_pct)
    }
    if (nrow(dominant) > 0L) {
      setorder(dominant, -query_coverage_unique_pct, -model_coverage_pct)
    }

    power_geom <- powered_group_geometry(active_powered, state)
    alu_copy_ambiguous <- if ("Alu" %in% names(state$mobile_potential))
      alu_possible_copy_reset(q, evidence, cfg) else FALSE

    event_detail <- "unresolved"
    event_class <- "none"
    event_coverage_pct <- 0
    model_coverage_pct <- NA_real_
    repeat_elements <- character()
    reason <- ""
    resolved <- FALSE

    h <- if ("HERVK_HML2" %in% names(state$mobile_potential)) {
      hervk_architecture(q, evidence, units, cfg)
    } else {
      list(event_detail = NA_character_, coverage_pct = 0, repeat_names = character(), architecture_power = FALSE)
    }

    l1 <- if ("L1" %in% names(state$mobile_potential)) {
      l1_inversion_like(q, evidence, units, cfg)
    } else {
      list(call = FALSE, coverage_pct = 0, min_3prime_distance_bp = NA_real_, anchor_strength = NA_character_)
    }

    sva <- if ("SVA" %in% names(state$mobile_potential)) {
      sva_composite_rescue(q, evidence, cfg)
    } else {
      list(call = FALSE, coverage_pct = 0, n_qualifying_fragments = 0L, repeat_names = character())
    }

    # 1) HML-2 architecture has highest priority.
    if (!is.na(h$event_detail)) {
      event_detail <- h$event_detail
      event_class <- "HERVK_HML2"
      event_coverage_pct <- h$coverage_pct
      repeat_elements <- h$repeat_names
      reason <- "recognized HML-2 architecture with event-defining evidence"
      resolved <- TRUE
    }

    # 2) L1 inversion-compatible geometry.  A 3' anchor is no longer obligatory.
    if (!resolved && isTRUE(l1$call)) {
      event_detail <- "L1_5prime_inversion_like"
      event_class <- "L1"
      event_coverage_pct <- l1$coverage_pct

      l1_rr <- evidence$rm[query == q & mobile_group == "L1" & effective_bp > 0]
      setorder(l1_rr, query_start, query_end, row_id)
      if (nrow(l1_rr) > 0L) {
        l1_run <- cumsum(c(TRUE, l1_rr$strand[-1] != l1_rr$strand[-nrow(l1_rr)]))
        repeat_elements <- vapply(
          split(seq_len(nrow(l1_rr)), l1_run),
          function(ii) l1_rr$repeat_name[ii[1]],
          character(1)
        )
      }

      reason <- paste0(
        "two coherent opposite-orientation L1 arms with independent model territory; ",
        "3-prime diagnostic=", l1$anchor_strength, " (distance=",
        ifelse(is.finite(l1$min_3prime_distance_bp), round(l1$min_3prime_distance_bp), "NA"), " bp)"
      )
      resolved <- TRUE
    }

    # 2b) High-confidence independent copies BEFORE a generic single rescue.
    # Copy-number evidence must take priority over Q80 coverage of one path.
    if (!resolved) {
      multi_same <- NULL

      if (nrow(near_full) >= 2L) {
        by_group <- split(seq_len(nrow(near_full)), near_full$mobile_group)
        for (mg in names(by_group)) {
          if (mg != "Alu") next
          gg <- near_full[by_group[[mg]]]
          if (nrow(gg) < 2L) next

          all_iv <- rbindlist(gg$query_chain, use.names = TRUE, fill = TRUE)
          union_iv <- reduce_intervals(all_iv)
          independent_units <- 0L

          for (i in seq_len(nrow(gg))) {
            this_iv <- gg$query_chain[[i]]
            other_idx <- setdiff(seq_len(nrow(gg)), i)
            other_iv <- if (length(other_idx) > 0L) {
              rbindlist(gg$query_chain[other_idx], use.names = TRUE, fill = TRUE)
            } else {
              data.table(start = numeric(), end = numeric())
            }
            other_iv <- reduce_intervals(other_iv)
            unique_bp <- interval_width(this_iv) - interval_intersection_width(this_iv, other_iv)
            copy_min_bp <- if (mg == "Alu") cfg$alu_multicopy_min_unique_bp else 1L
            if (unique_bp >= copy_min_bp) independent_units <- independent_units + 1L
          }

          group_bp <- if (mg %in% names(power_geom$group_bp)) power_geom$group_bp[[mg]] else 0
          group_coverage_pct <- 100 * group_bp / qlen

          if (group_coverage_pct >= cfg$query_threshold && independent_units >= 2L) {
            multi_same <- list(
              group = mg, coverage_pct = group_coverage_pct,
              union_coverage_pct = 100 * interval_width(union_iv) / qlen,
              repeat_names = ordered_unit_repeat_names(gg)
            )
            break
          }
        }
      }

      if (!is.null(multi_same)) {
        event_detail <- "multi_same_class"
        event_class <- multi_same$group
        event_coverage_pct <- multi_same$coverage_pct
        repeat_elements <- multi_same$repeat_names
        reason <- ">=2 independent near-full powered units of same target class jointly explain >= query threshold"
        resolved <- TRUE
      }
    }

    # 3) Generic single: only family-powered units are eligible.
    if (!resolved && nrow(dominant) > 0L) {
      u <- dominant[1]
      if (u$mobile_group != "HERVK_HML2" &&
          !(u$mobile_group == "Alu" && alu_copy_ambiguous)) {
        event_class <- u$mobile_group
        event_coverage_pct <- u$query_coverage_unique_pct
        model_coverage_pct <- u$model_coverage_pct
        repeat_elements <- u$repeat_name

        if (isTRUE(u$near_full)) {
          event_detail <- "single_near_full"
          reason <- "one coherent unit explains >= query threshold and >= model threshold"
        } else {
          event_detail <- "single_partial"
          reason <- paste0(
            "one coherent unit explains >= query threshold, is < model threshold, ",
            "and passes family class-power rule [", u$class_power_reason, "]"
          )
        }
        resolved <- TRUE
      }
    }

    # 3b) Alu colinear two-hit same-consensus rescue across RM IDs.
    #     Different Alu subfamily consensuses are not interchangeable and
    #     remain architecture-unresolved through the family pending path.
    if (!resolved && "Alu" %in% names(state$mobile_potential)) {
      alu_pair <- alu_same_model_two_piece(q, evidence, cfg)
      if (isTRUE(alu_pair$call) &&
          as.numeric(named_value(state$mobile_unique, "Alu", 0)) /
              max(as.numeric(state$mobile_bp), 1) >= cfg$group_dominance_threshold) {
        event_detail <- if (alu_pair$near_full) "single_near_full" else "single_partial"
        event_class <- "Alu"
        event_coverage_pct <- alu_pair$coverage_pct
        model_coverage_pct <- alu_pair$model_pct
        repeat_elements <- alu_pair$names
        reason <- "colinear two-piece Alu reconstructed across RepeatMasker IDs on the same consensus"
        resolved <- TRUE
      }
    }

    # 4) SVA fragmentation rescue using only >=100 bp fragments.
    if (!resolved && isTRUE(sva$call)) {
      event_detail <- "SVA_composite"
      event_class <- "SVA"
      event_coverage_pct <- sva$coverage_pct
      repeat_elements <- sva$repeat_names
      reason <- paste0(
        "SVA rescue: >=", cfg$sva_composite_min_fragments,
        " fragments each >=", cfg$sva_composite_min_fragment_bp,
        " bp jointly explain >= query threshold"
      )
      resolved <- TRUE
    }

    # 5) Retain legacy order for non-Alu multicopy architectures.
    if (!resolved) {
      multi_same <- NULL

      if (nrow(near_full) >= 2L) {
        by_group <- split(seq_len(nrow(near_full)), near_full$mobile_group)
        for (mg in names(by_group)) {
          if (mg == "Alu") next  # Alu was assessed before the single rule.
          gg <- near_full[by_group[[mg]]]
          if (nrow(gg) < 2L) next

          all_iv <- rbindlist(gg$query_chain, use.names = TRUE, fill = TRUE)
          union_iv <- reduce_intervals(all_iv)
          independent_units <- 0L

          for (i in seq_len(nrow(gg))) {
            this_iv <- gg$query_chain[[i]]
            other_idx <- setdiff(seq_len(nrow(gg)), i)
            other_iv <- if (length(other_idx) > 0L) {
              rbindlist(gg$query_chain[other_idx], use.names = TRUE, fill = TRUE)
            } else {
              data.table(start = numeric(), end = numeric())
            }
            other_iv <- reduce_intervals(other_iv)
            unique_bp <- interval_width(this_iv) - interval_intersection_width(this_iv, other_iv)
            if (unique_bp > 0) independent_units <- independent_units + 1L
          }

          group_bp <- if (mg %in% names(power_geom$group_bp)) power_geom$group_bp[[mg]] else 0
          group_coverage_pct <- 100 * group_bp / qlen

          if (group_coverage_pct >= cfg$query_threshold && independent_units >= 2L) {
            multi_same <- list(
              group = mg, coverage_pct = group_coverage_pct,
              union_coverage_pct = 100 * interval_width(union_iv) / qlen,
              repeat_names = ordered_unit_repeat_names(gg)
            )
            break
          }
        }
      }

      if (!is.null(multi_same)) {
        event_detail <- "multi_same_class"
        event_class <- multi_same$group
        event_coverage_pct <- multi_same$coverage_pct
        repeat_elements <- multi_same$repeat_names
        reason <- ">=2 independent near-full powered units of same target class jointly explain >= query threshold"
        resolved <- TRUE
      }
    }

    # 6) Mosaic: sub-threshold microfragments cannot establish a second class.
    if (!resolved &&
        power_geom$coverage_pct >= cfg$query_threshold &&
        power_geom$n_unique_groups >= 2L &&
        power_geom$max_fraction < cfg$group_dominance_threshold &&
        uniqueN(near_full$mobile_group) >= 2L) {
      event_detail <- "multi_class_mosaic"
      event_class <- "multi_active"
      event_coverage_pct <- power_geom$coverage_pct
      mosaic_units <- active_powered[query_coverage_unique_pct > 0]
      repeat_elements <- ordered_unit_repeat_names(mosaic_units)
      reason <- "target-mobile mosaic using only family-powered evidence, with >=1 near-full unit and no dominant powered class"
      resolved <- TRUE
    }

    # 7) Embedded remains intentionally permissive for near-full powered units.
    if (!resolved && nrow(near_full) == 1L &&
        near_full$query_coverage_unique_pct[1L] < cfg$query_threshold &&
        !(near_full$mobile_group[1L] == "Alu" && alu_copy_ambiguous)) {
      event_detail <- "embedded_single"
      event_class <- near_full$mobile_group[1]
      event_coverage_pct <- power_geom$coverage_pct
      repeat_elements <- near_full$repeat_name[1]
      reason <- "one near-full event-defining mobile unit embedded in a larger insertion"
      resolved <- TRUE
    } else if (!resolved && nrow(near_full) >= 2L && power_geom$coverage_pct < cfg$query_threshold) {
      classes <- sort(unique(near_full$mobile_group))
      event_detail <- "embedded_multi"
      event_class <- if (length(classes) == 1L) classes[1] else "multi_active"
      event_coverage_pct <- power_geom$coverage_pct
      repeat_elements <- ordered_unit_repeat_names(near_full)
      reason <- "multiple near-full event-defining mobile units embedded in a larger insertion"
      resolved <- TRUE
    }

    # 8) Explicit unresolved outcome when dominant evidence fails class power.
    if (!resolved) {
      blocked <- dominant_raw[class_power_pass == FALSE]

      if (nrow(blocked) > 0L) {
        setorder(blocked, -query_coverage_unique_pct, -model_coverage_pct)
        u <- blocked[1]
        event_class <- "unresolved"
        event_coverage_pct <- u$query_coverage_unique_pct
        model_coverage_pct <- u$model_coverage_pct
        repeat_elements <- u$repeat_name

        if (u$mobile_group == "Alu") {
          event_detail <- "Alu_microfragment"
          reason <- paste0("Alu evidence <", cfg$alu_min_partial_bp, " coherent bp; insufficient class power")
        } else if (u$mobile_group == "SVA") {
          event_detail <- "SVA_microfragment"
          reason <- paste0("SVA evidence <", cfg$sva_min_partial_bp, " coherent bp; insufficient class power")
        } else if (u$mobile_group == "L1") {
          event_detail <- "L1_internal_fragment"
          reason <- paste0("L1 partial below ", cfg$l1_min_partial_bp,
            " bp without short-unit 3-prime support")
        } else if (u$mobile_group == "HERVK_HML2") {
          event_detail <- "HERVK_fragment"
          reason <- "HML-2 fragment lacks solo-LTR or provirus architecture; insufficient event-defining power"
        } else {
          reason <- "dominant mobile evidence lacks family-specific class power"
        }

      } else if (length(power_geom$groups) > 0L && power_geom$total_bp > 0) {
        powered_unique <- state$mobile_unique[names(state$mobile_unique) %in% power_geom$groups]

        if (length(powered_unique) > 0L && max(powered_unique) > 0) {
          best_group <- names(which.max(powered_unique))
          best_lower_bp <- max(powered_unique)
          if (best_lower_bp / power_geom$total_bp >= cfg$group_dominance_threshold) {
            event_class <- best_group
            reason <- "powered target-mobile evidence dominated but topology unresolved"
          } else {
            event_class <- if (length(power_geom$groups) >= 2L) "multi_active" else power_geom$groups[1]
            reason <- "powered target-mobile evidence insufficient for defined architecture"
          }
        } else {
          event_class <- if (length(power_geom$groups) >= 2L) "multi_active" else power_geom$groups[1]
          reason <- "powered target-mobile evidence is ambiguous"
        }
      } else {
        event_class <- if (state$mobile_bp > 0) {
          "unresolved"
        } else if (length(state$inter_potential) > 0L) {
          "other_repeat"
        } else {
          "none"
        }
        reason <- if (state$mobile_bp > 0) {
          "mobile RepeatMasker evidence present but below family-specific class power"
        } else {
          "no target polymorphic-retroelement evidence"
        }
      }
    }

    event <- coarse_event(event_detail)
    completeness <- event_completeness(event_detail, cfg$explicit_not_applicable)

    # The legacy classification is NEVER overwritten here.  Fragment assembly
    # still sees the original unresolved status.  The aggregator may promote a
    # family-supported candidate AFTER compound-fragment reconstruction.
    pending <- family_pending_gate(q, evidence, cfg)
    architecture_unresolved <- identical(event, "unresolved")
    legacy_eligible <- !architecture_unresolved &&
      !event_class %in% c("unresolved", "none", "other_repeat")
    family_pending <- architecture_unresolved && isTRUE(pending$pass)
    retro_delete <- !legacy_eligible && !family_pending
    architecture <- if (architecture_unresolved) "unresolved" else event_detail
    delete_reason <- if (!retro_delete) {
      if (family_pending) "FAMILY_SUPPORTED_PENDING" else "RESOLVED_LEGACY"
    } else {
      pending$reason
    }

    n_out <- n_out + 1L
    out[[n_out]] <- data.table(
      query = q,
      query_length = qlen,
      event_class = event_class,
      event = event,
      completeness = completeness,
      event_coverage_pct = round(event_coverage_pct, 3),
      model_coverage_pct = if (is.na(model_coverage_pct)) NA_real_ else round(model_coverage_pct, 3),
      mobile_coverage_pct = round(state$mobile_coverage_pct, 3),
      interspersed_TE_coverage_pct = round(state$interspersed_TE_coverage_pct, 3),
      repeatmasker_coverage_pct = round(state$repeatmasker_coverage_pct, 3),
      n_mobile_groups = length(state$mobile_potential),
      n_near_full_units = nrow(near_full),
      star_only_bp = state$star_only_bp,
      repeat_elements = format_repeat_elements(repeat_elements, cfg$max_repeat_elements),
      reason = reason,
      event_detail = event_detail,
      ARCHITECTURE = architecture,
      FAMILY_EVIDENCE_CLASS = pending$family,
      FAMILY_EVIDENCE_BP = round(pending$family_bp, 0),
      FAMILY_EVIDENCE_QUERY_PCT = round(pending$family_pct, 3),
      FAMILY_EVIDENCE_DOMINANCE = round(pending$dominance, 4),
      FAMILY_EVIDENCE_HITS = pending$n_hits,
      FAMILY_EVIDENCE_PRIMARY_MAX_BP = round(pending$primary_fragment_bp, 0),
      FAMILY_EVIDENCE_STAR_ONLY_FRACTION = round(pending$star_only_fraction, 4),
      FAMILY_EVIDENCE_L1_DISTINCT_PIECES = pending$l1_distinct_pieces,
      FAMILY_EVIDENCE_BLOCKS = pending$blocks,
      FAMILY_EVIDENCE_ORIENTATION_RUNS = pending$orientation_runs,
      FAMILY_EVIDENCE_MAX_GAP_BP = pending$max_gap_bp,
      FAMILY_EVIDENCE_MODELS = format_repeat_elements(pending$models, cfg$max_repeat_elements),
      FAMILY_EVIDENCE_RULE = pending$reason,
      FAMILY_EVIDENCE_ALU_POSSIBLE_MULTICOPY = pending$alu_possible_multicopy,
      FAMILY_EVIDENCE_ALU_SHORT_RESCUE = pending$alu_short_rescue,
      FAMILY_EVIDENCE_COMPLEX_ORIENTATION =
        isTRUE(pending$family == "L1") &&
        pending$orientation_runs > cfg$family_pending_l1_max_orientation_runs,
      FAMILY_SUPPORTED_PENDING = family_pending,
      DELETE = retro_delete,
      DELETE_REASON = delete_reason
    )
  }

  ans <- rbindlist(out, use.names = TRUE, fill = TRUE)
  setorder(ans, query)
  ans
}


# ==============================================================================
# 9. QUERY METRICS EXPORT
# ==============================================================================
# 9. QUERY METRICS EXPORT
# ==============================================================================

query_metrics_table <- function(evidence) {
  queries <- ls(
    evidence$query_state,
    all.names = TRUE
  )

  out <- lapply(queries, function(q) {
    s <- get(
      q,
      envir = evidence$query_state,
      inherits = FALSE
    )

    data.table(
      query = q,
      query_length = s$query_length,
      repeatmasker_coverage_pct =
        s$repeatmasker_coverage_pct,
      interspersed_TE_coverage_pct =
        s$interspersed_TE_coverage_pct,
      mobile_coverage_pct =
        s$mobile_coverage_pct,
      star_only_bp =
        s$star_only_bp,
      mobile_groups =
        paste(
          names(s$mobile_potential),
          collapse = ","
        ),
      interspersed_groups =
        paste(
          names(s$inter_potential),
          collapse = ","
        )
    )
  })

  ans <- rbindlist(
    out,
    use.names = TRUE,
    fill = TRUE
  )

  setorder(ans, query)

  ans
}


# ==============================================================================
# 10. EXPORT UNITS
# ==============================================================================

units_export_table <- function(units) {
  out <- copy(units)

  out[, query_chain_txt := vapply(
    query_chain,
    interval_string,
    character(1)
  )]

  out[, model_chain_txt := vapply(
    model_chain,
    interval_string,
    character(1)
  )]

  out[, chain_row_ids_txt := vapply(
    chain_row_ids,
    function(x) paste(x, collapse = ","),
    character(1)
  )]

  out[, c(
    "query_chain",
    "model_chain",
    "chain_row_ids"
  ) := NULL]

  out
}


# ==============================================================================
# 11. MAIN RUNNER
# ==============================================================================

run_te_classifier <- function(
    repeatmasker_out,
    output_prefix = "TE",
    output_dir = ".",
    query_threshold = 80,
    model_threshold = 80,
    group_dominance_threshold = 0.80,
    explicit_not_applicable = TRUE,
    max_repeat_elements = 10L,
    write_segments = TRUE,
    config_overrides = list(),
    query_ids = NULL,
    quiet = TRUE
) {
  cfg <- te_config(
    query_threshold = query_threshold,
    model_threshold = model_threshold,
    group_dominance_threshold =
      group_dominance_threshold,
    explicit_not_applicable =
      explicit_not_applicable,
    max_repeat_elements = max_repeat_elements
  )

  unknown <- setdiff(names(config_overrides), names(cfg))
  if (length(unknown)) stop("Unknown classifier settings: ", paste(unknown, collapse=", "))
  cfg <- utils::modifyList(cfg, config_overrides)

  dir.create(
    output_dir,
    showWarnings = FALSE,
    recursive = TRUE
  )

  if (!quiet) message("1/5 Reading RepeatMasker .out ...")

  rm <- read_repeatmasker_out(
    repeatmasker_out,
    cfg = cfg
  )

  if (!is.null(query_ids)) {
    rm <- rm[query %in% query_ids]
    # Chains use row_id as a physical index into rm. Filtering removes rows,
    # so restore that contract before building any evidence or chain paths.
    rm[, row_id := .I]
  }
  if (!nrow(rm)) return(list(calls=data.table::data.table(), config=cfg))

  if (!quiet) message(
    "   rows: ", nrow(rm),
    " | queries: ", uniqueN(rm$query),
    " | starred rows: ", sum(rm$star)
  )

  if (!quiet) message("2/5 Building atomic effective evidence ...")

  evidence <- build_effective_evidence(
    rm,
    cfg = cfg,
    keep_segments = write_segments
  )

  if (!quiet) message("3/5 Reconstructing query+ID units ...")

  units <- reconstruct_units(
    evidence,
    cfg = cfg
  )

  units <- annotate_unit_class_power(
    units,
    cfg = cfg
  )

  if (!quiet) message("4/5 Classifying queries ...")

  calls <- classify_queries(
    evidence,
    units,
    cfg = cfg
  )

  qmetrics <- query_metrics_table(
    evidence
  )

  units_export <- units_export_table(
    units
  )

  if (!quiet) message("5/5 Writing outputs ...")

  calls_path <- file.path(
    output_dir,
    paste0(
      output_prefix,
      "_event_classification.csv"
    )
  )

  units_path <- file.path(
    output_dir,
    paste0(
      output_prefix,
      "_units.csv"
    )
  )

  qmetrics_path <- file.path(
    output_dir,
    paste0(
      output_prefix,
      "_query_metrics.csv"
    )
  )

  fwrite(calls, calls_path)
  fwrite(units_export, units_path)
  fwrite(qmetrics, qmetrics_path)

  if (write_segments) {
    segments_path <- file.path(
      output_dir,
      paste0(
        output_prefix,
        "_segments.csv"
      )
    )

    fwrite(
      evidence$segments,
      segments_path
    )
  }

  if (!quiet) message("")
  if (!quiet) message("EVENT counts:")
  if (!quiet) print(calls[, .N, by = event][order(-N)])

  if (!quiet) message("")
  if (!quiet) message("CLASS x EVENT for queries with mobile evidence:")
  if (!quiet) print(
    calls[
      mobile_coverage_pct > 0,
      .N,
      by = .(event_class, event)
    ][order(event_class, event)]
  )

  if (!quiet) message("")
  if (!quiet) message("Family-aware refinement rules:")
  if (!quiet) message("  Alu partial minimum             : ", cfg$alu_min_partial_bp, " bp")
  if (!quiet) message("  SVA partial minimum             : ", cfg$sva_min_partial_bp, " bp")
  if (!quiet) message("  SVA composite fragment minimum : ", cfg$sva_composite_min_fragment_bp, " bp")
  if (!quiet) message("  L1 coherent partial minimum    : ", cfg$l1_min_partial_bp, " bp")
  if (!quiet) message("  L1 3-prime QC threshold        : ", cfg$l1_3prime_tolerance_bp, " bp (diagnostic only)")
  if (!quiet) message("  L1 strong 3-prime QC           : ", cfg$l1_strong_3prime_tolerance_bp, " bp (diagnostic only)")
  if (!quiet) message("  HML-2 fragment evidence floor  : ", cfg$hml2_min_fragment_bp, " bp")

  invisible(
    list(
      calls = calls,
      units = units,
      query_metrics = qmetrics,
      segments = evidence$segments,
      repeatmasker = evidence$rm,
      config = cfg
    )
  )
}


# ==============================================================================
# 12. OPTIONAL SANITY CHECKS
# ==============================================================================

check_known_examples <- function(calls) {
  examples <- c(
    "Atlas_10004;1",
    "Atlas_940;1",
    "Atlas_10716;1",
    "Atlas_1020;1",
    "Atlas_935;1",
    "Atlas_9400;1",
    "Atlas_7121;1",
    "Atlas_9981;1",
    "Atlas_320;1",
    "Atlas_10825;1"
  )

  calls[
    query %in% examples,
    .(
      query,
      event_class,
      event,
      completeness,
      event_coverage_pct,
      model_coverage_pct,
      mobile_coverage_pct,
      n_near_full_units,
      star_only_bp,
      repeat_elements,
      reason
    )
  ][match(examples, query)]
}


run_sensitivity <- function(
    repeatmasker_out,
    output_dir = ".",
    thresholds = c(70, 80, 90),
    write_segments = FALSE
) {
  out <- list()

  for (thr in thresholds) {
    message(
      "\n=============================="
    )
    message(
      "Query threshold = ", thr
    )
    message(
      "=============================="
    )

    out[[as.character(thr)]] <- run_te_classifier(
      repeatmasker_out = repeatmasker_out,
      output_prefix = paste0("TE_q", thr),
      output_dir = output_dir,
      query_threshold = thr,
      model_threshold = 80,
      group_dominance_threshold = 0.80,
      explicit_not_applicable = TRUE,
      write_segments = write_segments
    )
  }

  invisible(out)
}


# ==============================================================================
# 13. EXAMPLE USAGE
# ==============================================================================

# Edit the path:
#
# rm_out <- "data/alignments_sample.fasta.out"
#
# result <- run_te_classifier(
#   repeatmasker_out = rm_out,
#   output_prefix = "sample_0039",
#   output_dir = "results/TE_classifier",
#   query_threshold = 80,
#   model_threshold = 80,
#   group_dominance_threshold = 0.80,
#   explicit_not_applicable = TRUE,
#   max_repeat_elements = 10,
#   write_segments = TRUE
# )
#
# # Main per-insertion calls:
# result$calls
#
# # Reconstructed units:
# result$units
#
# # Inspect known examples:
# check_known_examples(result$calls)
#
# # Distribution:
# result$calls[, .N, by = .(event_class, event, completeness)][order(-N)]
#
# # Only target mobile evidence:
# result$calls[
#   mobile_coverage_pct > 0
# ]
#
# # Sensitivity 70/80/90:
# sensitivity <- run_sensitivity(
#   repeatmasker_out = rm_out,
#   output_dir = "results/TE_classifier/sensitivity",
#   thresholds = c(70, 80, 90)
# )
#
# ==============================================================================
# END
# ==============================================================================
