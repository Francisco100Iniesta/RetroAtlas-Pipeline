#!/usr/bin/env Rscript

library(tidyverse)
library(data.table)
data.table::setDTthreads(1L)
library(glue)

# ==============================================================================
# Retroposon aggregate v5 — evidence-aware intra-sample reconstruction
# ==============================================================================
#
# The legacy locus-aware merge is preserved.  The new layer audits, and
# optionally reconstructs, super-fragmented INS signals before the existing
# intra-sample clustering.  SimpleRepeats is used as CONTEXT/RISK only; it never
# forces a merge.  Optional cuteSV/Severus/Sniffles2 VCFs can corroborate the
# size of a proposed compound allele but also never redefine locus identity.
#
# Required CLI:
#   --predata <HitSV pre-data.rds>
#   --repeatmasker <RepeatMasker .out>
#   --sample <sample prefix>
#   --regions <coverage regions>
#   --out <output directory>
#
# Recommended optional CLI:
#   --simple-repeats-bed <SimpleRepeats.bed>
#   --vcf-cutesv <*.vcf[.gz]>
#   --vcf-severus <*.vcf[.gz]>
#   --vcf-sniffles2 <*.vcf[.gz]>
#   --fragmentation-mode audit|highconf|off
#   --fragment-population-policy evidence_aware|annotate|exclude_shadow_complex|exclude_fragmented_suspect
#
# highconf reconstructs resolved compounds plus a conservative ANCHORED_DROPOUT path.
# Source-level complexity is caution only; it never poisons a clean subgroup.
# Default output is one *_retroatlas_intra_final.rds bundle (loci + fragment evidence).
# ==============================================================================

# 1. CLI / parameters -----------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

arg_value <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) return(default)
  if (i == length(args)) stop("Missing value after ", flag)
  args[[i + 1L]]
}

arg_num <- function(flag, default) {
  x <- arg_value(flag, NULL)
  if (is.null(x)) return(as.numeric(default))
  y <- suppressWarnings(as.numeric(x))
  if (is.na(y)) stop("Invalid numeric value for ", flag, ": ", x)
  y
}

arg_int <- function(flag, default) {
  y <- arg_num(flag, default)
  if (y != as.integer(y)) stop("Expected integer for ", flag, ": ", y)
  as.integer(y)
}

arg_bool <- function(flag, default = FALSE) {
  x <- arg_value(flag, NULL)
  if (is.null(x)) return(isTRUE(default))
  z <- tolower(trimws(as.character(x)))
  if (z %in% c("true", "t", "1", "yes", "y")) return(TRUE)
  if (z %in% c("false", "f", "0", "no", "n")) return(FALSE)
  stop("Invalid boolean value for ", flag, ": ", x)
}

# Pipeline interface: coverage is derived exclusively from regions BED.
predata <- arg_value("--predata")
RM_data <- arg_value("--repeatmasker", arg_value("--rm-data"))
sample <- arg_value("--sample")
Regions_data <- arg_value("--regions")
out <- arg_value("--out")
required_cli <- list(
  predata = predata,
  repeatmasker = RM_data,
  sample = sample,
  regions = Regions_data,
  out = out
)

missing_cli <- names(required_cli)[vapply(
  required_cli,
  function(x) is.null(x) || length(x) == 0L || is.na(x) || !nzchar(as.character(x)),
  logical(1)
)]

if (length(missing_cli)) {
  stop(
    "Missing required arguments: ", paste(missing_cli, collapse = ", "), "\nUse:\n",
    "  Rscript Retroposon_aggregate_v4_fragmentation.R \\\n",
    "    --predata PRE.rds --repeatmasker RM.out --sample SAMPLE \\\n",
    "    --regions REGIONS --out OUTDIR\n"
  )
}

for (pp in c(predata, RM_data, Regions_data)) {
  if (!file.exists(pp)) stop("Required input not found: ", pp)
}

dir.create(out, recursive = TRUE, showWarnings = FALSE)

classifier_script <- arg_value("--classifier-script", "te_event_classifier_refined_v7.R")
if (!file.exists(classifier_script)) stop("Classifier script not found: ", classifier_script)
simple_repeats_bed <- arg_value("--simple-repeats-bed", NA_character_)
vcf_cutesv <- arg_value("--vcf-cutesv", NA_character_)
vcf_severus <- arg_value("--vcf-severus", NA_character_)
vcf_sniffles2 <- arg_value("--vcf-sniffles2", NA_character_)
write_audits <- arg_bool("--write-audits", FALSE)
write_legacy_rds <- arg_bool("--write-legacy-rds", FALSE)

CFG <- list(
  # RepeatMasker classifier.
  query_threshold = arg_num("--query-threshold", 80),
  model_threshold = arg_num("--model-threshold", 80),
  group_dominance_threshold = arg_num("--group-dominance-threshold", 0.80),

  primary_min_support = arg_num("--primary-min-support", 3),
  min_svlen = arg_num("--min-svlen", 30),

  # Coverage gates.
  max_depth_ratio = arg_num("--max-depth-ratio", 3.0),
  min_depth = arg_num("--min-depth", 0),

  # Existing intra-sample merge v3 geometry.
  audit_dpos_max = arg_int("--audit-dpos-max", 500),
  exact_locus_dpos = arg_int("--exact-locus-dpos", 0),
  heteroallelic_near_dpos_max = arg_int("--heteroallelic-near-dpos-max", 5),
  very_near_dpos_max = arg_int("--very-near-dpos-max", 2),
  very_near_size_min = arg_num("--very-near-size-min", 0.85),
  near_dpos_max = arg_int("--near-dpos-max", 10),
  near_size_min = arg_num("--near-size-min", 0.90),
  mid_dpos_max = arg_int("--mid-dpos-max", 25),
  mid_size_min = arg_num("--mid-size-min", 0.93),
  far_dpos_max = arg_int("--far-dpos-max", 100),
  far_size_min = arg_num("--far-size-min", 0.95),
  exact_event_flex_size_min = arg_num("--exact-event-flex-size-min", 0.93),
  very_near_event_flex_dpos_max = arg_int("--very-near-event-flex-dpos-max", 5),
  very_near_event_flex_size_min = arg_num("--very-near-event-flex-size-min", 0.97),
  near_event_flex_dpos_max = arg_int("--near-event-flex-dpos-max", 10),
  near_event_flex_size_min = arg_num("--near-event-flex-size-min", 0.98),
  class_conflict_audit_dpos_max = arg_int("--class-conflict-audit-dpos-max", 2),
  class_conflict_audit_size_min = arg_num("--class-conflict-audit-size-min", 0.98),
  svlen_stable_min = arg_num("--svlen-stable-min", 0.95),
  svlen_moderate_min = arg_num("--svlen-moderate-min", 0.85),

  # Terminal poly(A/T) detector.
  polyat_primary_window = arg_int("--polyat-primary-window", 50),
  polyat_primary_exact_run = arg_int("--polyat-primary-exact-run", 9),
  polyat_primary_tolerant_run = arg_int("--polyat-primary-tolerant-run", 15),
  polyat_primary_mismatch_tolerance = arg_int("--polyat-primary-mismatch-tolerance", 2),
  polyat_rescue_window = arg_int("--polyat-rescue-window", 100),
  polyat_rescue_exact_run = arg_int("--polyat-rescue-exact-run", 10),
  polyat_rescue_tolerant_run = arg_int("--polyat-rescue-tolerant-run", 15),
  polyat_rescue_mismatch_tolerance = arg_int("--polyat-rescue-mismatch-tolerance", 3)
)

FRAG <- list(
  mode = tolower(arg_value("--fragmentation-mode", "highconf")),
  population_policy = tolower(arg_value("--fragment-population-policy", "evidence_aware")),

  # Local fragmentation signal.
  local_gap_bp = arg_int("--fragment-local-gap-bp", 150),
  audit_window_bp = arg_int("--fragment-audit-window-bp", 1000),
  auto_max_components = arg_int("--fragment-auto-max-components", 3),
  auto_max_span_bp = arg_int("--fragment-auto-max-span-bp", 500),
  complex_min_components = arg_int("--fragment-complex-min-components", 4),
  complex_span_bp = arg_int("--fragment-complex-span-bp", 1000),
  source_complex_min_components = arg_int("--fragment-source-complex-min-components", 6),
  source_complex_span_bp = arg_int("--fragment-source-complex-span-bp", 2000),
  shadow_max_fraction = arg_num("--fragment-shadow-max-fraction", 0.50),

  # Context only: repeat overlap raises risk but never forces collapse.
  repeat_near_bp = arg_int("--fragment-repeat-near-bp", 0),

  # Optional independent VCF corroboration of COMPOUND_SVLEN.
  vcf_window_bp = arg_int("--fragment-vcf-window-bp", 1000),
  vcf_min_svlen = arg_num("--fragment-vcf-min-svlen", 20),
  vcf_min_similarity = arg_num("--fragment-vcf-min-similarity", 0.80),
  vcf_strong_similarity = arg_num("--fragment-vcf-strong-similarity", 0.90),
  auto_vcf_min_callers = arg_int("--fragment-auto-vcf-min-callers", 2),
  auto_require_vcf = arg_bool("--fragment-auto-require-vcf", FALSE),
  vcf_require_pass = arg_bool("--fragment-vcf-require-pass", FALSE),

  # Conservative rescue of classifier-dropout compounds. Unresolved fragments
  # never define class/event; they may only contribute size around a resolved
  # anchor and require independent VCF corroboration.
  anchored_dropout_enabled = arg_bool("--fragment-anchored-dropout", TRUE),
  anchored_min_component_lr = arg_num("--fragment-anchored-min-component-lr", 5),
  anchored_min_anchor_fraction = arg_num("--fragment-anchored-min-anchor-fraction", 0.50),
  anchored_max_unresolved = arg_int("--fragment-anchored-max-unresolved", 2),
  anchored_vcf_min_callers = arg_int("--fragment-anchored-vcf-min-callers", 2),
  anchored_vcf_similarity = arg_num("--fragment-anchored-vcf-similarity", 0.90),
  anchored_require_same_event = arg_bool("--fragment-anchored-require-same-event", TRUE)
)

if (!FRAG$mode %in% c("off", "audit", "highconf")) {
  stop("--fragmentation-mode must be off, audit, or highconf")
}
if (!FRAG$population_policy %in% c(
  "annotate",
  "evidence_aware",
  "exclude_shadow_complex",
  "exclude_fragmented_suspect"
)) {
  stop(
    "--fragment-population-policy must be annotate, evidence_aware, ",
    "exclude_shadow_complex, or exclude_fragmented_suspect"
  )
}

message("Retroposon aggregate v4 fragmentation-aware")
message("Sample: ", sample)
message("Fragmentation mode: ", FRAG$mode)
message("Population policy: ", FRAG$population_policy)

# 2. Preprocesamiento de datos --------------------------------------
detect_terminal_polyAT <- function(
    seq,

    # ----------------------------------------------------------
    # Argumentos legacy:
    # se mantienen para que las llamadas antiguas no fallen,
    # pero YA NO controlan la detección.
    # ----------------------------------------------------------
    fraction = NULL,
    min_run = NULL,
    mismatch_tolerance = NULL,

    bases = c("A", "T"),

    # ----------------------------------------------------------
    # Parámetros calibrados en HG002 TPRT
    # ----------------------------------------------------------

    # Ventana primaria
    primary_window = 50L,
    primary_exact_run = 9L,
    primary_tolerant_run = 15L,
    primary_mismatch_tolerance = 2L,

    # Ventana rescue
    rescue_window = 100L,
    rescue_exact_run = 10L,
    rescue_tolerant_run = 15L,
    rescue_mismatch_tolerance = 3L
) {

  # ============================================================
  # 1. Normalizar secuencia
  # ============================================================

  seq <- toupper(gsub("\\s+", "", seq))

  if (nchar(seq) == 0) {
    stop("La secuencia está vacía.")
  }

  bases <- unique(toupper(bases))

  if (!all(bases %in% c("A", "C", "G", "T"))) {
    stop("'bases' solo puede contener A, C, G o T.")
  }


  # ============================================================
  # 2. Validación de parámetros
  # ============================================================

  integer_params <- c(
    primary_window,
    primary_exact_run,
    primary_tolerant_run,
    primary_mismatch_tolerance,
    rescue_window,
    rescue_exact_run,
    rescue_tolerant_run,
    rescue_mismatch_tolerance
  )

  if (
    any(is.na(integer_params)) ||
    any(integer_params < 0) ||
    any(integer_params != as.integer(integer_params))
  ) {
    stop("Los parámetros de ventana/run/mismatch deben ser enteros >= 0.")
  }

  if (primary_window < 1L) {
    stop("'primary_window' debe ser >= 1.")
  }

  if (rescue_window < primary_window) {
    stop("'rescue_window' debe ser >= 'primary_window'.")
  }

  if (
    primary_exact_run < 1L ||
    primary_tolerant_run < 1L ||
    rescue_exact_run < 1L ||
    rescue_tolerant_run < 1L
  ) {
    stop("Los tamaños mínimos de run deben ser >= 1.")
  }


  # ============================================================
  # 3. Información de longitud / ventanas
  # ============================================================

  seq_length <- nchar(seq)

  primary_window_used <- min(
    as.integer(primary_window),
    seq_length
  )

  rescue_window_used <- min(
    as.integer(rescue_window),
    seq_length
  )


  # ============================================================
  # 4. Buscar runs válidos dentro de UNA región
  #
  # Una señal pasa si cumple:
  #
  # A) exact:
  #      target_count >= exact_run
  #      mismatch_count == 0
  #
  # OR
  #
  # B) tolerant:
  #      target_count >= tolerant_run
  #      mismatch_count <= mismatch_tolerance
  #
  # ============================================================

  find_valid_runs <- function(
      region_seq,
      target_base,
      region_start,
      region_name,
      exact_run,
      tolerant_run,
      mismatch_tolerance
  ) {

    chars <- strsplit(
      region_seq,
      "",
      fixed = TRUE
    )[[1]]

    target_positions <- which(
      chars == target_base
    )

    minimum_required <- min(
      exact_run,
      tolerant_run
    )

    if (length(target_positions) < minimum_required) {
      return(NULL)
    }


    candidates <- list()
    candidate_id <- 0L


    # ----------------------------------------------------------
    # Cada candidato comienza y termina en la base objetivo.
    #
    # Como mismatch_count solamente puede aumentar al ampliar
    # un intervalo, podemos cortar cuando deja de ser posible
    # tanto el criterio exacto como el tolerante.
    # ----------------------------------------------------------

    for (i in seq_along(target_positions)) {

      start_pos <- target_positions[i]

      for (j in i:length(target_positions)) {

        end_pos <- target_positions[j]

        target_count <- j - i + 1L

        segment_length <-
          end_pos - start_pos + 1L

        mismatch_count <-
          segment_length - target_count


        qualifies_exact <-
          target_count >= exact_run &&
          mismatch_count == 0L


        qualifies_tolerant <-
          target_count >= tolerant_run &&
          mismatch_count <= mismatch_tolerance


        if (
          qualifies_exact ||
          qualifies_tolerant
        ) {

          candidate_id <- candidate_id + 1L

          candidates[[candidate_id]] <- data.frame(

            region = region_name,

            base = target_base,

            run = paste0(
              chars[start_pos:end_pos],
              collapse = ""
            ),

            run_length = segment_length,

            target_count = target_count,

            mismatch_count = mismatch_count,

            purity =
              target_count /
              segment_length,

            start =
              region_start +
              start_pos - 1L,

            end =
              region_start +
              end_pos - 1L,

            stringsAsFactors = FALSE
          )
        }


        # Si ya hay demasiados mismatches para el modo tolerante,
        # ampliar este mismo intervalo no puede volver a hacerlo válido.
        #
        # Pero solo podemos cortar después de alcanzar al menos
        # exact_run bases objetivo; los siguientes inicios se probarán
        # normalmente en la siguiente iteración de i.

        if (
          mismatch_count >
            mismatch_tolerance &&
          target_count >= exact_run
        ) {
          break
        }
      }
    }


    if (length(candidates) == 0) {
      return(NULL)
    }


    hits <- do.call(
      rbind,
      candidates
    )


    # Eliminar duplicados exactos
    hits <- unique(hits)


    # ----------------------------------------------------------
    # Conservar intervalos máximos:
    # eliminar candidatos completamente contenidos en otro.
    # ----------------------------------------------------------

    if (nrow(hits) > 1) {

      contained <- vapply(

        seq_len(nrow(hits)),

        function(i) {

          any(

            seq_len(nrow(hits)) != i &

            hits$start <= hits$start[i] &

            hits$end >= hits$end[i] &

            (
              hits$start < hits$start[i] |
              hits$end > hits$end[i]
            )
          )
        },

        logical(1)
      )

      hits <- hits[
        !contained,
        ,
        drop = FALSE
      ]
    }


    hits <- hits[
      order(
        hits$start,
        -hits$target_count,
        hits$mismatch_count,
        -hits$purity,
        -hits$run_length
      ),
      ,
      drop = FALSE
    ]

    rownames(hits) <- NULL

    hits
  }


  # ============================================================
  # 5. Buscar una combinación BASE + EXTREMO
  #
  # Primero:
  #   ventana primaria 50 bp
  #
  # Si NO encuentra nada:
  #   rescue hasta 100 bp
  #
  # IMPORTANTE:
  # el rescue examina hasta 100 bp completos.
  # Por tanto también puede rescatar dentro de los primeros
  # 50 bp un tracto que fallase primaria por tener 3 mismatches.
  # Esto reproduce la regla validada en HG002.
  # ============================================================

  search_one_end <- function(
      target_base,
      region_name
  ) {

    # ----------------------------------------------------------
    # PRIMARY
    # ----------------------------------------------------------

    if (region_name == "initial") {

      primary_seq <-
        substr(
          seq,
          1L,
          primary_window_used
        )

      primary_start <- 1L

    } else {

      primary_start <-
        seq_length -
        primary_window_used +
        1L

      primary_seq <-
        substr(
          seq,
          primary_start,
          seq_length
        )
    }


    primary_hits <- find_valid_runs(

      region_seq = primary_seq,

      target_base = target_base,

      region_start = primary_start,

      region_name = region_name,

      exact_run =
        primary_exact_run,

      tolerant_run =
        primary_tolerant_run,

      mismatch_tolerance =
        primary_mismatch_tolerance
    )


    if (!is.null(primary_hits)) {
      return(primary_hits)
    }


    # ----------------------------------------------------------
    # RESCUE
    # ----------------------------------------------------------

    if (region_name == "initial") {

      rescue_seq <-
        substr(
          seq,
          1L,
          rescue_window_used
        )

      rescue_start <- 1L

    } else {

      rescue_start <-
        seq_length -
        rescue_window_used +
        1L

      rescue_seq <-
        substr(
          seq,
          rescue_start,
          seq_length
        )
    }


    rescue_hits <- find_valid_runs(

      region_seq = rescue_seq,

      target_base = target_base,

      region_start = rescue_start,

      region_name = region_name,

      exact_run =
        rescue_exact_run,

      tolerant_run =
        rescue_tolerant_run,

      mismatch_tolerance =
        rescue_mismatch_tolerance
    )


    rescue_hits
  }


  # ============================================================
  # 6. Ejecutar para A/T × initial/final
  # ============================================================

  all_hits <- list()
  hit_id <- 0L


  for (base in bases) {

    initial_hits <- search_one_end(
      target_base = base,
      region_name = "initial"
    )

    if (!is.null(initial_hits)) {

      hit_id <- hit_id + 1L

      all_hits[[hit_id]] <-
        initial_hits
    }


    final_hits <- search_one_end(
      target_base = base,
      region_name = "final"
    )

    if (!is.null(final_hits)) {

      hit_id <- hit_id + 1L

      all_hits[[hit_id]] <-
        final_hits
    }
  }


  # ============================================================
  # 7. Construir tabla de hits
  #
  # EXACTAMENTE las mismas columnas que antes.
  # ============================================================

  if (length(all_hits) == 0) {

    hits <- data.frame(

      region = character(),

      base = character(),

      run = character(),

      run_length = integer(),

      target_count = integer(),

      mismatch_count = integer(),

      purity = numeric(),

      start = integer(),

      end = integer(),

      stringsAsFactors = FALSE
    )

  } else {

    hits <- do.call(
      rbind,
      all_hits
    )

    hits <- hits[
      order(
        hits$start,
        hits$base
      ),
      ,
      drop = FALSE
    ]

    rownames(hits) <- NULL
  }


  # ============================================================
  # 8. Salida
  #
  # Se mantienen los mismos nombres para no romper
  # el resto del pipeline.
  # ============================================================

  list(

    sequence_length =
      seq_length,

    # Compatibilidad con salida antigua.
    # Ya no se utiliza una fracción dinámica.
    fraction_examined =
      if (is.null(fraction)) {
        NA_real_
      } else {
        fraction
      },

    # Ahora representa la ventana máxima posible de búsqueda.
    window_size =
      rescue_window_used,

    # Compatibilidad con nombres antiguos:
    # reflejan el modo tolerante PRIMARY.
    min_run =
      primary_tolerant_run,

    mismatch_tolerance =
      primary_mismatch_tolerance,

    # Para mantener exactamente el tipo de información anterior,
    # exponemos las secuencias de la ventana máxima.
    initial_sequence =
      substr(
        seq,
        1L,
        rescue_window_used
      ),

    final_sequence =
      substr(
        seq,
        seq_length -
          rescue_window_used +
          1L,
        seq_length
      ),

    hits =
      hits,

    has_polyA =
      any(
        hits$base == "A"
      ),

    has_polyT =
      any(
        hits$base == "T"
      ),

    possible_direct_insertion =
      any(
        hits$region == "final" &
        hits$base == "A"
      ),

    possible_inverted_insertion =
      any(
        hits$region == "initial" &
        hits$base == "T"
      )
  )
}

extraer_campo <- function(columna_info, etiqueta) {
  # Buscamos la etiqueta, capturamos lo que sigue al '=' hasta el ';' o final
  pattern <- paste0(".*", etiqueta, "=([^;]+).*")
  # Si no encuentra la etiqueta, devuelve NA
  resultado <- sub(pattern, "\\1", columna_info)
  
  
  # Limpieza: si el campo no existía, sub devuelve la cadena original, así que lo corregimos
  resultado[resultado == columna_info] <- NA 
  return(resultado)
}

source(arg_value("--reference-helpers"))
mapping <- read_chromosome_map(arg_value("--chromosome-map"))
sex <- arg_value("--sex", "U")
procedence <- arg_value("--procedence", "UNKNOW")
assembly <- arg_value("--assembly")
panel <- read_panel(arg_value("--panel-bed"))
empty_result <- function(reason) {
  message("No reconstructed loci: ", reason)
  loci <- data.frame(IDs=character(), CHROM=character(), POS=numeric(), svlen=numeric(),
                     GT=character(), contig=character(), cover_status=character(), region_status=character())
  saveRDS(list(format_version="RetroAtlas_single_v1", sample_id=sample, loci=loci,
               fragment_evidence=data.frame(), rescue_candidates=loci, nocallablecandidates=loci, parameters=list(CFG=CFG, FRAG=FRAG),
               counts=list(n_loci=0L), empty_reason=reason),
          file.path(out, paste0(sample, "_retroatlas_intra_final.rds")))
  writeLines(character(), file.path(out, paste0(sample, "_samtools.bed")))
  writeLines(character(), file.path(out, paste0(sample, "_contig.fasta")))
  quit(save="no", status=0)
}
predata=readRDS(predata)
normalized <- standardize_loci(predata, mapping, sample, sex, procedence, assembly, panel)
coverage <- read_coverage(Regions_data, mapping, CFG$min_depth, CFG$max_depth_ratio, sex, panel)
predata <- annotate_coverage(normalized$kept, coverage)
excluded <- dplyr::bind_rows(normalized$excluded, predata[!is.na(predata$EXCLUSION_REASON),,drop=FALSE])
write_flat_tsv(excluded, file.path(out, paste0(sample, "_excluded_predata.tsv")))
# Reconstruct validated positions even when non-callable; partition final loci below.
if (!nrow(predata)) empty_result("No validated primary-chromosome insertions")
# Keep the classifier working only on eligible query IDs, including zero-hit cases.
rm_lines <- readLines(RM_data, warn=FALSE)
rm_rows <- rm_lines[grepl("^\\s*[0-9]+\\s+", rm_lines)]
if (!length(rm_rows)) {
  if (length(rm_lines) && !any(grepl("no repetitive|SW|score", rm_lines, ignore.case=TRUE))) stop("Malformed RepeatMasker output")
  empty_result("No RepeatMasker hits")
}

predata$contig=extraer_campo(predata$INFO,"CONTIG")



flank <- 0.10
umbral <- 100

# Inicializar columnas
predata$haspoly_T <- FALSE
predata$haspoly_A <- FALSE

predata$length_poly_T_5 <- NA_integer_
predata$length_poly_T_3 <- NA_integer_
predata$mismatches_poly_T_5 <- NA_integer_
predata$mismatches_poly_T_3 <- NA_integer_
predata$purity_poly_T_5 <- NA_real_
predata$purity_poly_T_3 <- NA_real_

predata$length_poly_A_5 <- NA_integer_
predata$length_poly_A_3 <- NA_integer_
predata$mismatches_poly_A_5 <- NA_integer_
predata$mismatches_poly_A_3 <- NA_integer_
predata$purity_poly_A_5 <- NA_real_
predata$purity_poly_A_3 <- NA_real_


for (i in seq_len(nrow(predata))) {
  
  sequence <- as.character(predata$ALT[i])
  
  # Saltar secuencias ausentes o vacías
  if (
    is.na(sequence) ||
    !nzchar(sequence)
  ) {
    next
  }
  
  sequence_length <- nchar(sequence)
  
  # Fracción dinámica, con máximo de 'umbral' nucleótidos
  current_fraction <- min(
    flank,
    umbral / sequence_length
  )
  
result <- detect_terminal_polyAT(
  seq = sequence,
  primary_window = CFG$polyat_primary_window,
  primary_exact_run = CFG$polyat_primary_exact_run,
  primary_tolerant_run = CFG$polyat_primary_tolerant_run,
  primary_mismatch_tolerance = CFG$polyat_primary_mismatch_tolerance,
  rescue_window = CFG$polyat_rescue_window,
  rescue_exact_run = CFG$polyat_rescue_exact_run,
  rescue_tolerant_run = CFG$polyat_rescue_tolerant_run,
  rescue_mismatch_tolerance = CFG$polyat_rescue_mismatch_tolerance
)
  
  predata$haspoly_T[i] <- result$has_polyT
  predata$haspoly_A[i] <- result$has_polyA
  
  df <- as.data.frame(result$hits)
  
  # No hay hits
  if (nrow(df) == 0) {
    next
  }
  
  # Elegir el mejor hit de cada combinación base-región
  df_best <- df %>%
    dplyr::group_by(base, region) %>%
    dplyr::arrange(
      dplyr::desc(target_count),
      mismatch_count,
      dplyr::desc(purity),
      dplyr::desc(run_length),
      .by_group = TRUE
    ) %>%
    dplyr::slice_head(n = 1) %>%
    dplyr::ungroup()
  
  # T en el extremo inicial, equivalente al extremo 5'
  hit <- df_best %>%
    dplyr::filter(
      base == "T",
      region == "initial"
    )
  
  if (nrow(hit) == 1) {
    predata$length_poly_T_5[i] <- hit$target_count[1]
    predata$mismatches_poly_T_5[i] <- hit$mismatch_count[1]
    predata$purity_poly_T_5[i] <- hit$purity[1]
  }
  
  # T en el extremo final, equivalente al extremo 3'
  hit <- df_best %>%
    dplyr::filter(
      base == "T",
      region == "final"
    )
  
  if (nrow(hit) == 1) {
    predata$length_poly_T_3[i] <- hit$target_count[1]
    predata$mismatches_poly_T_3[i] <- hit$mismatch_count[1]
    predata$purity_poly_T_3[i] <- hit$purity[1]
  }
  
  # A en el extremo inicial, equivalente al extremo 5'
  hit <- df_best %>%
    dplyr::filter(
      base == "A",
      region == "initial"
    )
  
  if (nrow(hit) == 1) {
    predata$length_poly_A_5[i] <- hit$target_count[1]
    predata$mismatches_poly_A_5[i] <- hit$mismatch_count[1]
    predata$purity_poly_A_5[i] <- hit$purity[1]
  }
  
  # A en el extremo final, equivalente al extremo 3'
  hit <- df_best %>%
    dplyr::filter(
      base == "A",
      region == "final"
    )
  
  if (nrow(hit) == 1) {
    predata$length_poly_A_3[i] <- hit$target_count[1]
    predata$mismatches_poly_A_3[i] <- hit$mismatch_count[1]
    predata$purity_poly_A_3[i] <- hit$purity[1]
  }
}



 
# 3. Enriquecimiento RepeatMasker--------------------------------------

source(classifier_script)

result <- run_te_classifier(
  repeatmasker_out = RM_data,
  output_prefix = sample,
  output_dir = out,
  query_threshold = CFG$query_threshold,
  model_threshold = CFG$model_threshold,
  group_dominance_threshold = CFG$group_dominance_threshold,
  write_segments = write_audits,
  config_overrides = yaml::read_yaml(arg_value("--classifier-config")),
  query_ids = predata$header
)

final_RM=result$calls
if (!nrow(final_RM)) empty_result("No eligible RepeatMasker hits")


# 4. Fusion----
data=predata%>%rename(query=header)
FinalTable=left_join(final_RM,data,by="query")
FinalTable <- FinalTable %>%
  separate_wider_delim(
    cols = GENOTYPE,
    delim = ":",
    names = c("GT", "GQ", "LR", "SR", "SL", "IF")
  )
FinalTable$FORMAT=NULL
FinalTable$length_poly_T_3=NULL
FinalTable$length_poly_A_5=NULL
FinalTable$mismatches_poly_A_5=NULL
FinalTable$mismatches_poly_T_3=NULL
FinalTable$purity_poly_T_3=NULL
FinalTable$purity_poly_A_5=NULL
FinalTable$mismatches_poly_T_3=NULL
FinalTable$haspoly_A=NULL
FinalTable$haspoly_T=NULL
FinalTable$ID2=NULL
FinalTable$anchor_id=NULL
FinalTable$terminal_polyA=NULL
FinalTable$terminal_polyA=NULL
FinalTable$terminal_polyA_length=NULL
FinalTable$terminal_polyA_purity=NULL
FinalTable=FinalTable%>%mutate(polyA=case_when(
  !is.na(length_poly_T_5) ~ "TRUE",
  !is.na(length_poly_A_3) ~ "TRUE",
  # Reemplaza "No es NA" por el valor o variable que quieras asignar
  TRUE                     ~ "FALSE"     # (Opcional) Equivalente al 'else' para capturar el resto de casos
),polyA_length=case_when(
  !is.na(length_poly_T_5) ~ length_poly_T_5,
  !is.na(length_poly_A_3) ~ length_poly_A_3,
  # Reemplaza "No es NA" por el valor o variable que quieras asignar
  TRUE                     ~ 0     # (Opcional) Equivalente al 'else' para capturar el resto de casos
))





checksize <- function(a, b) {

  a <- abs(as.numeric(a))
  b <- abs(as.numeric(b))

  if (
    is.na(a) ||
    is.na(b) ||
    a == 0 ||
    b == 0
  ) {
    return(Inf)
  }

  ((max(a, b) / min(a, b)) - 1) * 100
}


# ==============================================================================
# FRAGMENTATION-AWARE AUDIT HELPERS
# ==============================================================================

frag_resolved <- function(x) {
  z <- tolower(trimws(as.character(x)))
  !is.na(z) & !(z %in% c(
    "", ".", "na", "nan", "none", "null", "unknown",
    "unresolved", "not_resolved", "not resolved"
  ))
}

frag_safe_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

frag_size_similarity <- function(a, b) {
  a <- abs(frag_safe_num(a))
  b <- abs(frag_safe_num(b))
  n <- max(length(a), length(b))
  if (!length(a)) a <- rep(NA_real_, n)
  if (!length(b)) b <- rep(NA_real_, n)
  if (length(a) == 1L && n > 1L) a <- rep(a, n)
  if (length(b) == 1L && n > 1L) b <- rep(b, n)
  out <- rep(NA_real_, n)
  ok <- is.finite(a) & is.finite(b) & a > 0 & b > 0
  out[ok] <- pmin(a[ok], b[ok]) / pmax(a[ok], b[ok])
  out
}

frag_mode_chr <- function(x) {
  z <- as.character(x)
  z <- z[!is.na(z) & nzchar(z)]
  if (!length(z)) return(NA_character_)
  tt <- table(z)
  sort(names(tt)[tt == max(tt)])[1L]
}

frag_collapse_set <- function(x) {
  z <- sort(unique(as.character(x)))
  z <- z[!is.na(z) & nzchar(z)]
  if (!length(z)) return(NA_character_)
  paste(z, collapse = ";")
}

frag_normalize_gt <- function(x) {
  x <- as.character(x)
  x[x %in% c("1|0", "0|1", "1/0")] <- "0/1"
  x[x %in% c("1|1")] <- "1/1"
  x[x %in% c("0|0")] <- "0/0"
  x
}

frag_gt_alt <- function(x) {
  z <- frag_normalize_gt(x)
  z %in% c("0/1", "1/1")
}

fragment_source_id <- function(query) {
  q <- as.character(query)
  out <- sub(";.*$", "", q)
  out[is.na(q) | !nzchar(q)] <- NA_character_
  out
}

fragment_haplotype <- function(IF) {
  x <- as.character(IF)
  out <- stringr::str_extract(x, "^[^,]+,[^,]+")
  out[is.na(x) | !nzchar(x)] <- NA_character_
  out
}

load_simple_repeats <- function(path) {
  if (is.na(path) || !nzchar(path)) {
    message("SimpleRepeats BED: not supplied")
    return(data.table())
  }
  if (!file.exists(path)) stop("SimpleRepeats BED not found: ", path)

  message("Reading SimpleRepeats BED: ", path)
  x <- data.table::fread(
    path,
    header = FALSE,
    showProgress = FALSE
  )
  if (ncol(x) < 3L) stop("SimpleRepeats BED must contain at least 3 columns")

  while (ncol(x) < 5L) x[, (paste0("V", ncol(x) + 1L)) := NA_character_]
  data.table::setnames(x, names(x)[1:5], c("CHROM", "START0", "END0", "REPEAT_TYPE", "REPEAT_SCORE"))
  x[, CHROM := as.character(CHROM)]
  start0 <- suppressWarnings(as.numeric(x$START0))
  end0 <- suppressWarnings(as.numeric(x$END0))
  if (any(!is.finite(start0) | !is.finite(end0) | start0 < 0 | end0 <= start0 |
          start0 != floor(start0) | end0 != floor(end0) | end0 > .Machine$integer.max))
    stop("SimpleRepeats BED contains invalid coordinates: ", path)
  x[, START := as.integer(start0) + 1L]
  x[, END := as.integer(end0)]
  x <- x[!is.na(CHROM) & !is.na(START) & !is.na(END) & END >= START]
  x[, repeat_id := .I]
  data.table::setkey(x, CHROM, START, END)
  x[]
}

annotate_simple_repeats <- function(df, repeats) {
  d <- data.table::as.data.table(data.table::copy(df))
  d[, .fragment_row_id := .I]
  d[, SIMPLE_REPEAT_OVERLAP := FALSE]
  d[, SIMPLE_REPEAT_N := 0L]
  d[, SIMPLE_REPEAT_IDS := NA_character_]
  d[, DIST_TO_SIMPLE_REPEAT_BP := NA_real_]

  if (!nrow(d) || !nrow(repeats)) return(as.data.frame(d))

  pts <- d[
    !is.na(CHROM) & !is.na(POS),
    .(
      .fragment_row_id,
      CHROM = as.character(CHROM),
      START = as.integer(round(frag_safe_num(POS))),
      END = as.integer(round(frag_safe_num(POS)))
    )
  ]
  data.table::setkey(pts, CHROM, START, END)

  ov <- data.table::foverlaps(
    pts,
    repeats,
    by.x = c("CHROM", "START", "END"),
    by.y = c("CHROM", "START", "END"),
    type = "within",
    nomatch = 0L
  )

  if (nrow(ov)) {
    ovs <- ov[
      ,
      .(
        SIMPLE_REPEAT_N = .N,
        SIMPLE_REPEAT_IDS = paste(sort(unique(repeat_id)), collapse = ";")
      ),
      by = .fragment_row_id
    ]
    d[ovs, on = ".fragment_row_id", `:=`(
      SIMPLE_REPEAT_OVERLAP = TRUE,
      SIMPLE_REPEAT_N = i.SIMPLE_REPEAT_N,
      SIMPLE_REPEAT_IDS = i.SIMPLE_REPEAT_IDS,
      DIST_TO_SIMPLE_REPEAT_BP = 0
    )]
  }

  # Nearest-repeat distance for non-overlapping points.  prefix max END handles
  # nested/overlapping BED intervals without requiring interval reduction.
  for (chr_now in unique(as.character(d$CHROM))) {
    if (is.na(chr_now)) next
    rr <- repeats[CHROM == chr_now]
    if (!nrow(rr)) next
    data.table::setorder(rr, START, END)
    starts <- as.numeric(rr$START)
    prefix_end <- cummax(as.numeric(rr$END))

    idx_rows <- which(as.character(d$CHROM) == chr_now & !is.na(d$POS))
    if (!length(idx_rows)) next
    pos <- frag_safe_num(d$POS[idx_rows])
    ii <- findInterval(pos, starts)

    dist <- rep(Inf, length(pos))
    has_prev <- ii > 0L
    if (any(has_prev)) {
      prev_end <- prefix_end[ii[has_prev]]
      dist[has_prev] <- pmax(0, pos[has_prev] - prev_end)
    }
    has_next <- ii < length(starts)
    if (any(has_next)) {
      next_start <- starts[ii[has_next] + 1L]
      dist[has_next] <- pmin(dist[has_next], pmax(0, next_start - pos[has_next]))
    }
    d$DIST_TO_SIMPLE_REPEAT_BP[idx_rows] <- ifelse(
      d$SIMPLE_REPEAT_OVERLAP[idx_rows],
      0,
      dist
    )
  }

  as.data.frame(d)
}

frag_extract_info <- function(info, key) {
  pat <- paste0("(?:^|;)", key, "=([^;]+)")
  m <- regexec(pat, as.character(info), perl = TRUE)
  z <- regmatches(as.character(info), m)
  vapply(z, function(x) if (length(x) >= 2L) x[[2L]] else NA_character_, character(1))
}

frag_extract_format <- function(fmt, sample_field, key) {
  vapply(seq_along(fmt), function(i) {
    f <- as.character(fmt[i])
    s <- as.character(sample_field[i])
    if (is.na(f) || is.na(s)) return(NA_character_)
    fn <- strsplit(f, ":", fixed = TRUE)[[1L]]
    sv <- strsplit(s, ":", fixed = TRUE)[[1L]]
    k <- match(key, fn)
    if (is.na(k) || k > length(sv)) return(NA_character_)
    sv[[k]]
  }, character(1))
}

read_fragment_vcf <- function(path, caller) {
  if (is.na(path) || !nzchar(path) || !file.exists(path)) return(data.table())
  message("Fragmentation VCF corroboration [", caller, "]: ", path)

  cmd <- if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    sprintf("zcat %s | grep -v '^##'", shQuote(path))
  } else {
    sprintf("cat %s | grep -v '^##'", shQuote(path))
  }

  v <- tryCatch(
    data.table::fread(cmd = cmd, sep = "\t", header = TRUE, fill = TRUE, showProgress = FALSE),
    error = function(e) {
      warning("Could not read fragmentation VCF [", caller, "]: ", e$message)
      NULL
    }
  )
  if (is.null(v) || !nrow(v)) return(data.table())

  chrom_col <- if ("#CHROM" %in% names(v)) "#CHROM" else if ("CHROM" %in% names(v)) "CHROM" else NA_character_
  if (is.na(chrom_col) || !all(c("POS", "INFO") %in% names(v))) return(data.table())

  svtype <- frag_extract_info(v$INFO, "SVTYPE")
  svlen_txt <- frag_extract_info(v$INFO, "SVLEN")
  svlen <- vapply(strsplit(as.character(svlen_txt), ",", fixed = TRUE), function(z) {
    suppressWarnings(as.numeric(z[[1L]]))
  }, numeric(1))

  miss <- is.na(svlen)
  if (any(miss) && all(c("REF", "ALT") %in% names(v))) {
    nonsym <- miss & !grepl("^<", as.character(v$ALT))
    svlen[nonsym] <- nchar(as.character(v$ALT[nonsym])) - nchar(as.character(v$REF[nonsym]))
  }

  GT <- rep(NA_character_, nrow(v))
  if ("FORMAT" %in% names(v)) {
    fixed <- c(chrom_col, "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT")
    sample_cols <- setdiff(names(v), fixed)
    if (length(sample_cols)) {
      GT <- frag_extract_format(v$FORMAT, v[[sample_cols[[1L]]]], "GT")
    }
  }

  out <- data.table(
    caller = caller,
    CHROM = as.character(v[[chrom_col]]),
    POS = frag_safe_num(v$POS),
    SVTYPE = as.character(svtype),
    SVLEN = abs(frag_safe_num(svlen)),
    GT = GT,
    ALT_GT = frag_gt_alt(GT),
    FILTER = if ("FILTER" %in% names(v)) as.character(v$FILTER) else NA_character_,
    ID = if ("ID" %in% names(v)) as.character(v$ID) else NA_character_
  )

  out <- out[
    !is.na(CHROM) & !is.na(POS) & !is.na(SVLEN) &
      SVTYPE == "INS" & SVLEN >= FRAG$vcf_min_svlen
  ]
  if (FRAG$vcf_require_pass) {
    out <- out[is.na(FILTER) | FILTER %in% c("PASS", ".")]
  }
  out[]
}

load_fragment_vcfs <- function() {
  specs <- list(cutesv = vcf_cutesv, severus = vcf_severus, sniffles2 = vcf_sniffles2)
  zz <- lapply(names(specs), function(nm) read_fragment_vcf(specs[[nm]], nm))
  zz <- zz[vapply(zz, nrow, integer(1)) > 0L]
  if (!length(zz)) return(data.table())
  data.table::rbindlist(zz, fill = TRUE)
}

build_fragment_group_audit <- function(prefusion_all, fragment_vcf) {
  d <- data.table::as.data.table(data.table::copy(prefusion_all))
  if (!nrow(d)) return(data.table())

  if (!"structural_pass" %in% names(d)) d[, structural_pass := NA]
  d[, .source_event_id := fragment_source_id(query)]
  d[, .fragment_haplotype := fragment_haplotype(IF)]
  d[, .fragment_gt := frag_normalize_gt(GT)]
  d[, .fragment_svlen := abs(frag_safe_num(svlen))]
  d[, .fragment_resolved := frag_resolved(event_class) & frag_resolved(event)]

  z <- d[
    !is.na(CHROM) & !is.na(POS) & !is.na(.fragment_svlen) & .fragment_svlen > 0 &
      !is.na(.source_event_id) & nzchar(.source_event_id) &
      !is.na(.fragment_haplotype) & nzchar(.fragment_haplotype)
  ]
  if (!nrow(z)) return(data.table())

  groups <- z[
    ,
    {
      pp <- sort(unique(frag_safe_num(POS)))
      gaps <- if (length(pp) >= 2L) diff(pp) else numeric()
      cls <- unique(as.character(event_class[.fragment_resolved]))
      cls <- cls[frag_resolved(cls)]
      evt <- unique(as.character(event[.fragment_resolved]))
      evt <- evt[frag_resolved(evt)]
      gtset <- unique(as.character(.fragment_gt))
      gtset <- gtset[!is.na(gtset) & nzchar(gtset)]
      rep_dist <- frag_safe_num(DIST_TO_SIMPLE_REPEAT_BP)
      rep_dist <- rep_dist[is.finite(rep_dist)]
      lr_vals <- if ("LR" %in% names(.SD)) frag_safe_num(LR) else rep(NA_real_, .N)
      lr_ok <- lr_vals[is.finite(lr_vals)]
      resolved_sv <- .fragment_svlen[.fragment_resolved & is.finite(.fragment_svlen)]
      resolved_sum <- if (length(resolved_sv)) sum(resolved_sv, na.rm = TRUE) else 0
      anchor_max <- if (length(resolved_sv)) max(resolved_sv, na.rm = TRUE) else NA_real_
      compound_sum <- sum(.fragment_svlen, na.rm = TRUE)

      list(
        N_COMPONENTS = .N,
        N_RESOLVED_COMPONENTS = sum(.fragment_resolved, na.rm = TRUE),
        N_UNRESOLVED_COMPONENTS = sum(!.fragment_resolved, na.rm = TRUE),
        N_STRUCTURAL_PASS_FALSE = if (all(is.na(structural_pass))) {
          0L
        } else {
          sum(!is.na(structural_pass) & !as.logical(structural_pass), na.rm = TRUE)
        },
        POS_MIN_COMPONENT = min(frag_safe_num(POS), na.rm = TRUE),
        POS_MAX_COMPONENT = max(frag_safe_num(POS), na.rm = TRUE),
        COMPONENT_SPAN_BP = max(frag_safe_num(POS), na.rm = TRUE) - min(frag_safe_num(POS), na.rm = TRUE),
        MIN_COMPONENT_GAP_BP = if (length(gaps)) min(gaps) else NA_real_,
        MAX_COMPONENT_GAP_BP = if (length(gaps)) max(gaps) else NA_real_,
        COMPOUND_SVLEN = compound_sum,
        RESOLVED_SVLEN_SUM = resolved_sum,
        ANCHOR_MAX_SVLEN = anchor_max,
        ANCHOR_FRACTION = if (is.finite(anchor_max) && compound_sum > 0) anchor_max / compound_sum else NA_real_,
        MIN_COMPONENT_LR = if (length(lr_ok) == .N) min(lr_ok) else NA_real_,
        COMPONENT_POS_VECTOR = paste(frag_safe_num(POS)[order(frag_safe_num(POS))], collapse = ";"),
        COMPONENT_SVLEN_VECTOR = paste(.fragment_svlen[order(frag_safe_num(POS))], collapse = ";"),
        COMPONENT_QUERY_VECTOR = paste(as.character(query)[order(frag_safe_num(POS))], collapse = ";"),
        COMPONENT_ID_VECTOR = paste(as.character(ID)[order(frag_safe_num(POS))], collapse = ";"),
        COMPONENT_GT_SET = frag_collapse_set(gtset),
        RESOLVED_CLASS_SET = frag_collapse_set(cls),
        RESOLVED_EVENT_SET = frag_collapse_set(evt),
        SAME_RESOLVED_CLASS = length(cls) <= 1L && length(cls) >= 1L,
        SAME_RESOLVED_EVENT = length(evt) <= 1L && length(evt) >= 1L,
        ANY_SIMPLE_REPEAT = any(SIMPLE_REPEAT_OVERLAP %in% TRUE, na.rm = TRUE),
        MIN_DIST_SIMPLE_REPEAT_BP = if (length(rep_dist)) min(rep_dist) else NA_real_
      )
    },
    by = .(CHROM, .source_event_id, .fragment_haplotype)
  ]

  groups[, FRAGMENT_GROUP_ID := paste(CHROM, .source_event_id, .fragment_haplotype, sep = "||")]
  groups[, GROUP_POS := (POS_MIN_COMPONENT + POS_MAX_COMPONENT) / 2]

  # VCF corroboration is evaluated against the SUM of components, not each
  # individual fragment.  A caller can contribute at most once per group.
  groups[, `:=`(
    BEST_VCF_ALT_SIMILARITY = NA_real_,
    BEST_VCF_ALT_DPOS = NA_real_,
    BEST_VCF_ALT_SVLEN = NA_real_,
    BEST_VCF_ALT_CALLER = NA_character_,
    N_VCF_ALT_CALLERS_GE_MIN = 0L,
    N_VCF_ALT_CALLERS_GE_STRONG = 0L
  )]

  if (nrow(fragment_vcf)) {
    gi <- groups[
      N_COMPONENTS >= 2L,
      .(
        FRAGMENT_GROUP_ID,
        CHROM,
        START = as.integer(floor(GROUP_POS - FRAG$vcf_window_bp)),
        END = as.integer(ceiling(GROUP_POS + FRAG$vcf_window_bp)),
        GROUP_POS,
        COMPOUND_SVLEN
      )
    ]
    vv <- data.table::copy(fragment_vcf)
    vv[, START := as.integer(round(POS))]
    vv[, END := START]

    mt <- data.table()
    if (nrow(gi)) {
      data.table::setkey(gi, CHROM, START, END)
      data.table::setkey(vv, CHROM, START, END)

      mt <- data.table::foverlaps(
        vv,
        gi,
        by.x = c("CHROM", "START", "END"),
        by.y = c("CHROM", "START", "END"),
        type = "within",
        nomatch = 0L
      )
    }

    if (nrow(mt)) {
      mt[, dPOS := abs(POS - GROUP_POS)]
      mt[, compound_similarity := frag_size_similarity(SVLEN, COMPOUND_SVLEN)]
      mt <- mt[ALT_GT %in% TRUE & !is.na(compound_similarity)]

      if (nrow(mt)) {
        data.table::setorder(mt, FRAGMENT_GROUP_ID, caller, -compound_similarity, dPOS)
        best_caller <- mt[, .SD[1L], by = .(FRAGMENT_GROUP_ID, caller)]
        sums <- best_caller[
          ,
          {
            ii <- order(-compound_similarity, dPOS)[1L]
            list(
              BEST_VCF_ALT_SIMILARITY = compound_similarity[ii],
              BEST_VCF_ALT_DPOS = dPOS[ii],
              BEST_VCF_ALT_SVLEN = SVLEN[ii],
              BEST_VCF_ALT_CALLER = caller[ii],
              N_VCF_ALT_CALLERS_GE_MIN = uniqueN(caller[compound_similarity >= FRAG$vcf_min_similarity]),
              N_VCF_ALT_CALLERS_GE_STRONG = uniqueN(caller[compound_similarity >= FRAG$vcf_strong_similarity])
            )
          },
          by = FRAGMENT_GROUP_ID
        ]
        groups[sums, on = "FRAGMENT_GROUP_ID", `:=`(
          BEST_VCF_ALT_SIMILARITY = i.BEST_VCF_ALT_SIMILARITY,
          BEST_VCF_ALT_DPOS = i.BEST_VCF_ALT_DPOS,
          BEST_VCF_ALT_SVLEN = i.BEST_VCF_ALT_SVLEN,
          BEST_VCF_ALT_CALLER = i.BEST_VCF_ALT_CALLER,
          N_VCF_ALT_CALLERS_GE_MIN = i.N_VCF_ALT_CALLERS_GE_MIN,
          N_VCF_ALT_CALLERS_GE_STRONG = i.N_VCF_ALT_CALLERS_GE_STRONG
        )]
      }
    }
  }

  groups[, REPEAT_CONTEXT := (
    ANY_SIMPLE_REPEAT |
      (!is.na(MIN_DIST_SIMPLE_REPEAT_BP) & MIN_DIST_SIMPLE_REPEAT_BP <= FRAG$repeat_near_bp)
  )]

  groups[, FRAGMENTATION_RISK := data.table::fcase(
    N_COMPONENTS <= 1L, "LOW",
    REPEAT_CONTEXT | (!is.na(MIN_COMPONENT_GAP_BP) & MIN_COMPONENT_GAP_BP <= FRAG$local_gap_bp), "HIGH",
    N_COMPONENTS >= 2L, "MEDIUM",
    default = "LOW"
  )]

  groups[, AUTO_LOCAL_GEOMETRY_OK := (
    N_COMPONENTS >= 2L &
      N_COMPONENTS <= FRAG$auto_max_components &
      COMPONENT_SPAN_BP <= FRAG$auto_max_span_bp &
      !is.na(MAX_COMPONENT_GAP_BP) & MAX_COMPONENT_GAP_BP <= FRAG$local_gap_bp
  )]

  # A wider fragmented representation can still be high-confidence when >=2
  # independent callers observe an ALT insertion whose size matches the SUM.
  groups[, AUTO_VCF_GEOMETRY_OK := (
    N_COMPONENTS >= 2L &
      N_COMPONENTS <= FRAG$auto_max_components &
      COMPONENT_SPAN_BP <= FRAG$auto_max_span_bp &
      N_VCF_ALT_CALLERS_GE_STRONG >= FRAG$auto_vcf_min_callers
  )]

  groups[, AUTO_GEOMETRY_OK := AUTO_LOCAL_GEOMETRY_OK | AUTO_VCF_GEOMETRY_OK]

  groups[, AUTO_ANNOTATION_OK := (
    N_UNRESOLVED_COMPONENTS == 0L &
      N_STRUCTURAL_PASS_FALSE == 0L &
      SAME_RESOLVED_CLASS
  )]

  groups[, ANCHORED_DROPOUT_CANDIDATE := (
    FRAG$anchored_dropout_enabled &
      N_COMPONENTS >= 2L &
      N_COMPONENTS <= FRAG$auto_max_components &
      N_RESOLVED_COMPONENTS >= 1L &
      N_UNRESOLVED_COMPONENTS >= 1L &
      N_UNRESOLVED_COMPONENTS <= FRAG$anchored_max_unresolved &
      N_STRUCTURAL_PASS_FALSE == 0L &
      SAME_RESOLVED_CLASS &
      (!FRAG$anchored_require_same_event | SAME_RESOLVED_EVENT) &
      COMPONENT_SPAN_BP <= FRAG$auto_max_span_bp &
      !is.na(ANCHOR_FRACTION) & ANCHOR_FRACTION >= FRAG$anchored_min_anchor_fraction &
      !is.na(MIN_COMPONENT_LR) & MIN_COMPONENT_LR >= FRAG$anchored_min_component_lr &
      N_VCF_ALT_CALLERS_GE_STRONG >= FRAG$anchored_vcf_min_callers &
      !is.na(BEST_VCF_ALT_SIMILARITY) &
      BEST_VCF_ALT_SIMILARITY >= FRAG$anchored_vcf_similarity
  )]

  groups[, AUTO_VCF_OK := (
    !FRAG$auto_require_vcf |
      N_VCF_ALT_CALLERS_GE_STRONG >= 1L
  )]

  groups[, FRAGMENT_AUTO_CANDIDATE := AUTO_GEOMETRY_OK & AUTO_ANNOTATION_OK & AUTO_VCF_OK]
  groups[, FRAGMENT_RECONSTRUCT_CANDIDATE := FRAGMENT_AUTO_CANDIDATE | ANCHORED_DROPOUT_CANDIDATE]

  groups[, FRAGMENTATION_CLASS := data.table::fcase(
    N_COMPONENTS <= 1L, "SINGLE_COMPONENT",
    FRAGMENT_AUTO_CANDIDATE, "COMPOUND_HIGHCONF",
    ANCHORED_DROPOUT_CANDIDATE, "COMPOUND_ANCHORED_DROPOUT",
    N_COMPONENTS >= FRAG$complex_min_components |
      COMPONENT_SPAN_BP >= FRAG$complex_span_bp |
      (REPEAT_CONTEXT & N_COMPONENTS >= 3L), "SOURCE_COMPLEX",
    N_UNRESOLVED_COMPONENTS > 0L, "COMPOUND_WITH_CLASSIFIER_DROPOUT",
    default = "COMPOUND_REVIEW"
  )]

  groups[, FRAGMENT_Q := data.table::fcase(
    FRAGMENTATION_CLASS == "COMPOUND_HIGHCONF" & N_VCF_ALT_CALLERS_GE_STRONG >= 2L, 5L,
    FRAGMENTATION_CLASS == "COMPOUND_HIGHCONF" & (N_VCF_ALT_CALLERS_GE_STRONG >= 1L | ANY_SIMPLE_REPEAT), 4L,
    FRAGMENTATION_CLASS == "COMPOUND_HIGHCONF", 3L,
    FRAGMENTATION_CLASS == "COMPOUND_ANCHORED_DROPOUT", 4L,
    FRAGMENTATION_CLASS %in% c("COMPOUND_REVIEW", "COMPOUND_WITH_CLASSIFIER_DROPOUT"), 2L,
    FRAGMENTATION_CLASS == "SOURCE_COMPLEX", 1L,
    default = 5L
  )]

  groups[]
}


build_fragment_source_audit <- function(prefusion_all, group_audit) {
  d <- data.table::as.data.table(data.table::copy(prefusion_all))
  if (!nrow(d)) return(data.table())

  d[, .source_event_id := fragment_source_id(query)]
  d[, .fragment_haplotype := fragment_haplotype(IF)]
  d[, .fragment_svlen := abs(frag_safe_num(svlen))]

  z <- d[
    !is.na(CHROM) & !is.na(POS) & !is.na(.fragment_svlen) & .fragment_svlen > 0 &
      !is.na(.source_event_id) & nzchar(.source_event_id)
  ]
  if (!nrow(z)) return(data.table())

  src <- z[
    ,
    {
      hp <- unique(as.character(.fragment_haplotype))
      hp <- hp[!is.na(hp) & nzchar(hp)]
      dd <- frag_safe_num(DIST_TO_SIMPLE_REPEAT_BP)
      dd <- dd[is.finite(dd)]
      list(
        N_SOURCE_COMPONENTS = .N,
        N_SOURCE_HAPLOTYPES = length(hp),
        SOURCE_POS_MIN = min(frag_safe_num(POS), na.rm = TRUE),
        SOURCE_POS_MAX = max(frag_safe_num(POS), na.rm = TRUE),
        SOURCE_SPAN_BP = max(frag_safe_num(POS), na.rm = TRUE) - min(frag_safe_num(POS), na.rm = TRUE),
        SOURCE_SVLEN_SUM = sum(.fragment_svlen, na.rm = TRUE),
        SOURCE_ANY_SIMPLE_REPEAT = any(SIMPLE_REPEAT_OVERLAP %in% TRUE, na.rm = TRUE),
        SOURCE_MIN_DIST_SIMPLE_REPEAT_BP = if (length(dd)) min(dd) else NA_real_,
        SOURCE_QUERY_VECTOR = paste(as.character(query)[order(frag_safe_num(POS))], collapse = ";")
      )
    },
    by = .(CHROM, .source_event_id)
  ]

  if (nrow(group_audit)) {
    gs <- group_audit[
      ,
      .(
        N_FRAGMENT_GROUPS = .N,
        N_MULTI_COMPONENT_GROUPS = sum(N_COMPONENTS >= 2L),
        N_GROUP_SOURCE_COMPLEX = sum(FRAGMENTATION_CLASS == "SOURCE_COMPLEX"),
        N_GROUP_COMPOUND_HIGHCONF = sum(FRAGMENTATION_CLASS == "COMPOUND_HIGHCONF")
      ),
      by = .(CHROM, .source_event_id)
    ]
    src <- merge(src, gs, by = c("CHROM", ".source_event_id"), all.x = TRUE)
  } else {
    src[, `:=`(
      N_FRAGMENT_GROUPS = 0L,
      N_MULTI_COMPONENT_GROUPS = 0L,
      N_GROUP_SOURCE_COMPLEX = 0L,
      N_GROUP_COMPOUND_HIGHCONF = 0L
    )]
  }

  for (nm in c(
    "N_FRAGMENT_GROUPS", "N_MULTI_COMPONENT_GROUPS",
    "N_GROUP_SOURCE_COMPLEX", "N_GROUP_COMPOUND_HIGHCONF"
  )) {
    src[is.na(get(nm)), (nm) := 0L]
  }

  src[, SOURCE_REPEAT_CONTEXT := (
    SOURCE_ANY_SIMPLE_REPEAT |
      (!is.na(SOURCE_MIN_DIST_SIMPLE_REPEAT_BP) &
         SOURCE_MIN_DIST_SIMPLE_REPEAT_BP <= FRAG$repeat_near_bp)
  )]

  # Source-level complexity is intentionally conservative.  A normal biallelic
  # event (one component per haplotype) is never complex merely because N=2.
  src[, SOURCE_COMPLEX_FLAG := (
    N_GROUP_SOURCE_COMPLEX > 0L |
      N_SOURCE_COMPONENTS >= FRAG$source_complex_min_components |
      (
        N_SOURCE_COMPONENTS >= FRAG$complex_min_components &
          SOURCE_SPAN_BP >= FRAG$source_complex_span_bp &
          SOURCE_REPEAT_CONTEXT
      )
  )]

  src[]
}

apply_source_complexity_to_groups <- function(group_audit, source_audit) {
  if (!nrow(group_audit) || !nrow(source_audit)) return(group_audit)
  x <- data.table::copy(group_audit)
  ss <- source_audit[, .(CHROM, .source_event_id, SOURCE_COMPLEX_FLAG)]
  x[ss, on = .(CHROM, .source_event_id), SOURCE_COMPLEX_GLOBAL := i.SOURCE_COMPLEX_FLAG]
  x[is.na(SOURCE_COMPLEX_GLOBAL), SOURCE_COMPLEX_GLOBAL := FALSE]

  # IMPORTANT v5: source-level complexity is an audit/caution flag only.
  # It must never poison an otherwise clean haplotype/group.  The old behaviour
  # promoted every subgroup of a globally complex Atlas/source to SOURCE_COMPLEX,
  # which deleted exact Q5 counterparts in the twin benchmark.
  x[]
}

attach_fragment_source_annotations <- function(prefusion_all, source_audit) {
  d <- data.table::as.data.table(data.table::copy(prefusion_all))
  if (!".source_event_id" %in% names(d)) d[, .source_event_id := fragment_source_id(query)]
  d[, SOURCE_COMPLEX_GLOBAL := FALSE]
  d[, SOURCE_N_COMPONENTS := NA_integer_]
  d[, SOURCE_SPAN_BP := NA_real_]

  if (!nrow(source_audit)) return(as.data.frame(d))
  ss <- source_audit[, .(
    CHROM, .source_event_id,
    SOURCE_COMPLEX_FLAG,
    N_SOURCE_COMPONENTS,
    SOURCE_SPAN_BP
  )]
  d[ss, on = .(CHROM, .source_event_id), `:=`(
    SOURCE_COMPLEX_GLOBAL = i.SOURCE_COMPLEX_FLAG,
    SOURCE_N_COMPONENTS = i.N_SOURCE_COMPONENTS,
    SOURCE_SPAN_BP = i.SOURCE_SPAN_BP
  )]
  d[is.na(SOURCE_COMPLEX_GLOBAL), SOURCE_COMPLEX_GLOBAL := FALSE]
  as.data.frame(d)
}

attach_fragment_group_annotations <- function(prefusion_all, group_audit) {
  d <- data.table::as.data.table(data.table::copy(prefusion_all))
  d[, .source_event_id := fragment_source_id(query)]
  d[, .fragment_haplotype := fragment_haplotype(IF)]
  d[, .fragment_group_id := paste(CHROM, .source_event_id, .fragment_haplotype, sep = "||")]
  d[is.na(.source_event_id) | is.na(.fragment_haplotype), .fragment_group_id := NA_character_]

  if (!nrow(group_audit)) {
    # Singletons and unknown haplotypes still need a typed evidence schema.
    d[, FRAGMENT_GROUP_N_COMPONENTS := 1L]
    return(as.data.frame(d))
  }

  cols <- c(
    "FRAGMENT_GROUP_ID", "FRAGMENTATION_RISK", "FRAGMENTATION_CLASS", "FRAGMENT_Q",
    "N_COMPONENTS", "COMPOUND_SVLEN", "COMPONENT_SPAN_BP", "MIN_COMPONENT_GAP_BP",
    "MAX_COMPONENT_GAP_BP", "ANY_SIMPLE_REPEAT", "REPEAT_CONTEXT", "MIN_DIST_SIMPLE_REPEAT_BP",
    "BEST_VCF_ALT_SIMILARITY", "BEST_VCF_ALT_DPOS", "BEST_VCF_ALT_SVLEN",
    "BEST_VCF_ALT_CALLER", "N_VCF_ALT_CALLERS_GE_MIN", "N_VCF_ALT_CALLERS_GE_STRONG",
    "RESOLVED_SVLEN_SUM", "ANCHOR_MAX_SVLEN", "ANCHOR_FRACTION", "MIN_COMPONENT_LR",
    "FRAGMENT_AUTO_CANDIDATE", "ANCHORED_DROPOUT_CANDIDATE", "FRAGMENT_RECONSTRUCT_CANDIDATE"
  )
  ga <- group_audit[, ..cols]
  data.table::setnames(ga, "FRAGMENT_GROUP_ID", ".fragment_group_id")
  d[ga, on = ".fragment_group_id", `:=`(
    FRAGMENTATION_RISK = i.FRAGMENTATION_RISK,
    FRAGMENTATION_CLASS = i.FRAGMENTATION_CLASS,
    FRAGMENT_Q = i.FRAGMENT_Q,
    FRAGMENT_GROUP_N_COMPONENTS = i.N_COMPONENTS,
    FRAGMENT_COMPOUND_SVLEN = i.COMPOUND_SVLEN,
    FRAGMENT_COMPONENT_SPAN_BP = i.COMPONENT_SPAN_BP,
    FRAGMENT_MIN_GAP_BP = i.MIN_COMPONENT_GAP_BP,
    FRAGMENT_MAX_GAP_BP = i.MAX_COMPONENT_GAP_BP,
    FRAGMENT_ANY_SIMPLE_REPEAT = i.ANY_SIMPLE_REPEAT,
    FRAGMENT_REPEAT_CONTEXT = i.REPEAT_CONTEXT,
    FRAGMENT_MIN_DIST_SIMPLE_REPEAT_BP = i.MIN_DIST_SIMPLE_REPEAT_BP,
    FRAGMENT_BEST_VCF_ALT_SIMILARITY = i.BEST_VCF_ALT_SIMILARITY,
    FRAGMENT_BEST_VCF_ALT_DPOS = i.BEST_VCF_ALT_DPOS,
    FRAGMENT_BEST_VCF_ALT_SVLEN = i.BEST_VCF_ALT_SVLEN,
    FRAGMENT_BEST_VCF_ALT_CALLER = i.BEST_VCF_ALT_CALLER,
    FRAGMENT_N_VCF_ALT_CALLERS_GE_MIN = i.N_VCF_ALT_CALLERS_GE_MIN,
    FRAGMENT_N_VCF_ALT_CALLERS_GE_STRONG = i.N_VCF_ALT_CALLERS_GE_STRONG,
    FRAGMENT_RESOLVED_SVLEN_SUM = i.RESOLVED_SVLEN_SUM,
    FRAGMENT_ANCHOR_MAX_SVLEN = i.ANCHOR_MAX_SVLEN,
    FRAGMENT_ANCHOR_FRACTION = i.ANCHOR_FRACTION,
    FRAGMENT_MIN_COMPONENT_LR = i.MIN_COMPONENT_LR,
    FRAGMENT_AUTO_CANDIDATE = i.FRAGMENT_AUTO_CANDIDATE,
    FRAGMENT_ANCHORED_DROPOUT_CANDIDATE = i.ANCHORED_DROPOUT_CANDIDATE,
    FRAGMENT_RECONSTRUCT_CANDIDATE = i.FRAGMENT_RECONSTRUCT_CANDIDATE
  )]
  as.data.frame(d)
}

collapse_reconstructable_fragment_groups <- function(prefusion_all, group_audit) {
  if (FRAG$mode != "highconf" || !nrow(prefusion_all) || !nrow(group_audit)) return(prefusion_all)

  targets <- group_audit[FRAGMENT_RECONSTRUCT_CANDIDATE %in% TRUE, FRAGMENT_GROUP_ID]
  if (!length(targets)) return(prefusion_all)

  d <- data.table::as.data.table(data.table::copy(prefusion_all))
  if (!".fragment_group_id" %in% names(d)) return(prefusion_all)

  synthetic <- list()
  remove_rows <- integer()

  for (gid in targets) {
    idx <- which(d$.fragment_group_id == gid)
    if (length(idx) < 2L) next
    g <- d[idx]
    meta <- group_audit[FRAGMENT_GROUP_ID == gid][1L]
    if (!nrow(meta)) next

    resolved_mask <- frag_resolved(g$event_class) & frag_resolved(g$event)
    cls <- unique(as.character(g$event_class[resolved_mask]))
    cls <- cls[frag_resolved(cls)]
    evt <- unique(as.character(g$event[resolved_mask]))
    evt <- evt[frag_resolved(evt)]
    hap <- unique(as.character(g$.fragment_haplotype))
    hap <- hap[!is.na(hap) & nzchar(hap)]

    anchored <- isTRUE(meta$ANCHORED_DROPOUT_CANDIDATE[[1L]])
    standard <- isTRUE(meta$FRAGMENT_AUTO_CANDIDATE[[1L]])

    # Identity must be anchored by exactly one resolved class.  In the special
    # dropout route, resolved event must also be unique; unresolved components
    # contribute size/support only and never annotation identity.
    if (length(cls) != 1L || length(hap) != 1L) next
    if (anchored && length(evt) != 1L) next
    if (standard && any(!resolved_mask)) next

    rep_pool <- which(resolved_mask)
    if (!length(rep_pool)) next
    if ("QUAL" %in% names(g)) {
      qq <- frag_safe_num(g$QUAL[rep_pool])
      if (all(is.na(qq))) {
        ss <- abs(frag_safe_num(g$svlen[rep_pool]))
        rep_i <- rep_pool[which.max(replace(ss, is.na(ss), -Inf))[1L]]
      } else {
        rep_i <- rep_pool[which.max(replace(qq, is.na(qq), -Inf))[1L]]
      }
    } else {
      ss <- abs(frag_safe_num(g$svlen[rep_pool]))
      rep_i <- rep_pool[which.max(replace(ss, is.na(ss), -Inf))[1L]]
    }
    rep <- data.table::copy(g[rep_i])

    component_pos <- frag_safe_num(g$POS)
    component_sv <- abs(frag_safe_num(g$svlen))
    ord <- order(component_pos)

    rep[, POS := stats::median(component_pos, na.rm = TRUE)]
    rep[, svlen := sum(component_sv, na.rm = TRUE)]
    if ("query_length" %in% names(rep)) rep[, query_length := as.numeric(svlen)]
    rep[, event_class := cls[[1L]]]
    rep[, event := if (length(evt)) frag_mode_chr(evt) else frag_mode_chr(g$event)]
    rep[, completeness := frag_mode_chr(g$completeness[resolved_mask])]
    rep[, query := paste0(g$.source_event_id[1L], ";COMPOUND_", gsub(",", "_", hap[[1L]]))]

    # A synthetic compound has no literal assembled sequence; do not borrow one
    # component sequence for downstream TSD/BLAST.
    if ("ALT" %in% names(rep)) rep[, ALT := NA_character_]
    if ("contig" %in% names(rep)) rep[, contig := NA_character_]

    if ("SIMPLE_REPEAT_OVERLAP" %in% names(rep)) {
      rep[, SIMPLE_REPEAT_OVERLAP := any(g$SIMPLE_REPEAT_OVERLAP %in% TRUE, na.rm = TRUE)]
    }
    if ("DIST_TO_SIMPLE_REPEAT_BP" %in% names(rep)) {
      dd <- frag_safe_num(g$DIST_TO_SIMPLE_REPEAT_BP)
      dd <- dd[is.finite(dd)]
      rep[, DIST_TO_SIMPLE_REPEAT_BP := if (length(dd)) min(dd) else NA_real_]
    }
    if ("polyA" %in% names(rep)) {
      rep[, polyA := if (any(toupper(as.character(g$polyA)) == "TRUE", na.rm = TRUE)) "TRUE" else "FALSE"]
    }
    if ("polyA_length" %in% names(rep)) {
      pl <- frag_safe_num(g$polyA_length)
      rep[, polyA_length := if (any(is.finite(pl))) max(pl, na.rm = TRUE) else 0]
    }

    rep[, FRAGMENT_COMPOUND_SYNTHETIC := TRUE]
    rep[, FRAGMENT_RECONSTRUCTION_MODE := if (anchored) "ANCHORED_DROPOUT" else "HIGHCONF_RESOLVED"]
    rep[, FRAGMENT_RECONSTRUCTION_CAUTION := anchored]
    rep[, FRAGMENT_RECONSTRUCTION_ANCHOR_CLASS := cls[[1L]]]
    rep[, FRAGMENT_RECONSTRUCTION_VCF_SIM := meta$BEST_VCF_ALT_SIMILARITY[[1L]]]
    rep[, FRAGMENT_RECONSTRUCTION_VCF_CALLERS := meta$N_VCF_ALT_CALLERS_GE_STRONG[[1L]]]
    rep[, FRAGMENT_RECONSTRUCTION_MIN_LR := meta$MIN_COMPONENT_LR[[1L]]]
    rep[, FRAGMENT_RECONSTRUCTION_ANCHOR_FRACTION := meta$ANCHOR_FRACTION[[1L]]]
    rep[, FRAGMENT_COMPONENT_QUERIES := paste(as.character(g$query)[ord], collapse = ";")]
    rep[, FRAGMENT_COMPONENT_POS_VECTOR := paste(component_pos[ord], collapse = ";")]
    rep[, FRAGMENT_COMPONENT_SVLEN_VECTOR := paste(component_sv[ord], collapse = ";")]
    rep[, FRAGMENT_ORIGINAL_N := nrow(g)]
    rep[, FRAGMENTATION_CLASS := if (anchored) "COMPOUND_RECONSTRUCTED_ANCHORED" else "COMPOUND_RECONSTRUCTED_HIGHCONF"]
    rep[, FRAGMENT_Q := if (anchored) 4L else meta$FRAGMENT_Q[[1L]]]

    synthetic[[length(synthetic) + 1L]] <- rep
    remove_rows <- c(remove_rows, idx)
  }

  if (!length(synthetic)) return(prefusion_all)

  kept <- d[-unique(remove_rows)]
  out <- data.table::rbindlist(
    c(list(kept), synthetic),
    fill = TRUE,
    use.names = TRUE,
    ignore.attr = TRUE
  )
  data.table::setorder(out, CHROM, POS, svlen)
  message(
    "Fragment compounds reconstructed: ", length(synthetic),
    " | anchored_dropout=", sum(vapply(synthetic, function(x) isTRUE(x$FRAGMENT_RECONSTRUCTION_CAUTION[[1L]]), logical(1)))
  )
  as.data.frame(out)
}


# Coverage was annotated and filtered before poly(A/T) and reconstruction.
prefusion_all <- as.data.frame(FinalTable)

# ------------------------------------------------------------------------------
# Fragmentation audit: repeat context + same-source/same-haplotype components +
# optional VCF corroboration of the SUM of components.
# ------------------------------------------------------------------------------
if (FRAG$mode == "off") {
  simple_repeats <- data.table()
  fragment_vcf <- data.table()
  prefusion_all <- annotate_simple_repeats(prefusion_all, simple_repeats)
  fragment_group_audit <- data.table()
  fragment_source_audit <- data.table()
  prefusion_all <- attach_fragment_group_annotations(prefusion_all, fragment_group_audit)
  prefusion_all <- attach_fragment_source_annotations(prefusion_all, fragment_source_audit)
} else {
  simple_repeats <- load_simple_repeats(simple_repeats_bed)
  prefusion_all <- annotate_simple_repeats(prefusion_all, simple_repeats)
  fragment_vcf <- load_fragment_vcfs()
  message("Fragmentation VCF INS records: ", nrow(fragment_vcf))
  fragment_group_audit <- build_fragment_group_audit(prefusion_all, fragment_vcf)
  fragment_source_audit <- build_fragment_source_audit(prefusion_all, fragment_group_audit)
  fragment_group_audit <- apply_source_complexity_to_groups(fragment_group_audit, fragment_source_audit)
  prefusion_all <- attach_fragment_group_annotations(prefusion_all, fragment_group_audit)
  prefusion_all <- attach_fragment_source_annotations(prefusion_all, fragment_source_audit)
}

message(
  "Fragment groups: ", nrow(fragment_group_audit),
  if (nrow(fragment_group_audit)) paste0(
    " | multi-component=", sum(fragment_group_audit$N_COMPONENTS >= 2L),
    " | high-risk=", sum(fragment_group_audit$FRAGMENTATION_RISK == "HIGH"),
    " | highconf-compound=", sum(fragment_group_audit$FRAGMENT_AUTO_CANDIDATE %in% TRUE),
    " | source-complex-groups=", sum(fragment_group_audit$FRAGMENTATION_CLASS == "SOURCE_COMPLEX"),
    " | source-complex-events=", if (nrow(fragment_source_audit)) sum(fragment_source_audit$SOURCE_COMPLEX_FLAG %in% TRUE) else 0L
  ) else ""
)

# Reconstruct before filtering unresolved classifier dropouts.  This is required
# for ANCHORED_DROPOUT: unresolved components can contribute size only after a
# resolved anchor + VCF + LR gates have established a conservative hypothesis.
prefusion_reconstructed <- collapse_reconstructable_fragment_groups(prefusion_all, fragment_group_audit)

"strong filter MEI"
prefusion <- prefusion_reconstructed %>%
  filter(
    event_class != "unresolved",
    event != "unresolved"
  )


if (!nrow(prefusion)) empty_result("No resolved mobile events")

"Processing - locus-aware aggregate intrasample merge v3 + fragmentation v4"

library(dplyr)
library(stringr)
library(tibble)
library(purrr)

# ============================================================
# FUSION INTRA-MUESTRA V3: LOCUS != REPRESENTACION/ALELO
#
# Principios:
#   1) Atlas_ID es provenance, nunca una frontera biologica.
#   2) Primero se decide si dos reconstrucciones pertenecen al
#      mismo LOCUS estructural.
#   3) El breakpoint exacto es la evidencia espacial mas fuerte.
#      Si POS es exactamente igual:
#        a) same class + HET en haplotipos distintos -> mismo locus,
#           aunque SVLEN/event/completeness difieran.
#        b) same class + same event -> mismo locus, aunque SVLEN
#           fluctue mucho; se interpreta como variabilidad de
#           reconstruccion/representacion del mismo breakpoint.
#   4) Regla heteroalelica validada manualmente:
#      dPOS <= 5 bp + same class + same event + HET en haplotipos
#      distintos -> mismo locus, sin usar SVLEN/completeness como veto.
#      El desplazamiento pequeno se interpreta como incertidumbre de
#      breakpoint/reconstruccion, no como un segundo locus.
#   5) Fuera de esa excepcion, SVLEN vuelve a ser gate y aumenta
#      progresivamente su exigencia con dPOS.
#   6) event/completeness son anotaciones descriptivas y pueden
#      fluctuar; event_class sigue siendo un gate fuerte.
#   7) Los conflictos de event_class NO se auto-fusionan en v3:
#      se auditan para evitar sobreagregacion no validada.
#   8) No hay chaining: complete-linkage exige que TODOS los pares
#      dentro de un agregado sean compatibles.
#   9) Fase/GT se interpretan DESPUES de definir el locus.
#  10) RECONSTRUCTION_QUALITY_SCORE (1-5) describe consistencia de
#      la reconstruccion agregada. NO es calidad de presencia del locus.
# ============================================================


# ------------------------------------------------------------
# Parametros del clasificador
# ------------------------------------------------------------
AUDIT_DPOS_MAX <- CFG$audit_dpos_max

# POS exacta: regla especial de locus.
EXACT_LOCUS_DPOS <- CFG$exact_locus_dpos

# Heteroalelia de locus validada por inspeccion IGV:
# mismo tipo/evento en ambos haplotipos, breakpoint desplazado <=5 bp.
# En esta regla SVLEN y completeness son descriptivos, no gates.
HETEROALLELIC_NEAR_DPOS_MAX <- CFG$heteroallelic_near_dpos_max

# Fuera de POS exacta / excepcion heteroalelica:
# identidad estructural tipo same-allele.
VERY_NEAR_DPOS_MAX <- CFG$very_near_dpos_max
VERY_NEAR_SIZE_MIN <- CFG$very_near_size_min

NEAR_DPOS_MAX <- CFG$near_dpos_max
NEAR_SIZE_MIN <- CFG$near_size_min

MID_DPOS_MAX <- CFG$mid_dpos_max
MID_SIZE_MIN <- CFG$mid_size_min

FAR_DPOS_MAX <- CFG$far_dpos_max
FAR_SIZE_MIN <- CFG$far_size_min

# Event puede fluctuar, pero fuera de la regla heteroalelica exacta
# exigimos geometria/tamano fuertes.
EXACT_EVENT_FLEX_SIZE_MIN <- CFG$exact_event_flex_size_min
VERY_NEAR_EVENT_FLEX_DPOS_MAX <- CFG$very_near_event_flex_dpos_max
VERY_NEAR_EVENT_FLEX_SIZE_MIN <- CFG$very_near_event_flex_size_min
NEAR_EVENT_FLEX_DPOS_MAX <- CFG$near_event_flex_dpos_max
NEAR_EVENT_FLEX_SIZE_MIN <- CFG$near_event_flex_size_min

# event_class conflict: solo auditoria en v3.
CLASS_CONFLICT_AUDIT_DPOS_MAX <- CFG$class_conflict_audit_dpos_max
CLASS_CONFLICT_AUDIT_SIZE_MIN <- CFG$class_conflict_audit_size_min

# Descriptivo: NO decide fusion.
SVLEN_STABLE_MIN <- CFG$svlen_stable_min
SVLEN_MODERATE_MIN <- CFG$svlen_moderate_min


# ============================================================
# Utilidades
# ============================================================

safe_num <- function(x) {
  suppressWarnings(as.numeric(x))
}

size_similarity <- function(a, b) {
  a <- abs(safe_num(a))
  b <- abs(safe_num(b))

  if (
    length(a) != 1L ||
    length(b) != 1L ||
    is.na(a) ||
    is.na(b) ||
    a <= 0 ||
    b <= 0
  ) {
    return(NA_real_)
  }

  min(a, b) / max(a, b)
}

extract_haplotype <- function(x) {
  x <- as.character(x)
  out <- stringr::str_extract(x, "^[^,]+,[^,]+")
  out[is.na(x) | !nzchar(x)] <- NA_character_
  out
}

is_het_gt <- function(x) {
  as.character(x) %in% c("0/1", "1/0", "0|1", "1|0")
}

is_hom_alt_gt <- function(x) {
  as.character(x) %in% c("1/1", "1|1")
}

normalize_gt <- function(x) {
  x <- as.character(x)
  dplyr::case_when(
    x %in% c("1/1", "1|1") ~ "1/1",
    x %in% c("0/1", "1/0", "0|1", "1|0") ~ "0/1",
    x %in% c("0/0", "0|0") ~ "0/0",
    TRUE ~ x
  )
}

collapse_set <- function(x) {
  x <- unique(as.character(x))
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x) == 0) return(NA_character_)
  paste(sort(x), collapse = ";")
}

collapse_counts <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x) == 0) return(NA_character_)

  tt <- sort(table(x), decreasing = TRUE)
  paste0(names(tt), ":", as.integer(tt), collapse = ";")
}

safe_median <- function(x) {
  x <- safe_num(x)
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  stats::median(x)
}

safe_iqr <- function(x) {
  x <- safe_num(x)
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  stats::IQR(x)
}

safe_min <- function(x) {
  x <- safe_num(x)
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  min(x)
}

safe_max <- function(x) {
  x <- safe_num(x)
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  max(x)
}

safe_max_qual_index <- function(x) {
  q <- safe_num(x)
  if (length(q) == 0) return(1L)
  if (all(is.na(q))) return(1L)
  which.max(replace(q, is.na(q), -Inf))[1]
}

representative_position <- function(x) {
  x <- safe_num(x)
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)

  tt <- table(x)
  top <- as.numeric(names(tt)[tt == max(tt)])

  if (length(top) == 1L) return(top[1])

  med <- stats::median(x)
  top[which.min(abs(top - med))][1]
}

minmax_similarity <- function(x) {
  x <- abs(safe_num(x))
  x <- x[is.finite(x) & x > 0]
  if (length(x) == 0) return(NA_real_)
  if (length(x) == 1L) return(1)
  min(x) / max(x)
}

collapse_svlen_by_haplotype <- function(haplotype, svlen) {
  h <- as.character(haplotype)
  s <- safe_num(svlen)
  ok <- !is.na(h) & nzchar(h) & is.finite(s)

  if (!any(ok)) return(NA_character_)

  d <- tibble(h = h[ok], s = s[ok]) %>%
    group_by(h) %>%
    summarise(
      med = stats::median(s),
      min = min(s),
      max = max(s),
      n = n(),
      .groups = "drop"
    ) %>%
    arrange(h)

  paste0(
    d$h, "=", d$med,
    "[", d$min, ",", d$max, ";n=", d$n, "]",
    collapse = ";"
  )
}


# ============================================================
# Clasificador modular de pares
#
# mergeable = mismo LOCUS, no necesariamente mismo alelo ALT.
# relation diferencia identidad alelica de heteroalelia estructural.
# ============================================================

classify_locus_pair <- function(a, b) {

  chr_a <- as.character(a$CHROM[1])
  chr_b <- as.character(b$CHROM[1])

  pos_a <- safe_num(a$POS[1])
  pos_b <- safe_num(b$POS[1])

  size_a <- safe_num(a$svlen[1])
  size_b <- safe_num(b$svlen[1])

  cls_a <- as.character(a$event_class[1])
  cls_b <- as.character(b$event_class[1])

  evt_a <- as.character(a$event[1])
  evt_b <- as.character(b$event[1])

  cmp_a <- as.character(a$completeness[1])
  cmp_b <- as.character(b$completeness[1])

  hap_a <- as.character(a$.haplotype[1])
  hap_b <- as.character(b$.haplotype[1])

  same_chr <-
    !is.na(chr_a) && !is.na(chr_b) && chr_a == chr_b

  dpos <- if (
    same_chr &&
    !is.na(pos_a) &&
    !is.na(pos_b)
  ) {
    abs(pos_a - pos_b)
  } else {
    Inf
  }

  sim <- size_similarity(size_a, size_b)

  same_class <-
    !is.na(cls_a) && !is.na(cls_b) && cls_a == cls_b

  same_event <-
    !is.na(evt_a) && !is.na(evt_b) && evt_a == evt_b

  same_completeness <-
    !is.na(cmp_a) && !is.na(cmp_b) && cmp_a == cmp_b

  different_haplotype <-
    !is.na(hap_a) && !is.na(hap_b) && hap_a != hap_b

  same_haplotype <-
    !is.na(hap_a) && !is.na(hap_b) && hap_a == hap_b

  both_het <- isTRUE(a$.is_het[1]) && isTRUE(b$.is_het[1])

  relation <- "SEPARATE"
  rule <- "NO_RULE"
  mergeable <- FALSE
  same_allele_like <- FALSE
  annotation_flex <- FALSE
  rule_rank <- 99L

  if (!same_chr) {
    rule <- "DIFFERENT_CHROM"

  } else if (!is.finite(dpos) || is.na(sim)) {
    rule <- "MISSING_GEOMETRY"

  # ----------------------------------------------------------
  # POS EXACTA + HAPLOTIPOS DISTINTOS:
  # misma clase define el mismo locus aunque los dos ALT tengan
  # longitudes/event/completeness distintos.
  # ----------------------------------------------------------
  } else if (
    dpos == EXACT_LOCUS_DPOS &&
    same_class &&
    different_haplotype &&
    both_het
  ) {
    relation <- "SAME_LOCUS_HETEROALLELIC"
    rule <- "EXACT_POS_HETEROALLELIC_SIZE_FREE"
    mergeable <- TRUE
    same_allele_like <- same_event && sim >= VERY_NEAR_SIZE_MIN
    annotation_flex <- !same_event || !same_completeness
    rule_rank <- 1L

  # ----------------------------------------------------------
  # POS EXACTA + misma clase/event:
  # mismo locus aunque SVLEN fluctue mucho. Si es same-phase,
  # se trata como reconstruccion/representacion variable.
  # ----------------------------------------------------------
  } else if (
    dpos == EXACT_LOCUS_DPOS &&
    same_class &&
    same_event
  ) {
    relation <- if (same_haplotype) {
      "SAME_LOCUS_RECONSTRUCTION_VARIABLE"
    } else {
      "SAME_LOCUS_EXACT_CLASS_EVENT"
    }
    rule <- "EXACT_POS_CLASS_EVENT_SIZE_FREE"
    mergeable <- TRUE
    same_allele_like <- sim >= VERY_NEAR_SIZE_MIN
    annotation_flex <- !same_completeness
    rule_rank <- 2L

  # ----------------------------------------------------------
  # POS EXACTA, misma clase, cambia event:
  # solo si los tamanos siguen siendo muy coherentes.
  # ----------------------------------------------------------
  } else if (
    dpos == EXACT_LOCUS_DPOS &&
    same_class &&
    !same_event &&
    sim >= EXACT_EVENT_FLEX_SIZE_MIN
  ) {
    relation <- "SAME_LOCUS_ANNOTATION_FLEX"
    rule <- "EXACT_POS_EVENT_FLEX"
    mergeable <- TRUE
    same_allele_like <- TRUE
    annotation_flex <- TRUE
    rule_rank <- 3L

  # ----------------------------------------------------------
  # BREAKPOINT CASI EXACTO + DOS HAPLOTIPOS:
  # regla validada manualmente en LSA0040 y contrastada en 300-buffy.
  # Si ambos HET representan el mismo tipo/evento en haplotipos
  # diferentes y dPOS <=5 bp, se trata como un mismo locus aunque
  # SVLEN sea muy discordante. completeness puede fluctuar por
  # truncacion/reconstruccion diferencial entre alelos.
  # ----------------------------------------------------------
  } else if (
    dpos > EXACT_LOCUS_DPOS &&
    dpos <= HETEROALLELIC_NEAR_DPOS_MAX &&
    same_class &&
    same_event &&
    different_haplotype &&
    both_het
  ) {
    relation <- "SAME_LOCUS_HETEROALLELIC_NEAR"
    rule <- "NEAR_POS_HETEROALLELIC_SIZE_FREE"
    mergeable <- TRUE
    same_allele_like <- !is.na(sim) && sim >= NEAR_SIZE_MIN
    annotation_flex <- !same_completeness
    rule_rank <- 4L

  # ----------------------------------------------------------
  # Conflicto de clase: NO auto-merge. Solo auditoria.
  # ----------------------------------------------------------
  } else if (
    dpos <= CLASS_CONFLICT_AUDIT_DPOS_MAX &&
    !same_class &&
    sim >= CLASS_CONFLICT_AUDIT_SIZE_MIN
  ) {
    relation <- "EXACT_CLASS_CONFLICT_AUDIT"
    rule <- "CLASS_CONFLICT_REQUIRES_REVIEW"

  # ----------------------------------------------------------
  # A partir de aqui POS ya no es exacta: SVLEN vuelve a ser gate.
  # ----------------------------------------------------------
  } else if (
    dpos <= VERY_NEAR_DPOS_MAX &&
    same_class &&
    same_event &&
    sim >= VERY_NEAR_SIZE_MIN
  ) {
    relation <- "SAME_ALLELE_LIKE"
    rule <- "VERY_NEAR_CLASS_EVENT"
    mergeable <- TRUE
    same_allele_like <- TRUE
    rule_rank <- 5L

  } else if (
    dpos <= VERY_NEAR_EVENT_FLEX_DPOS_MAX &&
    same_class &&
    !same_event &&
    same_completeness &&
    sim >= VERY_NEAR_EVENT_FLEX_SIZE_MIN
  ) {
    relation <- "SAME_ALLELE_ANNOTATION_FLEX"
    rule <- "VERY_NEAR_EVENT_FLEX"
    mergeable <- TRUE
    same_allele_like <- TRUE
    annotation_flex <- TRUE
    rule_rank <- 6L

  } else if (
    dpos <= NEAR_DPOS_MAX &&
    same_class &&
    same_event &&
    sim >= NEAR_SIZE_MIN
  ) {
    relation <- "SAME_ALLELE_LIKE"
    rule <- "NEAR_CLASS_EVENT"
    mergeable <- TRUE
    same_allele_like <- TRUE
    rule_rank <- 7L

  } else if (
    dpos <= NEAR_EVENT_FLEX_DPOS_MAX &&
    same_class &&
    !same_event &&
    same_completeness &&
    sim >= NEAR_EVENT_FLEX_SIZE_MIN
  ) {
    relation <- "SAME_ALLELE_ANNOTATION_FLEX"
    rule <- "NEAR_EVENT_FLEX"
    mergeable <- TRUE
    same_allele_like <- TRUE
    annotation_flex <- TRUE
    rule_rank <- 8L

  } else if (
    dpos <= MID_DPOS_MAX &&
    same_class &&
    same_event &&
    sim >= MID_SIZE_MIN
  ) {
    relation <- "SAME_ALLELE_LIKE"
    rule <- "MID_CLASS_EVENT"
    mergeable <- TRUE
    same_allele_like <- TRUE
    rule_rank <- 9L

  } else if (
    dpos <= FAR_DPOS_MAX &&
    same_class &&
    same_event &&
    sim >= FAR_SIZE_MIN
  ) {
    relation <- "SAME_ALLELE_LIKE"
    rule <- "STRICT_100BP_CLASS_EVENT"
    mergeable <- TRUE
    same_allele_like <- TRUE
    rule_rank <- 10L

  } else if (dpos <= AUDIT_DPOS_MAX) {
    relation <- "LOCAL_NONMATCH"
    rule <- "LOCAL_BUT_NOT_COMPATIBLE"
  }

  # Complete-linkage: compatible <1, incompatible >1.
  merge_distance <- if (mergeable) {
    base <- rule_rank / 10
    geom <- min(dpos / max(FAR_DPOS_MAX, 1), 1) * 0.03
    size_penalty <- if (is.na(sim)) 0.03 else (1 - sim) * 0.03
    min(base + geom + size_penalty, 0.95)
  } else {
    2
  }

  tibble(
    same_chr = same_chr,
    dPOS = dpos,
    size_similarity = sim,
    same_class = same_class,
    same_event = same_event,
    same_completeness = same_completeness,
    same_haplotype = same_haplotype,
    different_haplotype = different_haplotype,
    both_het = both_het,
    relation = relation,
    merge_rule = rule,
    mergeable = mergeable,
    same_allele_like = same_allele_like,
    annotation_flex = annotation_flex,
    rule_rank = rule_rank,
    merge_distance = merge_distance
  )
}


# ============================================================
# Preparar filas de trabajo
# ============================================================

prefusion_work <- prefusion %>%
  mutate(
    .merge_rowid = seq_len(n()),
    .haplotype = extract_haplotype(IF),
    .GT_norm = normalize_gt(GT),
    .is_het = is_het_gt(GT),
    .is_hom_alt = is_hom_alt_gt(GT)
  )


# ============================================================
# Pares locales + complete-linkage por cromosoma
# ============================================================

pair_audit_list <- list()
cluster_assignments <- integer(nrow(prefusion_work))
next_cluster_id <- 0L

chr_levels <- unique(as.character(prefusion_work$CHROM))
chr_levels <- chr_levels[!is.na(chr_levels)]

for (chr_now in chr_levels) {

  idx <- which(as.character(prefusion_work$CHROM) == chr_now)
  sub <- prefusion_work[idx, , drop = FALSE]
  nsub <- nrow(sub)

  if (nsub == 0) next

  if (nsub == 1) {
    next_cluster_id <- next_cluster_id + 1L
    cluster_assignments[idx] <- next_cluster_id
    next
  }

  D <- matrix(2, nrow = nsub, ncol = nsub)
  diag(D) <- 0

  audit_chr <- list()
  audit_k <- 0L

  ord <- order(safe_num(sub$POS), na.last = TRUE)

  for (oo_i in seq_len(nsub - 1L)) {

    i <- ord[oo_i]
    pos_i <- safe_num(sub$POS[i])
    if (is.na(pos_i)) next

    for (oo_j in (oo_i + 1L):nsub) {

      j <- ord[oo_j]
      pos_j <- safe_num(sub$POS[j])
      if (is.na(pos_j)) next

      if ((pos_j - pos_i) > AUDIT_DPOS_MAX) break

      rel <- classify_locus_pair(
        sub[i, , drop = FALSE],
        sub[j, , drop = FALSE]
      )

      D[i, j] <- rel$merge_distance
      D[j, i] <- rel$merge_distance

      audit_k <- audit_k + 1L
      audit_chr[[audit_k]] <- tibble(
        rowid_A = sub$.merge_rowid[i],
        rowid_B = sub$.merge_rowid[j],
        source_ID_A = as.character(sub$ID[i]),
        source_ID_B = as.character(sub$ID[j]),
        query_A = if ("query" %in% names(sub)) as.character(sub$query[i]) else NA_character_,
        query_B = if ("query" %in% names(sub)) as.character(sub$query[j]) else NA_character_,
        CHROM = chr_now,
        POS_A = safe_num(sub$POS[i]),
        POS_B = safe_num(sub$POS[j]),
        SVLEN_A = safe_num(sub$svlen[i]),
        SVLEN_B = safe_num(sub$svlen[j]),
        class_A = as.character(sub$event_class[i]),
        class_B = as.character(sub$event_class[j]),
        event_A = as.character(sub$event[i]),
        event_B = as.character(sub$event[j]),
        completeness_A = as.character(sub$completeness[i]),
        completeness_B = as.character(sub$completeness[j]),
        haplotype_A = as.character(sub$.haplotype[i]),
        haplotype_B = as.character(sub$.haplotype[j]),
        GT_A = as.character(sub$.GT_norm[i]),
        GT_B = as.character(sub$.GT_norm[j]),
        dPOS = rel$dPOS,
        size_similarity = rel$size_similarity,
        same_class = rel$same_class,
        same_event = rel$same_event,
        same_completeness = rel$same_completeness,
        same_haplotype = rel$same_haplotype,
        different_haplotype = rel$different_haplotype,
        both_het = rel$both_het,
        relation = rel$relation,
        merge_rule = rel$merge_rule,
        mergeable = rel$mergeable,
        same_allele_like = rel$same_allele_like,
        annotation_flex = rel$annotation_flex,
        merge_distance = rel$merge_distance
      )
    }
  }

  if (length(audit_chr) > 0) {
    pair_audit_list[[length(pair_audit_list) + 1L]] <- bind_rows(audit_chr)
  }

  # A~B y B~C no basta: complete-linkage exige A~C.
  hc <- hclust(as.dist(D), method = "complete")
  cl_local <- cutree(hc, h = 1)

  for (local_id in unique(cl_local)) {
    next_cluster_id <- next_cluster_id + 1L
    members <- which(cl_local == local_id)
    cluster_assignments[idx[members]] <- next_cluster_id
  }
}

# CHROM ausente: singleton.
missing_cluster <- which(cluster_assignments == 0L)
if (length(missing_cluster) > 0) {
  for (ii in missing_cluster) {
    next_cluster_id <- next_cluster_id + 1L
    cluster_assignments[ii] <- next_cluster_id
  }
}

prefusion_work$.structural_cluster <- cluster_assignments

pair_audit <- if (length(pair_audit_list) > 0) {
  bind_rows(pair_audit_list)
} else {
  tibble()
}


# ============================================================
# Calidad de reconstruccion (1-5)
#
# IMPORTANTE:
#   - mide CONSISTENCIA de la representacion/reconstruccion agregada;
#   - NO mide evidencia de presencia del locus;
#   - un score bajo NO invalida el locus si la geometria/fase lo soportan.
#
# Se usa un esquema conservador por ejes y el score final es el minimo:
#   POS:  5 exacta; 4 <=2 bp; 3 <=5; 2 <=25; 1 >25
#   SIZE: 5 >=.95; 4 >=.85; 3 >=.50; 2 >=.25; 1 <.25
#   ANNOT:5 consistente; 4 solo completeness; 3 event; 1 class conflict
#   PHASE:5 resuelta/singleton; 3 no resuelta
#
# Esto evita que una dimension excelente oculte una reconstruccion
# extremadamente inestable en otra dimension.
# ============================================================

reconstruction_position_grade <- function(pos_span) {
  if (is.na(pos_span)) return(3L)
  if (pos_span == 0) return(5L)
  if (pos_span <= 2) return(4L)
  if (pos_span <= 5) return(3L)
  if (pos_span <= 25) return(2L)
  1L
}

reconstruction_size_grade <- function(sim, n_rep) {
  if (n_rep <= 1L) return(5L)
  if (is.na(sim)) return(3L)
  if (sim >= 0.95) return(5L)
  if (sim >= 0.85) return(4L)
  if (sim >= 0.50) return(3L)
  if (sim >= 0.25) return(2L)
  1L
}

reconstruction_annotation_grade <- function(
    class_conflict,
    event_conflict,
    completeness_conflict
) {
  if (isTRUE(class_conflict)) return(1L)
  if (isTRUE(event_conflict)) return(3L)
  if (isTRUE(completeness_conflict)) return(4L)
  5L
}

reconstruction_phase_grade <- function(
    n_rep,
    n_haplotypes,
    any_hom_alt
) {
  if (n_rep <= 1L) return(5L)
  if (isTRUE(any_hom_alt)) return(5L)
  if (!is.na(n_haplotypes) && n_haplotypes >= 1L) return(5L)
  3L
}

quality_label <- function(score) {
  dplyr::case_when(
    score >= 5L ~ "VERY_HIGH",
    score == 4L ~ "HIGH",
    score == 3L ~ "MODERATE",
    score == 2L ~ "LOW",
    TRUE ~ "VERY_LOW"
  )
}

quality_reasons <- function(
    score,
    pos_grade,
    size_grade,
    annotation_grade,
    phase_grade,
    pos_span,
    sv_sim,
    class_conflict,
    event_conflict,
    completeness_conflict,
    n_rep
) {
  if (n_rep <= 1L) {
    return("SINGLE_PASS_RECONSTRUCTION")
  }

  reasons <- character(0)

  if (pos_grade == score) {
    reasons <- c(
      reasons,
      if (is.na(pos_span)) {
        "POSITION_UNKNOWN"
      } else if (pos_span == 0) {
        "POSITION_EXACT"
      } else {
        paste0("POSITION_SPREAD_", format(pos_span, trim = TRUE), "BP")
      }
    )
  }

  if (size_grade == score) {
    reasons <- c(
      reasons,
      if (is.na(sv_sim)) {
        "SVLEN_SIMILARITY_UNKNOWN"
      } else if (sv_sim < 0.25) {
        "EXTREME_SVLEN_VARIATION"
      } else if (sv_sim < 0.50) {
        "HIGH_SVLEN_VARIATION"
      } else if (sv_sim < 0.85) {
        "MODERATE_SVLEN_VARIATION"
      } else if (sv_sim < 0.95) {
        "MILD_SVLEN_VARIATION"
      } else {
        "SVLEN_STABLE"
      }
    )
  }

  if (annotation_grade == score) {
    if (isTRUE(class_conflict)) {
      reasons <- c(reasons, "CLASS_CONFLICT")
    } else if (isTRUE(event_conflict)) {
      reasons <- c(reasons, "EVENT_CONFLICT")
    } else if (isTRUE(completeness_conflict)) {
      reasons <- c(reasons, "COMPLETENESS_VARIATION")
    } else {
      reasons <- c(reasons, "ANNOTATION_CONSISTENT")
    }
  }

  if (phase_grade == score && phase_grade < 5L) {
    reasons <- c(reasons, "PHASE_UNRESOLVED")
  }

  reasons <- unique(reasons)
  if (length(reasons) == 0) reasons <- "CONSISTENT_RECONSTRUCTION"
  paste(reasons, collapse = ";")
}


# ============================================================
# Agregar UN locus
# ============================================================

aggregate_locus_group <- function(g, pair_audit_table) {

  stopifnot(nrow(g) >= 1L)

  rep_idx <- safe_max_qual_index(g$QUAL)
  rep <- g[rep_idx, , drop = FALSE]

  haplotypes <- unique(as.character(g$.haplotype))
  haplotypes <- haplotypes[!is.na(haplotypes) & nzchar(haplotypes)]
  n_haplotypes <- length(haplotypes)

  het_haplotypes <- unique(as.character(g$.haplotype[g$.is_het]))
  het_haplotypes <- het_haplotypes[
    !is.na(het_haplotypes) & nzchar(het_haplotypes)
  ]
  n_het_haplotypes <- length(het_haplotypes)

  gt_norm <- unique(as.character(g$.GT_norm))
  gt_norm <- gt_norm[!is.na(gt_norm) & nzchar(gt_norm)]

  any_hom_alt <- any(g$.is_hom_alt, na.rm = TRUE)
  any_het <- any(g$.is_het, na.rm = TRUE)

  GT_final <- if (any_hom_alt) {
    "1/1"
  } else if (n_het_haplotypes >= 2L) {
    # Compatibilidad con el pipeline bialelico:
    # dos haplotipos ALT -> dosage 2. GT_STATUS informa si son
    # formas estructurales concordantes o heteroalelicas.
    "1/1"
  } else if (any_het) {
    "0/1"
  } else if (length(gt_norm) == 1L) {
    gt_norm[1]
  } else if (length(gt_norm) > 1L) {
    as.character(rep$.GT_norm[1])
  } else {
    NA_character_
  }

  ALT_DOSAGE <- dplyr::case_when(
    GT_final == "1/1" ~ 2L,
    GT_final == "0/1" ~ 1L,
    GT_final == "0/0" ~ 0L,
    TRUE ~ NA_integer_
  )

  classes <- unique(as.character(g$event_class))
  classes <- classes[!is.na(classes) & nzchar(classes)]

  events <- unique(as.character(g$event))
  events <- events[!is.na(events) & nzchar(events)]

  completes <- unique(as.character(g$completeness))
  completes <- completes[!is.na(completes) & nzchar(completes)]

  class_conflict <- length(classes) > 1L
  event_conflict <- length(events) > 1L
  completeness_conflict <- length(completes) > 1L

  member_ids <- g$.merge_rowid
  merge_rules <- character(0)
  merge_relations <- character(0)

  if (nrow(g) > 1L && nrow(pair_audit_table) > 0) {
    pa <- pair_audit_table %>%
      filter(
        rowid_A %in% member_ids,
        rowid_B %in% member_ids,
        mergeable
      )

    if (nrow(pa) > 0) {
      merge_rules <- unique(as.character(pa$merge_rule))
      merge_relations <- unique(as.character(pa$relation))
    }
  }

  source_ids <- unique(as.character(g$ID))
  source_ids <- source_ids[!is.na(source_ids) & nzchar(source_ids)]

  sv_min <- safe_min(g$svlen)
  sv_med <- safe_median(g$svlen)
  sv_max <- safe_max(g$svlen)
  sv_sim <- minmax_similarity(g$svlen)

  pos_min <- safe_min(g$POS)
  pos_med <- safe_median(g$POS)
  pos_max <- safe_max(g$POS)
  pos_span <- if (
    is.finite(pos_min) &&
    is.finite(pos_max)
  ) {
    pos_max - pos_min
  } else {
    NA_real_
  }

  size_variation_status <- if (nrow(g) == 1L) {
    "SINGLE"
  } else if (!is.na(sv_sim) && sv_sim >= SVLEN_STABLE_MIN) {
    "STABLE"
  } else if (!is.na(sv_sim) && sv_sim >= SVLEN_MODERATE_MIN) {
    "MODERATE_VARIATION"
  } else {
    "HIGH_VARIATION"
  }

  # Estado locus/genotipo.
  locus_status <- if (nrow(g) == 1L) {
    "UNIQUE"
  } else if (n_het_haplotypes >= 2L || any_hom_alt) {
    "BIALLELIC_ALT_LOCUS"
  } else if (n_haplotypes == 1L) {
    "SAME_HAPLOTYPE_RECONSTRUCTION_AGGREGATE"
  } else {
    "LOCUS_PHASE_UNRESOLVED"
  }

  GT_STATUS <- if (any_hom_alt && nrow(g) == 1L) {
    "HOM_ALT_REPORTED"
  } else if (n_het_haplotypes >= 2L) {
    if (
      size_variation_status == "HIGH_VARIATION" ||
      event_conflict ||
      completeness_conflict
    ) {
      "ALT_ALT_HETEROALLELIC"
    } else {
      "ALT_ALT_CONCORDANT"
    }
  } else if (GT_final == "0/1" && nrow(g) > 1L) {
    if (size_variation_status == "HIGH_VARIATION") {
      "HET_RECONSTRUCTION_VARIABLE"
    } else {
      "HET_REPRESENTATION_DUPLICATE"
    }
  } else if (GT_final == "0/1") {
    "HET_SINGLE"
  } else if (GT_final == "1/1") {
    "HOM_ALT"
  } else {
    "UNRESOLVED"
  }

  annotation_status <- if (class_conflict) {
    "CLASS_CONFLICT"
  } else if (event_conflict && completeness_conflict) {
    "EVENT_AND_COMPLETENESS_FLEX"
  } else if (event_conflict) {
    "EVENT_FLEX"
  } else if (completeness_conflict) {
    "COMPLETENESS_FLEX"
  } else {
    "CONSISTENT"
  }

  # ----------------------------------------------------------
  # Calidad de RECONSTRUCCION, no de presencia del locus.
  # El score final queda limitado por el eje mas debil.
  # ----------------------------------------------------------
  rq_position <- reconstruction_position_grade(pos_span)
  rq_size <- reconstruction_size_grade(sv_sim, nrow(g))
  rq_annotation <- reconstruction_annotation_grade(
    class_conflict,
    event_conflict,
    completeness_conflict
  )
  rq_phase <- reconstruction_phase_grade(
    nrow(g),
    n_haplotypes,
    any_hom_alt
  )

  reconstruction_quality_score <- min(
    c(rq_position, rq_size, rq_annotation, rq_phase),
    na.rm = TRUE
  )

  reconstruction_quality_label <- quality_label(
    reconstruction_quality_score
  )

  reconstruction_quality_reason <- quality_reasons(
    score = reconstruction_quality_score,
    pos_grade = rq_position,
    size_grade = rq_size,
    annotation_grade = rq_annotation,
    phase_grade = rq_phase,
    pos_span = pos_span,
    sv_sim = sv_sim,
    class_conflict = class_conflict,
    event_conflict = event_conflict,
    completeness_conflict = completeness_conflict,
    n_rep = nrow(g)
  )

  out <- rep
  component_support <- primary_alt_support(g)
  out$ALT_SUPPORT_MAX <- if(any(is.finite(component_support))) max(component_support[is.finite(component_support)]) else NA_real_

  # POS final = breakpoint observado mas representado; en empate,
  # el mas cercano a la mediana. POS_MEDIAN conserva el resumen.
  out$POS <- representative_position(g$POS)
  out$POS_MEDIAN <- pos_med

  # SVLEN escalar para compatibilidad: mediana.
  # MIN/MEDIAN/MAX/VECTOR conservan la variacion alelica real.
  out$svlen <- sv_med
  out$SVLEN_MEDIAN <- sv_med
  out$GT <- GT_final

  out$GT_STATUS <- GT_STATUS
  out$ALT_DOSAGE <- ALT_DOSAGE
  out$LOCUS_STATUS <- locus_status
  out$SIZE_VARIATION_STATUS <- size_variation_status
  out$ANNOTATION_STATUS <- annotation_status

  out$RECONSTRUCTION_QUALITY_SCORE <- reconstruction_quality_score
  out$RECONSTRUCTION_QUALITY_LABEL <- reconstruction_quality_label
  out$RECONSTRUCTION_QUALITY_REASON <- reconstruction_quality_reason
  out$RQ_POSITION <- rq_position
  out$RQ_SIZE <- rq_size
  out$RQ_ANNOTATION <- rq_annotation
  out$RQ_PHASE <- rq_phase

  out$N_REPRESENTATIONS <- nrow(g)
  out$N_SOURCE_IDS <- length(source_ids)
  out$SOURCE_IDS <- if (length(source_ids) > 0) {
    paste(source_ids, collapse = ";")
  } else {
    NA_character_
  }

  out$N_HAPLOTYPES <- n_haplotypes
  out$N_ALT_HAPLOTYPES <- if (any_hom_alt) {
    2L
  } else if (n_het_haplotypes > 0L) {
    n_het_haplotypes
  } else {
    NA_integer_
  }

  out$HAPLOTYPE_SET <- if (n_haplotypes > 0) {
    paste(sort(haplotypes), collapse = ";")
  } else {
    NA_character_
  }

  out$SVLEN_BY_HAPLOTYPE <- collapse_svlen_by_haplotype(
    g$.haplotype,
    g$svlen
  )

  # Compatibilidad con versiones anteriores: numero de cromosomas ALT,
  # NO numero de reconstrucciones.
  out$STRUCTURAL_COPY_COUNT <- ALT_DOSAGE

  out$CLASS_SET <- collapse_set(g$event_class)
  out$CLASS_COUNTS <- collapse_counts(g$event_class)
  out$EVENT_SET <- collapse_set(g$event)
  out$EVENT_COUNTS <- collapse_counts(g$event)
  out$COMPLETENESS_SET <- collapse_set(g$completeness)
  out$COMPLETENESS_COUNTS <- collapse_counts(g$completeness)

  out$event_class_status <- if (class_conflict) "CONFLICT" else "CONSISTENT"
  out$event_status <- if (event_conflict) "CONFLICT" else "CONSISTENT"
  out$completeness_status <- if (completeness_conflict) "CONFLICT" else "CONSISTENT"

  out$MERGE_RULES <- if (length(merge_rules) > 0) {
    paste(sort(merge_rules), collapse = ";")
  } else {
    "NONE"
  }

  out$MERGE_RELATIONS <- if (length(merge_relations) > 0) {
    paste(sort(merge_relations), collapse = ";")
  } else {
    "NONE"
  }

  out$POS_MIN <- pos_min
  out$POS_MAX <- pos_max
  out$POS_SPAN <- pos_span
  out$POS_IQR <- safe_iqr(g$POS)

  out$SVLEN_MIN <- sv_min
  out$SVLEN_MAX <- sv_max
  out$SVLEN_IQR <- safe_iqr(g$svlen)
  out$SVLEN_MINMAX_SIMILARITY <- sv_sim

  # Provenance completa.
  out$Pos_vector <- list(safe_num(g$POS))
  out$Size_Vector <- list(safe_num(g$svlen))
  out$Class_vector <- list(as.character(g$event_class))
  out$Event_vector <- list(as.character(g$event))
  out$Completeness_vector <- list(as.character(g$completeness))
  out$IF_vector <- list(as.character(g$IF))
  out$Alele_vector <- list(as.character(g$.haplotype))
  out$GT_vector <- list(as.character(g$.GT_norm))
  out$Source_ID_vector <- list(as.character(g$ID))

  if ("query" %in% names(g)) {
    out$Query_vector <- list(as.character(g$query))
  } else {
    out$Query_vector <- list(character(0))
  }

  if ("ALT" %in% names(g)) {
    out$ALT_list <- list(as.character(g$ALT))
  } else {
    out$ALT_list <- list(character(0))
  }

  if ("polyA" %in% names(g)) {
    out$polyA_vector <- list(as.character(g$polyA))
  }

  if ("polyA_length" %in% names(g)) {
    out$polyA_length_vector <- list(safe_num(g$polyA_length))
  }

  # Fragmentation-aware provenance / audit summary.  These fields do not alter
  # locus identity; they make source fragmentation auditable downstream.
  if (".source_event_id" %in% names(g)) {
    src_evt <- unique(as.character(g$.source_event_id))
    src_evt <- src_evt[!is.na(src_evt) & nzchar(src_evt)]
    out$SOURCE_EVENT_IDS <- if (length(src_evt)) paste(sort(src_evt), collapse = ";") else NA_character_
    out$Source_Event_ID_vector <- list(as.character(g$.source_event_id))
  } else {
    out$SOURCE_EVENT_IDS <- NA_character_
    out$Source_Event_ID_vector <- list(character(0))
  }

  if (".fragment_group_id" %in% names(g)) {
    fg <- unique(as.character(g$.fragment_group_id))
    fg <- fg[!is.na(fg) & nzchar(fg)]
    out$FRAGMENT_GROUP_IDS <- if (length(fg)) paste(sort(fg), collapse = ";") else NA_character_
    out$Fragment_Group_ID_vector <- list(as.character(g$.fragment_group_id))
  } else {
    out$FRAGMENT_GROUP_IDS <- NA_character_
    out$Fragment_Group_ID_vector <- list(character(0))
  }

  if ("FRAGMENTATION_RISK" %in% names(g)) {
    risk_rank <- c(LOW = 1L, MEDIUM = 2L, HIGH = 3L)
    rr <- risk_rank[as.character(g$FRAGMENTATION_RISK)]
    rr <- rr[!is.na(rr)]
    out$FRAGMENTATION_RISK <- if (length(rr)) names(risk_rank)[which(risk_rank == max(rr))[1L]] else "LOW"
  } else {
    out$FRAGMENTATION_RISK <- "LOW"
  }

  if ("FRAGMENTATION_CLASS" %in% names(g)) {
    fc <- unique(as.character(g$FRAGMENTATION_CLASS))
    fc <- fc[!is.na(fc) & nzchar(fc)]
    out$FRAGMENTATION_CLASS_SET <- if (length(fc)) paste(sort(fc), collapse = ";") else "SINGLE_COMPONENT"
  } else {
    out$FRAGMENTATION_CLASS_SET <- "SINGLE_COMPONENT"
  }

  if ("FRAGMENT_Q" %in% names(g)) {
    fq <- safe_num(g$FRAGMENT_Q)
    fq <- fq[is.finite(fq)]
    out$FRAGMENT_Q <- if (length(fq)) min(fq) else NA_real_
  } else {
    out$FRAGMENT_Q <- NA_real_
  }

  out$SIMPLE_REPEAT_OVERLAP <- (
    ("SIMPLE_REPEAT_OVERLAP" %in% names(g) && any(g$SIMPLE_REPEAT_OVERLAP %in% TRUE, na.rm = TRUE)) ||
      ("FRAGMENT_ANY_SIMPLE_REPEAT" %in% names(g) && any(g$FRAGMENT_ANY_SIMPLE_REPEAT %in% TRUE, na.rm = TRUE))
  )

  out$DIST_TO_SIMPLE_REPEAT_BP <- {
    dd <- numeric(0)
    if ("DIST_TO_SIMPLE_REPEAT_BP" %in% names(g)) {
      dd <- c(dd, safe_num(g$DIST_TO_SIMPLE_REPEAT_BP))
    }
    if ("FRAGMENT_MIN_DIST_SIMPLE_REPEAT_BP" %in% names(g)) {
      dd <- c(dd, safe_num(g$FRAGMENT_MIN_DIST_SIMPLE_REPEAT_BP))
    }
    dd <- dd[is.finite(dd)]
    if (length(dd)) min(dd) else NA_real_
  }

  out$FRAGMENT_GROUP_N_COMPONENTS_MAX <- if ("FRAGMENT_GROUP_N_COMPONENTS" %in% names(g)) {
    zz <- safe_num(g$FRAGMENT_GROUP_N_COMPONENTS)
    zz <- zz[is.finite(zz)]
    if (length(zz)) max(zz) else 1
  } else {
    1
  }

  out$FRAGMENT_COMPOUND_SYNTHETIC <- if ("FRAGMENT_COMPOUND_SYNTHETIC" %in% names(g)) {
    any(g$FRAGMENT_COMPOUND_SYNTHETIC %in% TRUE, na.rm = TRUE)
  } else {
    FALSE
  }

  out$FRAGMENT_RECONSTRUCTION_MODE <- if ("FRAGMENT_RECONSTRUCTION_MODE" %in% names(g)) {
    vv <- unique(as.character(g$FRAGMENT_RECONSTRUCTION_MODE))
    vv <- vv[!is.na(vv) & nzchar(vv)]
    if (length(vv)) paste(sort(vv), collapse = ";") else NA_character_
  } else NA_character_

  out$FRAGMENT_RECONSTRUCTION_CAUTION <- if ("FRAGMENT_RECONSTRUCTION_CAUTION" %in% names(g)) {
    any(g$FRAGMENT_RECONSTRUCTION_CAUTION %in% TRUE, na.rm = TRUE)
  } else FALSE

  for (nm in c(
    "FRAGMENT_RECONSTRUCTION_ANCHOR_CLASS",
    "FRAGMENT_RECONSTRUCTION_VCF_SIM",
    "FRAGMENT_RECONSTRUCTION_VCF_CALLERS",
    "FRAGMENT_RECONSTRUCTION_MIN_LR",
    "FRAGMENT_RECONSTRUCTION_ANCHOR_FRACTION",
    "FRAGMENT_COMPONENT_QUERIES",
    "FRAGMENT_COMPONENT_POS_VECTOR",
    "FRAGMENT_COMPONENT_SVLEN_VECTOR"
  )) {
    if (nm %in% names(g)) {
      vv <- unique(as.character(g[[nm]]))
      vv <- vv[!is.na(vv) & nzchar(vv)]
      out[[nm]] <- if (length(vv)) paste(vv, collapse = "|") else NA_character_
    }
  }

  helper_cols <- intersect(
    c(
      ".merge_rowid",
      ".haplotype",
      ".GT_norm",
      ".is_het",
      ".is_hom_alt",
      ".structural_cluster",
      ".source_event_id",
      ".fragment_haplotype",
      ".fragment_group_id"
    ),
    names(out)
  )

  out %>% select(-any_of(helper_cols))
}


# ============================================================
# Tabla final agregada
# ============================================================

cluster_split <- split(
  prefusion_work,
  prefusion_work$.structural_cluster
)

aggregated_list <- lapply(
  cluster_split,
  aggregate_locus_group,
  pair_audit_table = pair_audit
)

FinalRM_MERGE <- bind_rows(aggregated_list) %>%
  arrange(
    factor(CHROM, levels = unique(CHROM)),
    safe_num(POS),
    safe_num(svlen)
  )

FinalRM_MERGE$IDs <- as.character(seq_len(nrow(FinalRM_MERGE)))

if ("query" %in% names(FinalRM_MERGE)) {
  FinalRM_MERGE$query <- NULL
}

# ------------------------------------------------------------------------------
# Aggregate-level fragmentation interpretation.
# A same-source/same-haplotype group that still emitted multiple aggregates is
# a split source.  Small residual aggregates can be labelled SOURCE_SHADOW;
# SOURCE_COMPLEX is never auto-resolved here.
# ------------------------------------------------------------------------------
FinalRM_MERGE$SOURCE_GROUP_SPLIT_N_MAX <- 1L
FinalRM_MERGE$SOURCE_GROUP_SPLIT <- FALSE
FinalRM_MERGE$SOURCE_SHADOW <- FALSE
FinalRM_MERGE$SOURCE_SHADOW_OF <- NA_character_

if ("FRAGMENT_GROUP_IDS" %in% names(FinalRM_MERGE) && nrow(FinalRM_MERGE)) {
  frag_map <- FinalRM_MERGE %>%
    select(IDs, FRAGMENT_GROUP_IDS, POS, svlen) %>%
    tidyr::separate_rows(FRAGMENT_GROUP_IDS, sep = ";") %>%
    filter(!is.na(FRAGMENT_GROUP_IDS), nzchar(FRAGMENT_GROUP_IDS))

  if (nrow(frag_map)) {
    frag_counts <- frag_map %>%
      count(FRAGMENT_GROUP_IDS, name = "N_AGG_OUTPUTS")

    split_stats <- frag_map %>%
      left_join(frag_counts, by = "FRAGMENT_GROUP_IDS") %>%
      group_by(IDs) %>%
      summarise(
        SOURCE_GROUP_SPLIT_N_MAX = max(N_AGG_OUTPUTS, na.rm = TRUE),
        SOURCE_GROUP_SPLIT = any(N_AGG_OUTPUTS > 1L),
        .groups = "drop"
      )

    FinalRM_MERGE <- FinalRM_MERGE %>%
      select(-SOURCE_GROUP_SPLIT_N_MAX, -SOURCE_GROUP_SPLIT) %>%
      left_join(split_stats, by = "IDs") %>%
      mutate(
        SOURCE_GROUP_SPLIT_N_MAX = replace_na(SOURCE_GROUP_SPLIT_N_MAX, 1L),
        SOURCE_GROUP_SPLIT = replace_na(SOURCE_GROUP_SPLIT, FALSE)
      )

    shadow_map <- list()
    split_groups <- frag_counts %>% filter(N_AGG_OUTPUTS > 1L) %>% pull(FRAGMENT_GROUP_IDS)

    for (gid in split_groups) {
      gg <- frag_map %>%
        filter(FRAGMENT_GROUP_IDS == gid) %>%
        distinct(IDs, POS, svlen) %>%
        mutate(
          POS = frag_safe_num(POS),
          svlen = abs(frag_safe_num(svlen))
        ) %>%
        filter(is.finite(POS), is.finite(svlen), svlen > 0)

      if (nrow(gg) < 2L) next
      main_i <- order(-gg$svlen, gg$POS, gg$IDs)[1L]
      main <- gg[main_i, , drop = FALSE]
      others <- gg[-main_i, , drop = FALSE]
      if (!nrow(others)) next

      is_shadow <- (
        others$svlen / main$svlen <= FRAG$shadow_max_fraction &
          abs(others$POS - main$POS) <= FRAG$audit_window_bp
      )
      if (any(is_shadow)) {
        shadow_map[[length(shadow_map) + 1L]] <- data.frame(
          IDs = others$IDs[is_shadow],
          SOURCE_SHADOW = TRUE,
          SOURCE_SHADOW_OF = main$IDs[[1L]],
          stringsAsFactors = FALSE
        )
      }
    }

    if (length(shadow_map)) {
      sh <- bind_rows(shadow_map) %>% distinct(IDs, .keep_all = TRUE)
      FinalRM_MERGE <- FinalRM_MERGE %>%
        select(-SOURCE_SHADOW, -SOURCE_SHADOW_OF) %>%
        left_join(sh, by = "IDs") %>%
        mutate(
          SOURCE_SHADOW = replace_na(SOURCE_SHADOW, FALSE),
          SOURCE_SHADOW_OF = as.character(SOURCE_SHADOW_OF)
        )
    }
  }
}

FinalRM_MERGE <- FinalRM_MERGE %>%
  mutate(
    SOURCE_COMPLEX = grepl("SOURCE_COMPLEX", dplyr::coalesce(as.character(FRAGMENTATION_CLASS_SET), ""), fixed = TRUE),
    FRAGMENT_UNRESOLVED_SUSPECT = (
      grepl("COMPOUND_WITH_CLASSIFIER_DROPOUT", dplyr::coalesce(as.character(FRAGMENTATION_CLASS_SET), ""), fixed = TRUE) |
        grepl("COMPOUND_REVIEW", dplyr::coalesce(as.character(FRAGMENTATION_CLASS_SET), ""), fixed = TRUE)
    ) & !(FRAGMENT_COMPOUND_SYNTHETIC %in% TRUE),
    POPULATION_FRAGMENT_STATUS = case_when(
      SOURCE_SHADOW ~ "SOURCE_SHADOW",
      FRAGMENT_RECONSTRUCTION_CAUTION %in% TRUE ~ "COMPOUND_RECONSTRUCTED_ANCHORED",
      FRAGMENT_COMPOUND_SYNTHETIC %in% TRUE ~ "COMPOUND_RECONSTRUCTED_HIGHCONF",
      SOURCE_COMPLEX ~ "SOURCE_COMPLEX_CAUTION",
      FRAGMENT_UNRESOLVED_SUSPECT ~ "FRAGMENT_SUSPECT_CAUTION",
      FRAGMENTATION_RISK == "HIGH" ~ "FRAGMENTATION_HIGH_RISK",
      FRAGMENTATION_RISK == "MEDIUM" ~ "FRAGMENTATION_MEDIUM_RISK",
      TRUE ~ "CLEAN"
    ),
    POPULATION_PRIMARY_ELIGIBLE = case_when(
      # v5 evidence-aware policy: only a demonstrated shadow is prevented from
      # founding a locus. Complex/dropout groups remain structurally visible so
      # exact cross-sample matches are not destroyed.
      FRAG$population_policy == "evidence_aware" & SOURCE_SHADOW ~ FALSE,
      FRAG$population_policy == "exclude_shadow_complex" &
        (SOURCE_SHADOW | SOURCE_COMPLEX) ~ FALSE,
      FRAG$population_policy == "exclude_fragmented_suspect" &
        (
          SOURCE_SHADOW | SOURCE_COMPLEX | FRAGMENT_UNRESOLVED_SUSPECT |
            (
              grepl("COMPOUND_HIGHCONF", dplyr::coalesce(as.character(FRAGMENTATION_CLASS_SET), ""), fixed = TRUE) &
                !(FRAGMENT_COMPOUND_SYNTHETIC %in% TRUE)
            )
        ) ~ FALSE,
      TRUE ~ TRUE
    ),
    # A cautious reconstruction can participate in identity/shared loci, but a
    # singleton built from unresolved/complex evidence is not allowed to become
    # a clean biological differential without independent rescue context.
    POPULATION_DIFFERENTIAL_ELIGIBLE = !(
      SOURCE_SHADOW | SOURCE_COMPLEX | FRAGMENT_UNRESOLVED_SUSPECT |
        (FRAGMENT_RECONSTRUCTION_CAUTION %in% TRUE)
    ),
    POPULATION_EVIDENCE_ONLY = !POPULATION_PRIMARY_ELIGIBLE,
    ROW_ROLE = case_when(
      !POPULATION_PRIMARY_ELIGIBLE ~ "EVIDENCE_ONLY_AGGREGATE",
      !POPULATION_DIFFERENTIAL_ELIGIBLE ~ "LOCUS_CANDIDATE_CAUTION",
      TRUE ~ "LOCUS_CANDIDATE"
    )
  )

fragmentation_ambiguous <- FinalRM_MERGE %>%
  filter(SOURCE_SHADOW | SOURCE_COMPLEX | FRAGMENTATION_RISK %in% c("HIGH", "MEDIUM"))

population_input_fragmentaware <- FinalRM_MERGE %>%
  filter(POPULATION_PRIMARY_ELIGIBLE)

# Raw component evidence is retained in the SAME final RDS.  These rows never
# found population loci; population may only use them after loci are frozen to
# test direct/compound presence hypotheses.
fragment_evidence_final <- prefusion_all %>%
  mutate(
    COMPONENT_RESOLVED = frag_resolved(event_class) & frag_resolved(event),
    COMPONENT_LR = frag_safe_num(LR),
    SOURCE_EVENT_ID = .source_event_id,
    FRAGMENT_HAPLOTYPE = .fragment_haplotype,
    FRAGMENT_GROUP_ID = .fragment_group_id,
    ROW_ROLE = "FRAGMENT_EVIDENCE"
  ) %>%
  filter(
    !is.na(FRAGMENT_GROUP_ID), nzchar(FRAGMENT_GROUP_ID),
    (
      replace_na(FRAGMENT_GROUP_N_COMPONENTS, 1L) >= 2L |
        SOURCE_COMPLEX_GLOBAL %in% TRUE |
        !COMPONENT_RESOLVED
    )
  )


# ============================================================
# Auditoria opcional.  Por defecto v5 escribe un único RDS final.
# ============================================================

if (write_audits) {
if (nrow(pair_audit) > 0) {
  data.table::fwrite(
    pair_audit,
    file = file.path(
      out,
      paste0(sample, "_intrasample_pair_audit.tsv")
    ),
    sep = "\t",
    quote = FALSE,
    na = "NA"
  )
}

aggregate_audit <- FinalRM_MERGE %>%
  transmute(
    IDs,
    CHROM,
    POS,
    POS_MEDIAN,
    svlen,
    SVLEN_MEDIAN,
    SVLEN_MIN,
    SVLEN_MAX,
    SVLEN_IQR,
    SVLEN_MINMAX_SIMILARITY,
    event_class,
    event,
    completeness,
    GT,
    GT_STATUS,
    ALT_DOSAGE,
    LOCUS_STATUS,
    SIZE_VARIATION_STATUS,
    ANNOTATION_STATUS,
    RECONSTRUCTION_QUALITY_SCORE,
    RECONSTRUCTION_QUALITY_LABEL,
    RECONSTRUCTION_QUALITY_REASON,
    RQ_POSITION,
    RQ_SIZE,
    RQ_ANNOTATION,
    RQ_PHASE,
    N_REPRESENTATIONS,
    N_SOURCE_IDS,
    SOURCE_IDS,
    N_HAPLOTYPES,
    N_ALT_HAPLOTYPES,
    HAPLOTYPE_SET,
    SVLEN_BY_HAPLOTYPE,
    STRUCTURAL_COPY_COUNT,
    CLASS_SET,
    CLASS_COUNTS,
    EVENT_SET,
    EVENT_COUNTS,
    COMPLETENESS_SET,
    COMPLETENESS_COUNTS,
    MERGE_RULES,
    MERGE_RELATIONS,
    POS_MIN,
    POS_MAX,
    POS_SPAN,
    POS_IQR,
    SOURCE_EVENT_IDS,
    FRAGMENT_GROUP_IDS,
    FRAGMENTATION_RISK,
    FRAGMENTATION_CLASS_SET,
    FRAGMENT_Q,
    SIMPLE_REPEAT_OVERLAP,
    DIST_TO_SIMPLE_REPEAT_BP,
    FRAGMENT_GROUP_N_COMPONENTS_MAX,
    FRAGMENT_COMPOUND_SYNTHETIC,
    SOURCE_GROUP_SPLIT_N_MAX,
    SOURCE_GROUP_SPLIT,
    SOURCE_SHADOW,
    SOURCE_SHADOW_OF,
    SOURCE_COMPLEX,
    POPULATION_FRAGMENT_STATUS,
    POPULATION_PRIMARY_ELIGIBLE
  )

data.table::fwrite(
  aggregate_audit,
  file = file.path(
    out,
    paste0(sample, "_intrasample_aggregate_audit.tsv")
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

# Fragmentation audit outputs.
fragment_component_audit <- prefusion_all %>%
  transmute(
    query,
    SOURCE_EVENT_ID = .source_event_id,
    FRAGMENT_HAPLOTYPE = .fragment_haplotype,
    FRAGMENT_GROUP_ID = .fragment_group_id,
    CHROM,
    POS,
    svlen,
    GT,
    LR,
    event_class,
    event,
    completeness,
    structural_pass = if ("structural_pass" %in% names(prefusion_all)) structural_pass else NA,
    cover_status,
    region_status,
    SIMPLE_REPEAT_OVERLAP,
    SIMPLE_REPEAT_N,
    DIST_TO_SIMPLE_REPEAT_BP,
    FRAGMENTATION_RISK = if ("FRAGMENTATION_RISK" %in% names(prefusion_all)) FRAGMENTATION_RISK else NA_character_,
    FRAGMENTATION_CLASS = if ("FRAGMENTATION_CLASS" %in% names(prefusion_all)) FRAGMENTATION_CLASS else NA_character_,
    FRAGMENT_Q = if ("FRAGMENT_Q" %in% names(prefusion_all)) FRAGMENT_Q else NA_real_,
    FRAGMENT_GROUP_N_COMPONENTS = if ("FRAGMENT_GROUP_N_COMPONENTS" %in% names(prefusion_all)) FRAGMENT_GROUP_N_COMPONENTS else NA_integer_,
    FRAGMENT_COMPOUND_SVLEN = if ("FRAGMENT_COMPOUND_SVLEN" %in% names(prefusion_all)) FRAGMENT_COMPOUND_SVLEN else NA_real_,
    FRAGMENT_BEST_VCF_ALT_SIMILARITY = if ("FRAGMENT_BEST_VCF_ALT_SIMILARITY" %in% names(prefusion_all)) FRAGMENT_BEST_VCF_ALT_SIMILARITY else NA_real_,
    FRAGMENT_N_VCF_ALT_CALLERS_GE_STRONG = if ("FRAGMENT_N_VCF_ALT_CALLERS_GE_STRONG" %in% names(prefusion_all)) FRAGMENT_N_VCF_ALT_CALLERS_GE_STRONG else NA_integer_
  )

data.table::fwrite(
  fragment_component_audit,
  file.path(out, paste0(sample, "_fragmentation_component_audit.tsv")),
  sep = "\t", quote = FALSE, na = "NA"
)

# Full pre-filter object keeps ALT/contig/classifier fields for forensic review.
saveRDS(
  prefusion_all,
  file.path(out, paste0(sample, "_fragmentation_components_full.rds"))
)

if (nrow(fragment_group_audit)) {
  data.table::fwrite(
    fragment_group_audit,
    file.path(out, paste0(sample, "_fragmentation_group_audit.tsv")),
    sep = "\t", quote = FALSE, na = "NA"
  )
  saveRDS(
    fragment_group_audit,
    file.path(out, paste0(sample, "_fragmentation_group_audit.rds"))
  )
}

if (nrow(fragment_source_audit)) {
  data.table::fwrite(
    fragment_source_audit,
    file.path(out, paste0(sample, "_fragmentation_source_audit.tsv")),
    sep = "\t", quote = FALSE, na = "NA"
  )
  saveRDS(
    fragment_source_audit,
    file.path(out, paste0(sample, "_fragmentation_source_audit.rds"))
  )
}

if (nrow(fragmentation_ambiguous)) {
  # Compact scalar audit; list columns remain in the RDS outputs.
  frag_ambig_tsv <- fragmentation_ambiguous %>%
    select(
      IDs, CHROM, POS, svlen, event_class, event, completeness,
      RECONSTRUCTION_QUALITY_SCORE, SOURCE_EVENT_IDS, FRAGMENT_GROUP_IDS,
      FRAGMENTATION_RISK, FRAGMENTATION_CLASS_SET, FRAGMENT_Q,
      SIMPLE_REPEAT_OVERLAP, DIST_TO_SIMPLE_REPEAT_BP,
      SOURCE_GROUP_SPLIT_N_MAX, SOURCE_SHADOW, SOURCE_SHADOW_OF, SOURCE_COMPLEX,
      POPULATION_FRAGMENT_STATUS, POPULATION_PRIMARY_ELIGIBLE
    )
  data.table::fwrite(
    frag_ambig_tsv,
    file.path(out, paste0(sample, "_fragmentation_ambiguous_aggregates.tsv")),
    sep = "\t", quote = FALSE, na = "NA"
  )
}

fragment_parameter_table <- data.frame(
  parameter = c(
    "fragmentation_mode", "fragment_population_policy",
    "simple_repeats_bed", "vcf_cutesv", "vcf_severus", "vcf_sniffles2",
    paste0("CFG.", names(CFG)), paste0("FRAG.", names(FRAG))
  ),
  value = c(
    FRAG$mode, FRAG$population_policy,
    as.character(simple_repeats_bed), as.character(vcf_cutesv),
    as.character(vcf_severus), as.character(vcf_sniffles2),
    vapply(CFG, function(x) paste(x, collapse = ";"), character(1)),
    vapply(FRAG, function(x) paste(x, collapse = ";"), character(1))
  ),
  stringsAsFactors = FALSE
)
data.table::fwrite(
  fragment_parameter_table,
  file.path(out, paste0(sample, "_fragmentation_parameters.tsv")),
  sep = "\t", quote = FALSE, na = "NA"
)

} # write_audits

# Re-evaluate the actual reconstructed breakpoint, not inherited component flags.
primary_sets <- primary_partition(FinalRM_MERGE, coverage, CFG$primary_min_support, CFG$min_svlen)
FinalRM_MERGE <- primary_sets$loci
write_flat_tsv(primary_sets$nocallablecandidates, file.path(out,paste0(sample,"_nocallablecandidates.tsv")))
write_flat_tsv(primary_sets$rescue_candidates, file.path(out,paste0(sample,"_rescue_candidates.tsv")))
# One canonical intra-sample object.  Population v5 reads $loci for locus
# construction and $fragment_evidence only after loci are frozen.
intra_final <- list(
  format_version = "Retroposon_intra_v5_fragment_rescue_2026-10-07",
  sample_id = sample,
  loci = FinalRM_MERGE,
  fragment_evidence = fragment_evidence_final,
  rescue_candidates = primary_sets$rescue_candidates,
  nocallablecandidates = primary_sets$nocallablecandidates,
  parameters = list(CFG = CFG, FRAG = FRAG),
  counts = list(
    n_loci = nrow(FinalRM_MERGE),
    n_population_primary_eligible = sum(FinalRM_MERGE$POPULATION_PRIMARY_ELIGIBLE %in% TRUE, na.rm = TRUE),
    n_population_differential_eligible = sum(FinalRM_MERGE$POPULATION_DIFFERENTIAL_ELIGIBLE %in% TRUE, na.rm = TRUE),
    n_evidence_only_aggregates = sum(FinalRM_MERGE$POPULATION_EVIDENCE_ONLY %in% TRUE, na.rm = TRUE),
    n_fragment_evidence_components = nrow(fragment_evidence_final),
    n_shadow = sum(FinalRM_MERGE$SOURCE_SHADOW %in% TRUE, na.rm = TRUE),
    n_complex_group = sum(FinalRM_MERGE$SOURCE_COMPLEX %in% TRUE, na.rm = TRUE),
    n_compound_synthetic = sum(FinalRM_MERGE$FRAGMENT_COMPOUND_SYNTHETIC %in% TRUE, na.rm = TRUE),
    n_compound_anchored = sum(FinalRM_MERGE$FRAGMENT_RECONSTRUCTION_CAUTION %in% TRUE, na.rm = TRUE)
  )
)

saveRDS(
  intra_final,
  file.path(out, paste0(sample, "_retroatlas_intra_final.rds"))
)

if (write_legacy_rds) {
  saveRDS(FinalRM_MERGE, file.path(out, paste0(sample, "_aggregate_full.rds")))
  saveRDS(FinalRM_MERGE, file.path(out, paste0(sample, "_population_input.rds")))
  saveRDS(population_input_fragmentaware, file.path(out, paste0(sample, "_population_input_fragmentaware.rds")))
  saveRDS(fragmentation_ambiguous, file.path(out, paste0(sample, "_fragmentation_ambiguous_aggregates.rds")))
}

message(
  "Fragment-aware population rows: ", nrow(population_input_fragmentaware),
  " / ", nrow(FinalRM_MERGE),
  " | shadows=", sum(FinalRM_MERGE$SOURCE_SHADOW, na.rm = TRUE),
  " | complex=", sum(FinalRM_MERGE$SOURCE_COMPLEX, na.rm = TRUE),
  " | compound_synthetic=", sum(FinalRM_MERGE$FRAGMENT_COMPOUND_SYNTHETIC, na.rm = TRUE)
)


# ============================================================
# SUBSET PROCESABLE PARA TSD / CONTIG / BLAST
# ============================================================

df_limpio <- FinalRM_MERGE %>%
  filter(
    # El representante del agregado debe conservar un contig real
    # para entrar en TSD / BLAST. La tabla completa ya se ha guardado
    # como *_aggregate_full.rds.
    !is.na(contig),
    nzchar(contig)
  )


# ============================================================
# TSD
# ============================================================

w <- arg_int("--tsd-flank-bp", 2000)
chrom_limits <- mapping$expected_length[match(df_limpio$CHROM, mapping$CHROM)]

chordenates <- df_limpio %>%
  mutate(
    position = glue("{CHROM}:{pmax(1, POS - w)}-{pmin(chrom_limits, POS + w)}")
  ) %>%
  select(
    IDs,
    position,
    svlen
  )

write.table(
  chordenates,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE,
  file = glue("{out}/{sample}_samtools.bed")
)


library(Biostrings)

objeto_contig <- DNAStringSet(df_limpio$contig)
names(objeto_contig) <- df_limpio$IDs

writeXStringSet(
  objeto_contig,
  filepath = file.path(
    out,
    paste0(sample, "_contig.fasta")
  ),
  format = "fasta"
)


if (write_legacy_rds) {
  saveRDS(
    df_limpio,
    glue("{out}/{sample}.rds")
  )
}
