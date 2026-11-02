library(GenomicRanges)
library(IRanges)
library(tidyverse)
library(glue)


# ============================================================
# ARGUMENTOS
# ============================================================

log_connection <- file(snakemake@log[[1]], open = "wt")
sink(log_connection)
sink(log_connection, type = "message")
cores <- as.integer(snakemake@params[["chunks"]])
file1 <- snakemake@input[["sniffles"]]
file2 <- snakemake@input[["cutesv"]]
file3 <- snakemake@input[["severus"]]
windows_size <- as.integer(snakemake@params[["window"]])

if (is.na(cores) || cores < 1) {
  stop("cores debe ser un entero >= 1")
}

if (is.na(windows_size) || windows_size < 0) {
  stop("windows_size debe ser un entero >= 0")
}

window <- windows_size


# ============================================================
# FUNCIONES AUXILIARES
# ============================================================


# ------------------------------------------------------------
# Extrae un campo de INFO.
#
# Ejemplo:
#
# INFO = "SVTYPE=DEL;SVLEN=-300;END=1200"
#
# extraer_campo(INFO, "SVLEN")
# -> "-300"
#
# Si el campo no existe devuelve NA.
# ------------------------------------------------------------

extraer_campo <- function(columna_info, etiqueta) {

  pattern <- paste0(
    "(?:^|;)",
    etiqueta,
    "=([^;]+)"
  )

  resultado <- stringr::str_match(
    as.character(columna_info),
    pattern
  )[, 2]

  return(resultado)
}


# ------------------------------------------------------------
# Lectura sencilla de VCF.
#
# Este script, igual que el original, asume:
#
# CHROM POS ID REF ALT QUAL FILTER INFO FORMAT GENOTYPE
#
# es decir, un VCF con una única muestra.
# ------------------------------------------------------------

read_vcf_simple <- function(file) {

  vcf_names <- c(
    "CHROM",
    "POS",
    "ID",
    "REF",
    "ALT",
    "QUAL",
    "FILTER",
    "INFO",
    "FORMAT",
    "GENOTYPE"
  )

  # read.table falla con un VCF sin registros. Mantener las diez columnas.
  connection <- if (grepl("\\.gz$", file)) gzfile(file, "rt") else base::file(file, "rt")
  records <- readLines(connection, warn = FALSE)
  close(connection)
  records <- records[nzchar(records) & !startsWith(records, "#")]
  if (!length(records)) {
    return(as_tibble(setNames(rep(list(character()), length(vcf_names)), vcf_names)))
  }
  x <- read.table(
    text = paste(records, collapse = "\n"),
    sep = "\t",
    comment.char = "#",
    quote = "",
    header = FALSE,
    stringsAsFactors = FALSE,
    fill = FALSE,
    colClasses = "character"
  )

  if (ncol(x) != length(vcf_names)) {

    stop(
      glue(
        "El VCF {file} tiene {ncol(x)} columnas, ",
        "pero este script espera 10 columnas."
      )
    )
  }

  colnames(x) <- vcf_names

  return(
    as_tibble(x)
  )
}


# ============================================================
# LECTURA DE LOS VCF
# ============================================================

message("Leyendo VCF1: ", file1)
vcf1 <- read_vcf_simple(file1)

message("Leyendo VCF2: ", file2)
vcf2 <- read_vcf_simple(file2)

message("Leyendo VCF3: ", file3)
vcf3 <- read_vcf_simple(file3)


# ============================================================
# FILTRADO PASS
# ============================================================

message("Filtrando variantes PASS...")

vcf1 <- vcf1 %>%
  filter(FILTER == "PASS")

vcf2 <- vcf2 %>%
  filter(FILTER == "PASS")

vcf3 <- vcf3 %>%
  filter(FILTER == "PASS")


message(
  "Variantes PASS: ",
  "VCF1=", nrow(vcf1),
  " | VCF2=", nrow(vcf2),
  " | VCF3=", nrow(vcf3)
)


# ============================================================
# COLLAPSE DE LOS 3 VCF
# ============================================================
#
# Filosofía:
#
# Cada variante genera:
#
# POS-window ---------------- POS ---------------- POS+window
#
# Todas las ventanas que se solapan se fusionan mediante
# GenomicRanges::reduce().
#
# IMPORTANTE:
#
# Esto permite deliberadamente cadenas:
#
# A solapa B
# B solapa C
#
# => A/B/C quedan dentro de la misma región Atlas.
#
# Posteriormente podremos separar eventos diferentes si es
# necesario.
#
# Además de generar la región Atlas, guardamos TODAS las
# variantes originales de cada caller que contribuyeron a ella.
# ============================================================


collapse_3vcf <- function(
    vcf1,
    vcf2,
    vcf3,
    window = 500,
    only_unique = FALSE
) {

  # ==========================================================
  # PREPARACIÓN DE CADA CALLER
  # ==========================================================

  preparar_eventos <- function(vcf, caller) {

    vcf %>%

      mutate(

        # ID interno dentro de cada caller
        caller_row = row_number(),

        CHROM = as.character(CHROM),
        POS   = suppressWarnings(as.integer(POS)),
        ID    = as.character(ID),
        REF   = as.character(REF),
        ALT   = as.character(ALT),
        QUAL  = as.character(QUAL),
        FILTER = as.character(FILTER),
        INFO   = as.character(INFO),
        FORMAT = as.character(FORMAT),
        GENOTYPE = as.character(GENOTYPE),

        # ----------------------------------------------------
        # Campos SV obtenidos de INFO
        # ----------------------------------------------------

        SVTYPE = extraer_campo(
          INFO,
          "SVTYPE"
        ),

        SVLEN_raw = extraer_campo(
          INFO,
          "SVLEN"
        ),

        END_raw = extraer_campo(
          INFO,
          "END"
        ),

        # SVLEN puede ocasionalmente tener varios valores.
        # Para approx_size usamos el primero.
        SVLEN = suppressWarnings(
          as.numeric(
            sub(
              ",.*$",
              "",
              SVLEN_raw
            )
          )
        ),

        END = suppressWarnings(
          as.integer(
            sub(
              ",.*$",
              "",
              END_raw
            )
          )
        ),

        # Si SVLEN=0 no aporta información útil de tamaño
        SVLEN = if_else(
          !is.na(SVLEN) & SVLEN != 0,
          SVLEN,
          NA_real_
        ),

        # ----------------------------------------------------
        # Tipo de representación de ALT
        # ----------------------------------------------------

        symbolic_alt = grepl(
          "^<.*>$",
          ALT
        ),

        breakend_alt = grepl(
          "[\\[\\]]",
          ALT
        ),

        # ----------------------------------------------------
        # Tamaño aproximado usando END
        # ----------------------------------------------------
        #
        # Por ejemplo:
        #
        # POS = 1000
        # END = 1300
        #
        # tamaño ~300
        #
        # Utilizamos una estimación porque diferentes callers
        # pueden tener convenciones ligeramente distintas.
        # ----------------------------------------------------

        size_from_end = case_when(

          !is.na(END) &
            !is.na(POS) &
            END != POS ~
            as.numeric(
              abs(END - POS)
            ),

          TRUE ~ NA_real_
        ),

        # ----------------------------------------------------
        # Tamaño derivado de REF/ALT
        #
        # Solamente para variantes sequence-resolved.
        #
        # No intentamos obtenerlo así para:
        #
        # <DEL>
        # <INS>
        # <DUP>
        # BND
        # etc.
        # ----------------------------------------------------

        size_from_sequence = case_when(

          !symbolic_alt &
            !breakend_alt &
            !is.na(REF) &
            !is.na(ALT) ~
            as.numeric(
              abs(
                nchar(ALT) -
                  nchar(REF)
              )
            ),

          TRUE ~ NA_real_
        ),

        size_from_sequence = if_else(
          !is.na(size_from_sequence) &
            size_from_sequence > 0,
          size_from_sequence,
          NA_real_
        ),

        # ----------------------------------------------------
        # Tamaño aproximado definitivo
        #
        # Prioridad:
        #
        # 1. SVLEN
        # 2. END - POS
        # 3. diferencia REF/ALT
        #
        # Si ninguno funciona -> NA
        #
        # IMPORTANTE:
        #
        # La variante NO se elimina si approx_size es NA.
        # Simplemente queda registrado que había una llamada,
        # aunque no conozcamos su tamaño.
        # ----------------------------------------------------

        approx_size = dplyr::coalesce(
          abs(SVLEN),
          size_from_end,
          size_from_sequence
        ),

        size_method = case_when(

          !is.na(SVLEN) ~
            "SVLEN",

          !is.na(size_from_end) ~
            "END-POS",

          !is.na(size_from_sequence) ~
            "REF-ALT",

          TRUE ~
            NA_character_
        )
      ) %>%

      # ------------------------------------------------------
      # Filtrado mínimo
      #
      # NO eliminamos ALT simbólicos como <DEL>/<INS>.
      # ------------------------------------------------------

      filter(

        !is.na(CHROM),

        !is.na(POS),

        !is.na(ALT),

        ALT != "",

        ALT != ".",

        ALT != "*"
      ) %>%

      mutate(

        caller = caller,

        event_id = paste0(
          caller,
          "_",
          caller_row
        ),

        # Ventana de detección
        start_window = pmax(
          1L,
          POS - as.integer(window)
        ),

        end_window =
          POS + as.integer(window)
      ) %>%

      select(

        event_id,

        caller,

        caller_row,

        CHROM,

        POS,

        ID,

        REF,

        ALT,

        QUAL,

        FILTER,

        INFO,

        FORMAT,

        GENOTYPE,

        SVTYPE,

        SVLEN,

        END,

        approx_size,

        size_method,

        symbolic_alt,

        breakend_alt,

        start_window,

        end_window
      )
  }


  # ==========================================================
  # PREPARAMOS LOS TRES VCF
  # ==========================================================

  eventos_vcf1 <- preparar_eventos(
    vcf1,
    "vcf1"
  )

  eventos_vcf2 <- preparar_eventos(
    vcf2,
    "vcf2"
  )

  eventos_vcf3 <- preparar_eventos(
    vcf3,
    "vcf3"
  )


  eventos <- bind_rows(

    eventos_vcf1,

    eventos_vcf2,

    eventos_vcf3

  ) %>%

    mutate(
      event_index = row_number()
    )


  # ==========================================================
  # CASO: NO HAY NINGÚN EVENTO
  # ==========================================================

  if (nrow(eventos) == 0) {

    empty_detail <- eventos %>%
      select(
        -event_index,
        -start_window,
        -end_window
      )

    resultado <- tibble(

      Id = character(),

      chr = character(),

      start_win = integer(),

      end_win = integer(),

      bed_start = integer(),

      bed_end = integer(),

      callers = character(),

      n_callers = integer(),

      has_vcf1 = logical(),

      has_vcf2 = logical(),

      has_vcf3 = logical(),

      n_vcf1 = integer(),

      n_vcf2 = integer(),

      n_vcf3 = integer(),

      sizes_vcf1 = list(),

      sizes_vcf2 = list(),

      sizes_vcf3 = list(),

      vcf1 = list(),

      vcf2 = list(),

      vcf3 = list()
    )

    return(resultado)
  }


  # ==========================================================
  # GRANGES CON LAS VENTANAS
  # ==========================================================

  gr <- GRanges(

    seqnames = eventos$CHROM,

    ranges = IRanges(

      start = eventos$start_window,

      end = eventos$end_window
    ),

    event_index =
      eventos$event_index
  )


  # ==========================================================
  # COLLAPSE DE TODAS LAS VENTANAS SOLAPANTES
  # ==========================================================

  regiones <- GenomicRanges::reduce(gr)


  # ==========================================================
  # RELACIONAR CADA EVENTO ORIGINAL CON SU REGIÓN REDUCIDA
  # ==========================================================

  hits <- findOverlaps(
    regiones,
    gr
  )


  mapa_eventos <- tibble(

    region = queryHits(hits),

    event_index =
      gr$event_index[
        subjectHits(hits)
      ]

  ) %>%

    left_join(
      eventos,
      by = "event_index"
    )


  # ==========================================================
  # TABLA CON LAS COORDENADAS DE CADA REGIÓN
  # ==========================================================

  regiones_df <- tibble(

    region =
      seq_along(regiones),

    chr =
      as.character(
        seqnames(regiones)
      ),

    start_win =
      start(regiones),

    end_win =
      end(regiones),

    # Coordenadas BED correctas:
    #
    # BED start = 0-based
    # BED end   = end-exclusive
    #
    # Para un intervalo 1-based inclusivo:
    #
    # start_win - 1
    # end_win
    bed_start =
      pmax(
        0L,
        start(regiones) - 1L
      ),

    bed_end =
      end(regiones)
  )


  # ==========================================================
  # COLUMNAS QUE QUEREMOS CONSERVAR DE CADA VARIANTE
  # ==========================================================

  detail_columns <- c(

    "event_id",

    "caller",

    "caller_row",

    "CHROM",

    "POS",

    "ID",

    "REF",

    "ALT",

    "QUAL",

    "FILTER",

    "INFO",

    "FORMAT",

    "GENOTYPE",

    "SVTYPE",

    "SVLEN",

    "END",

    "approx_size",

    "size_method",

    "symbolic_alt",

    "breakend_alt"
  )


  obtener_detalle <- function(df, caller_name) {

    df %>%

      filter(
        caller == caller_name
      ) %>%

      select(
        all_of(
          detail_columns
        )
      )
  }


  # ==========================================================
  # RESUMEN POR REGIÓN
  # ==========================================================
  #
  # split() evita tener que filtrar todo mapa_eventos
  # repetidamente para cada región.
  # ==========================================================

  eventos_por_region <- split(
    mapa_eventos,
    mapa_eventos$region
  )


  resumen_eventos <- purrr::imap_dfr(

    eventos_por_region,

    function(df, region_name) {

      sizes1 <- df$approx_size[
        df$caller == "vcf1"
      ]

      sizes2 <- df$approx_size[
        df$caller == "vcf2"
      ]

      sizes3 <- df$approx_size[
        df$caller == "vcf3"
      ]


      detalle1 <- obtener_detalle(
        df,
        "vcf1"
      )

      detalle2 <- obtener_detalle(
        df,
        "vcf2"
      )

      detalle3 <- obtener_detalle(
        df,
        "vcf3"
      )


      tibble(

        region =
          as.integer(region_name),

        callers =
          paste(
            sort(
              unique(
                df$caller
              )
            ),
            collapse = ","
          ),

        n_callers =
          n_distinct(
            df$caller
          ),

        # ----------------------------------------------------
        # Presencia / ausencia
        # ----------------------------------------------------

        has_vcf1 =
          nrow(detalle1) > 0,

        has_vcf2 =
          nrow(detalle2) > 0,

        has_vcf3 =
          nrow(detalle3) > 0,

        # ----------------------------------------------------
        # Número de llamadas por caller
        # ----------------------------------------------------

        n_vcf1 =
          nrow(detalle1),

        n_vcf2 =
          nrow(detalle2),

        n_vcf3 =
          nrow(detalle3),

        # ----------------------------------------------------
        # Vectores de tamaños
        #
        # Ejemplos:
        #
        # numeric(0)
        # -> ninguna llamada
        #
        # c(300)
        # -> una llamada ~300 bp
        #
        # c(300, 700)
        # -> dos llamadas
        #
        # c(300, NA, 700)
        # -> tres llamadas; una sin tamaño estimable
        # ----------------------------------------------------

        sizes_vcf1 =
          list(sizes1),

        sizes_vcf2 =
          list(sizes2),

        sizes_vcf3 =
          list(sizes3),

        # ----------------------------------------------------
        # Registro completo de las variantes
        #
        # Cada celda contiene un tibble.
        # ----------------------------------------------------

        vcf1 =
          list(detalle1),

        vcf2 =
          list(detalle2),

        vcf3 =
          list(detalle3)
      )
    }
  )


  # ==========================================================
  # UNIR COORDENADAS + INFORMACIÓN DE LOS CALLERS
  # ==========================================================

  resultado <- regiones_df %>%

    left_join(
      resumen_eventos,
      by = "region"
    )


  # ==========================================================
  # OPCIONAL: SOLO REGIONES DE UN CALLER
  # ==========================================================

  if (only_unique) {

    resultado <- resultado %>%

      filter(
        n_callers == 1
      )
  }


  # ==========================================================
  # ORDEN Y IDs
  # ==========================================================

  resultado <- resultado %>%

    arrange(
      chr,
      start_win,
      end_win
    ) %>%

    mutate(

      Id = paste0(
        "Atlas_",
        seq_len(n())
      )
    ) %>%

    select(

      Id,

      chr,

      start_win,

      end_win,

      bed_start,

      bed_end,

      callers,

      n_callers,

      has_vcf1,

      has_vcf2,

      has_vcf3,

      n_vcf1,

      n_vcf2,

      n_vcf3,

      sizes_vcf1,

      sizes_vcf2,

      sizes_vcf3,

      vcf1,

      vcf2,

      vcf3
    )


  return(resultado)
}


# ============================================================
# GENERAR REGIONES ATLAS
# ============================================================

message("Generando regiones Atlas...")

vcf_info <- collapse_3vcf(

  vcf1 = vcf1,

  vcf2 = vcf2,

  vcf3 = vcf3,

  window = window,

  only_unique = FALSE
)


message(
  "Número total de regiones Atlas: ",
  nrow(vcf_info)
)


# ============================================================
# AÑADIR METADATOS AL OBJETO RDS
# ============================================================
#
# saveRDS conserva estos atributos.
#
# Se pueden consultar posteriormente con:
#
# attr(x, "window_size")
# attr(x, "source_files")
# attr(x, "size_definition")
# ============================================================

attr(
  vcf_info,
  "window_size"
) <- window


attr(
  vcf_info,
  "source_files"
) <- c(

  vcf1 = file1,

  vcf2 = file2,

  vcf3 = file3
)


attr(
  vcf_info,
  "size_definition"
) <- paste0(
  "approx_size priority: ",
  "abs(SVLEN) -> abs(END-POS) -> abs(length(ALT)-length(REF)). ",
  "NA means that a variant was present but its size could not be inferred."
)


# ============================================================
# GUARDAR RDS DE RESPALDO
# ============================================================

rds_file <- snakemake@output[["rds"]]


saveRDS(
  vcf_info,
  file = rds_file
)


message(
  "RDS guardado: ",
  rds_file
)


# ============================================================
# GENERAR BED PARA MOSDEPTH
# ============================================================
#
# BED:
#
# chr
# start  -> 0-based
# end    -> end-exclusive
# Id
# ============================================================

mosdepth_bed <- vcf_info %>%

  transmute(

    chr = chr,

    start = bed_start,

    end = bed_end,

    locus_id = Id
  ) %>%

  arrange(
    chr,
    start,
    end
  )


mosdepth_file <- snakemake@output[["bed"]]


write.table(

  mosdepth_bed,

  file = mosdepth_file,

  sep = "\t",

  quote = FALSE,

  row.names = FALSE,

  col.names = FALSE
)


message(
  "BED para mosdepth guardado: ",
  mosdepth_file
)


# ============================================================
# GENERAR CHUNKS
# ============================================================

path2 <- snakemake@output[["chunks"]]
# Eliminar únicamente chunks anteriores de esta regla en una reejecución.
unlink(list.files(path2, pattern = "^windows_[0-9]+\\.bed$", full.names = TRUE))


dir.create(

  path2,

  showWarnings = FALSE,

  recursive = TRUE
)


# ------------------------------------------------------------
# Para los chunks usamos también coordenadas BED correctas.
#
# En el código original los chunks utilizaban start_win,
# que es 1-based, aunque el fichero terminaba en .bed.
#
# Aquí ambos BED usan la misma convención:
#
# chr  0-based-start  end  Atlas_ID
# ------------------------------------------------------------

Final_chunks <- vcf_info %>%

  transmute(

    chr = chr,

    start = bed_start,

    end = bed_end,

    Id = Id
  )


if (nrow(Final_chunks) > 0) {

  # Nunca crear más chunks que regiones disponibles
  n_chunks <- min(
    cores,
    nrow(Final_chunks)
  )


  if (n_chunks == 1) {

    chunk_id <- rep(
      1L,
      nrow(Final_chunks)
    )

  } else {

    chunk_id <- cut(

      seq_len(
        nrow(Final_chunks)
      ),

      breaks = n_chunks,

      labels = FALSE
    )
  }


  chunks <- split(
    Final_chunks,
    chunk_id
  )


  for (i in seq_along(chunks)) {

    output_file <- sprintf(

      "%s/windows_%02d.bed",

      path2,

      i
    )


    write.table(

      chunks[[i]],

      file = output_file,

      sep = "\t",

      quote = FALSE,

      row.names = FALSE,

      col.names = FALSE
    )
  }


  message(
    "Chunks generados: ",
    length(chunks),
    " en ",
    path2
  )

} else {

  message(
    "No existen regiones Atlas; no se generan chunks."
  )
}


# ============================================================
# RESUMEN
# ============================================================

message("")
message("============================================")
message("PIPELINE DE GENERACIÓN DE VENTANAS TERMINADO")
message("============================================")
message("")
message("Regiones Atlas: ", nrow(vcf_info))
message("RDS:             ", rds_file)
message("BED mosdepth:    ", mosdepth_file)
message("Chunks:          ", path2)
message("")
sink(type = "message")
sink()
close(log_connection)
