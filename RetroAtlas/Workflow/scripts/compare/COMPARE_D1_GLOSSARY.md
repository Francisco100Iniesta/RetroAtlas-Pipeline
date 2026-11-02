# Pipeline Compare: rescate antes de callability

Esta adaptación usa la profundidad regional de `regions.bed` (mediana por ventana de Mosdepth con `-m`). Por defecto exige cobertura positiva y hasta 3 veces la mediana aplicable, sin mínimo de 20×. X/Y tienen medianas independientes en hombres; mujeres usan la autosómica para X y excluyen Y. Adaptive Sampling limita todas las estadísticas al panel. No se usan archivos thresholds.

Se evalúan VCF, secuencias softclip y fragmentos incluso en regiones no callable. RESCUE conserva la evidencia positiva; REVIEW conserva ambigüedades y duplicados. NOT_CALLABLE se asigna después de descartar rescate cuando no se puede evaluar ausencia; DIFFERENTIAL requiere evidencia negativa evaluable y región callable. Un VCF vacío válido cuenta como búsqueda completada.

Los detalles del algoritmo Q y las defensas frente a duplicados proceden del módulo aportado. Los parámetros efectivos del RDS y advance_compare.yaml prevalecen sobre los valores históricos descritos a continuación.

# RetroAtlas Compare D1 — Categories and Parameter Glossary

The following descriptions are written in simple scientific English, suitable for a methods document, results table, README, or supplementary materials. They correspond to the parameters used in your latest D1 comparison.

## 1. Result categories

| Category      | Short definition                                                                                                                                                                                      |
| ------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| PRIMARY_MATCH | An insertion independently reconstructed in both samples and matched based on genomic position, insertion size, structural annotations, and sequence evidence when evaluated.                         |
| RESCUE        | An insertion found as PRIMARY in one sample and supported by additional evidence in the other, such as VCF calls, soft-clipped reads, or structural fragments.                                        |
| REVIEW        | An insertion with ambiguous or conflicting evidence that cannot be confidently classified as shared or absent. Includes nearby competing insertions and possible alternative allelic representations. |
| NOT_CALLABLE  | An insertion whose presence or absence cannot be assessed reliably in the other sample because the genomic region does not meet the required coverage or quality criteria.                            |
| DIFFERENTIAL  | An insertion detected as PRIMARY in one sample but not supported in the other, despite adequate regional callability and the required negative evidence checks.                                       |

Important distinction: `REVIEW` is not a negative result. `NOT_CALLABLE` is not a confirmed absence. `RESCUE` does not always provide a reliable genotype.

### Complete PRIMARY matching quality scale (Q10–Q1)

The Q tier indicates the matching rule that was satisfied by independently reconstructed PRIMARY insertions. Tiers consider genomic breakpoint distance (ΔPOS), insertion-length similarity, and agreement of repeat class, event, and completeness. A higher Q generally reflects stricter structural agreement **within the same-class tiers**. Q2 and Q1 are special **different-class reconciliation** tiers requiring concordant aligned ALT sequences; they are not generic geometry-only matches.

| Quality | Maximum ΔPOS | Minimum SVLEN similarity | Required conditions |
| --- | ---: | ---: | --- |
| **Q10** | 5 bp | 98% | Same class, event and completeness |
| **Q9** | 10 bp | 95% | Same class, event and completeness |
| **Q8** | 25 bp | 90% | Same class and event; completeness may differ |
| **Q7** | 50 bp | 90% | Same class; event and completeness may differ |
| **Q6** | 100 bp | 90% | Same class; event and completeness may differ |
| **Q5** | 250 bp | 90% | Same class; event and completeness may differ |
| **Q4** | 500 bp | 85% | Same class; event and completeness may differ |
| **Q3** | 500 bp | 80% | Same class; event and completeness may differ |
| **Q2** | 100 bp | 90% | Different classes permitted **only with concordant ALT sequence alignment** |
| **Q1** | 500 bp | 80% | Different classes permitted **only with concordant ALT sequence alignment** |

**Sequence-evidence requirement:** for Q2 and Q1, concordant ALT sequences must satisfy the configured primary sequence-alignment thresholds (D1 defaults: at least **95% sequence identity** and **90% reciprocal sequence coverage**). Positional and size similarity alone never establish a Q2/Q1 match. Conversely, Q10–Q3 do **not** imply that the corresponding ALT sequences were aligned; when an alignment was not run, sequence identity must be reported as *not evaluated*.


## 2. Configuration parameters

### A. General settings and PRIMARY selection

| Parameter                       | D1 value          | Short English description                                         |
| ------------------------------- | ----------------- | ----------------------------------------------------------------- |
| `manifest`                      | TSV path          | File listing the two samples and their input data.                |
| `outdir`                        | Output path       | Directory where comparison results are saved.                     |
| `compare_id`                    | `TWINS_0039_0040` | Unique identifier assigned to the comparison.                     |
| `a`, `b`                        | Sample IDs        | Identifiers of the two samples being compared.                    |
| `reference_build`               | GRCh38            | Reference genome used for genomic coordinates.                    |
| `cores`                         | 18                | Maximum number of CPU cores used by the pipeline.                 |
| `io_cores`                      | 2                 | Workers available for parallel input loading.                     |
| `clip_cores`                    | 8                 | Workers available for soft-clipping analysis.                     |
| `primary_min_support`           | 3                 | Minimum supporting ALT reads required for a PRIMARY insertion.    |
| `min_svlen`                     | 50 bp             | Minimum insertion length considered in the analysis.              |
| `require_pass_cover`            | TRUE              | Requires a passing coverage-quality flag for PRIMARY selection.   |
| `require_pass_region`           | TRUE              | Requires a passing regional-quality flag.                         |
| `require_resolved_class`        | TRUE              | Requires an identified mobile-element class.                      |
| `require_resolved_event`        | TRUE              | Requires a resolved insertion-event annotation.                   |
| `require_resolved_completeness` | TRUE              | Requires a resolved completeness annotation.                      |
| `primary_window`                | 500 bp            | Maximum breakpoint-search distance for PRIMARY matching.          |
| `primary_min_size_sim`          | 0.80              | Minimum insertion-size similarity considered for PRIMARY pairing. |

### B. PRIMARY sequence comparison

| Parameter                    | D1 value    | Short English description                                                                      |
| ---------------------------- | ----------- | ---------------------------------------------------------------------------------------------- |
| `primary_seq_mode`           | `conflicts` | Aligns ALT sequences selectively when structural annotations conflict.                         |
| `primary_seq_identity`       | 95%         | Minimum sequence identity required for an accepted ALT alignment.                              |
| `primary_seq_coverage`       | 90%         | Minimum reciprocal sequence coverage required for ALT agreement.                               |
| `primary_seq_veto`           | TRUE        | Rejects a candidate PRIMARY match when an evaluated ALT alignment fails the required criteria. |
| `primary_seq_max_alignments` | 5,000       | Maximum number of PRIMARY sequence alignments permitted.                                       |
| `primary_seq_max_cells`      | 80 million  | Maximum alignment-matrix size allowed for an individual sequence comparison.                   |
| `q10_pos`–`q3_pos`           | 5–500 bp    | Explicit same-class Q10–Q3 position thresholds; see the Q2/Q1 ALT-aligned exception rules above.                                    |
| `q10_size`–`q3_size`         | 0.98–0.80   | Explicit same-class Q10–Q3 size thresholds; Q2/Q1 rules are described separately above.                                |

### C. Rescue and undermerge

| Parameter                 | D1 value | Short English description                                                       |
| ------------------------- | -------- | ------------------------------------------------------------------------------- |
| `rescue_window`           | 1,000 bp | Maximum distance used to search for evidence supporting an unmatched insertion. |
| `rescue_size_sim`         | 0.80     | Minimum insertion-size similarity for rescue candidates.                        |
| `undermerge_seq_identity` | 95%      | Minimum ALT sequence identity for linking displaced PRIMARY insertions.         |
| `undermerge_seq_coverage` | 90%      | Minimum reciprocal ALT coverage for a displaced PRIMARY link.                   |

Undermerge means that two reconstructions may represent the same insertion but were not paired during the initial PRIMARY matching because their reported positions differ too much.

### D. Duplicate guard and nearby-locus protection

| Parameter                             | D1 value | Short English description                                                                     |
| ------------------------------------- | -------- | --------------------------------------------------------------------------------------------- |
| `duplicate_guard_enabled`             | TRUE     | Prevents potentially redundant insertion evidence from being counted as an additional rescue. |
| `duplicate_window`                    | 1,000 bp | Maximum distance for searching related PRIMARY representations within the same sample.        |
| `duplicate_geometry_window`           | 100 bp   | Maximum distance for identifying particularly close, structurally similar representations.    |
| `duplicate_geometry_size_sim`         | 0.85     | Minimum size similarity required for a strong geometry-based link.                            |
| `duplicate_align_enabled`             | TRUE     | Enables selective ALT sequence comparison between nearby representations.                     |
| `duplicate_align_min_size_sim`        | 0.80     | Minimum size similarity required before attempting duplicate-related sequence alignment.      |
| `duplicate_seq_identity`              | 95%      | Minimum ALT sequence identity for a sequence-supported duplicate link.                        |
| `duplicate_seq_coverage`              | 90%      | Minimum reciprocal ALT coverage for a sequence-supported duplicate link.                      |
| `duplicate_max_alignments_per_source` | 4        | Maximum duplicate-related sequence alignments attempted per unmatched PRIMARY.                |

A nearby insertion is not automatically considered a duplicate. The guard looks for additional evidence, such as a shared original event identifier, concordant ALT sequence or strong structural similarity.

### E. VCF rescue

| Parameter                     | D1 value | Short English description                                                        |
| ----------------------------- | -------- | -------------------------------------------------------------------------------- |
| `vcf_window`                  | 500 bp   | Maximum distance for searching nearby VCF insertion calls.                       |
| `vcf_size_sim`                | 0.80     | Minimum size similarity between a VCF insertion and the source PRIMARY.          |
| `vcf_min_support`             | 2 reads  | Minimum ALT read support used by VCF evidence checks.                            |
| `vcf_require_pass`            | FALSE    | Determines whether only VCF records with a passing FILTER status are considered. |
| `vcf_allow_single_alt_caller` | TRUE     | Allows rescue using one compatible VCF caller when the other checks permit it.   |

RetroAtlas interrogates cuteSV, Severus and Sniffles2. Multiple supporting callers provide stronger corroboration than a single caller, but the pipeline also evaluates genotype and ALT sequence consistency.

### F. Soft-clipping rescue

| Parameter           | D1 value  | Short English description                                                   |
| ------------------- | --------- | --------------------------------------------------------------------------- |
| `clip_enabled`      | TRUE      | Enables rescue using soft-clipped sequencing reads.                         |
| `clip_input_mode`   | `preload` | Loads soft-clipping data before evaluating individual loci.                 |
| `clip_window`       | 5,000 bp  | Broad genomic window used to retrieve nearby soft-clipping evidence.        |
| `clip_core_window`  | 1,000 bp  | Main breakpoint window used for local soft-clipping evaluation.             |
| `clip_min_identity` | 90%       | Minimum sequence identity required for a compatible soft-clipped alignment. |
| `clip_min_coverage` | 90%       | Minimum required sequence coverage for soft-clipping evidence.              |
| `clip_min_clips`    | 3         | Minimum number of supporting soft-clipped reads.                            |
| `clip_min_cluster`  | 3         | Minimum number of clips supporting a breakpoint cluster.                    |
| `clip_min_union_bp` | 200 bp    | Minimum combined sequence span covered by supporting clips.                 |

Soft-clipping evidence can establish plausible insertion presence without determining whether the sample is heterozygous or homozygous.

### G. Fragment rescue

| Parameter                 | D1 value | Short English description                                                   |
| ------------------------- | -------- | --------------------------------------------------------------------------- |
| `fragment_enabled`        | TRUE     | Enables detection of insertions supported by multiple structural fragments. |
| `fragment_window`         | 1,000 bp | Maximum search distance for nearby fragments.                               |
| `fragment_max_components` | 4        | Maximum number of fragments considered in a candidate reconstruction.       |
| `fragment_max_span`       | 1,000 bp | Maximum genomic span allowed for a group of fragments.                      |
| `fragment_size_sim`       | 0.80     | Minimum size similarity required for a fragment-based rescue.               |

### H. Callability, sex chromosomes and outputs

| Parameter                          | D1 value   | Short English description                                                                    |
| ---------------------------------- | ---------- | -------------------------------------------------------------------------------------------- |
| `callable_max_depth_ratio`         | 3          | Maximum permitted regional depth relative to the expected sample depth.                      |
| `require_all_vcf_for_differential` | TRUE       | Requires all three configured VCF sources to be available for the strict differential check. |
| `include_sex_chromosomes`          | TRUE       | Includes eligible sex-chromosome regions using sex-aware rules.                              |
| `write_tsv`                        | TRUE       | Exports human-readable TSV results in addition to the RDS.                                   |
| `softclip_module`                  | R filename | R module containing the soft-clipping rescue functions.                                      |
| `sequence_module`                  | R filename | R module containing ALT sequence-comparison functions.                                       |
| `sex_module`                       | R filename | R module containing sex-aware ploidy functions.                                              |

## 3. REVIEW subtypes — useful output descriptors

These are particularly valuable for describing the improvements in D1.

| Descriptor                                          | Short English definition                                                                                                                 |
| --------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------- |
| `POSSIBLE_MULTIALLELIC_LOCUS_ALREADY_MATCHED`       | A nearby representation from the same sample is already matched to the other sample and may belong to the same biological locus.         |
| `POSSIBLE_MULTIALLELIC_OR_REDUNDANT_SOURCE_PRIMARY` | Nearby PRIMARY reconstructions in the same sample may represent different alleles, fragments, or redundant representations of one event. |
| `SECONDARY_EVIDENCE_AMBIGUOUS`                      | VCF, soft-clipping, or fragment evidence is present but insufficient or contradictory.                                                   |
| `OPPOSITE_PRIMARY_WITHIN_500BP`                     | A PRIMARY exists nearby in the opposite sample but could not be confidently paired.                                                      |
| `PLAUSIBLE_INTRASAMPLE_REPRESENTATION_LINK`         | A potential relationship between PRIMARY reconstructions within one sample prevents an automatic additional rescue.                      |
