#!/usr/bin/env Rscript

# Unified, paired-permutation calibration for the five additional trajectories
# reported in JHG-26-335, plus the five-seed ATTR stability analysis.
#
# This script uses already-derived endpoint-level MR results. It does not call
# OpenGWAS and does not rerun the underlying PWMR analyses.

project_dir <- normalizePath(getwd())
input_dir <- file.path(project_dir, "inputs")
output_dir <- file.path(project_dir, "rerun_outputs", "calibration")
if (!dir.exists(input_dir)) {
  stop("Run from the revision project root containing the inputs directory.")
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

primary_methods <- c("Inverse variance weighted", "Wald ratio")
n_perm <- 10000L
main_seed <- 42L

clean_endpoint <- function(filename, endpoint) {
  x <- read.csv(file.path(input_dir, filename), stringsAsFactors = FALSE,
                check.names = FALSE)
  x <- x[x$method %in% primary_methods, , drop = FALSE]
  if (anyDuplicated(x$id.exposure)) {
    stop(filename, " contains duplicate primary-method records.")
  }
  protein_name <- if ("protein_name" %in% names(x)) x$protein_name else x$exposure
  protein_name <- sub(" \\|\\| id:.*$", "", protein_name)
  out <- data.frame(
    id.exposure = as.character(x$id.exposure),
    protein_name = protein_name,
    beta = x$b,
    se = x$se,
    p = x$pval,
    method = x$method,
    nsnp = x$nsnp,
    stringsAsFactors = FALSE
  )
  names(out)[3:7] <- paste0(c("beta_", "se_", "p_", "method_", "nsnp_"), endpoint)
  out
}

lead_snp_map <- function() {
  iv <- readRDS(file.path(input_dir, "01_extraction_progress.rds"))
  iv <- iv[order(iv$id.exposure, iv$pval.exposure, na.last = TRUE),
           c("id.exposure", "SNP", "pval.exposure")]
  iv <- iv[!duplicated(iv$id.exposure), , drop = FALSE]
  iv[, c("id.exposure", "SNP")]
}

add_locus <- function(x) {
  x <- merge(x, lead_snp_map(), by = "id.exposure", all.x = TRUE, sort = FALSE)
  x$locus_id <- ifelse(is.na(x$SNP) | x$SNP == "", x$id.exposure, x$SNP)

  # The submitted analysis treated CA5A, TNS4, and STX10 as one PON1-region
  # locus because their instruments are identical or in complete LD.
  pon1_ids <- c("prot-a-332", "prot-a-3070", "prot-a-2883")
  x$locus_id[x$id.exposure %in% pon1_ids] <- "PON1_7q21"
  x
}

merge_trajectory <- function(endpoint_data, endpoints) {
  x <- endpoint_data[[endpoints[1L]]]
  for (endpoint in endpoints[-1L]) {
    y <- endpoint_data[[endpoint]]
    y <- y[, setdiff(names(y), "protein_name"), drop = FALSE]
    x <- merge(x, y, by = "id.exposure", sort = FALSE)
  }
  add_locus(x)
}

concordant <- function(x, endpoints, cutoff = 0.05) {
  p1 <- x[[paste0("p_", endpoints[1L])]]
  p2 <- x[[paste0("p_", endpoints[2L])]]
  p3 <- x[[paste0("p_", endpoints[3L])]]
  b1 <- x[[paste0("beta_", endpoints[1L])]]
  b2 <- x[[paste0("beta_", endpoints[2L])]]
  b3 <- x[[paste0("beta_", endpoints[3L])]]
  is.finite(p1) & is.finite(p2) & is.finite(p3) &
    is.finite(b1) & is.finite(b2) & is.finite(b3) &
    p1 < cutoff & p2 < cutoff & p3 < cutoff &
    sign(b1) != 0 & sign(b1) == sign(b2) & sign(b1) == sign(b3)
}

paired_permutation <- function(x, endpoints, fixed_endpoint,
                               cutoff = 0.05, seed = 42L,
                               n_iterations = 10000L) {
  n <- nrow(x)
  observed_hit <- concordant(x, endpoints, cutoff)
  observed_proteins <- sum(observed_hit)
  observed_loci <- length(unique(x$locus_id[observed_hit]))

  p_values <- lapply(endpoints, function(e) x[[paste0("p_", e)]])
  beta_values <- lapply(endpoints, function(e) x[[paste0("beta_", e)]])
  null <- integer(n_iterations)
  set.seed(seed)

  for (i in seq_len(n_iterations)) {
    indices <- lapply(endpoints, function(endpoint) {
      if (identical(endpoint, fixed_endpoint)) seq_len(n) else sample.int(n)
    })
    p1 <- p_values[[1L]][indices[[1L]]]
    p2 <- p_values[[2L]][indices[[2L]]]
    p3 <- p_values[[3L]][indices[[3L]]]
    b1 <- beta_values[[1L]][indices[[1L]]]
    b2 <- beta_values[[2L]][indices[[2L]]]
    b3 <- beta_values[[3L]][indices[[3L]]]

    # A single index is used for beta and P within each shuffled endpoint.
    # This preserves each MR result as one statistical unit.
    hit <- is.finite(p1) & is.finite(p2) & is.finite(p3) &
      is.finite(b1) & is.finite(b2) & is.finite(b3) &
      p1 < cutoff & p2 < cutoff & p3 < cutoff &
      sign(b1) != 0 & sign(b1) == sign(b2) & sign(b1) == sign(b3)
    null[i] <- length(unique(x$locus_id[hit]))
  }

  result <- data.frame(
    observed_proteins = observed_proteins,
    observed_loci = observed_loci,
    null_mean_loci = mean(null),
    null_median_loci = median(null),
    null_max_loci = max(null),
    null_ge_observed = sum(null >= observed_loci),
    empirical_p = (sum(null >= observed_loci) + 1) / (n_iterations + 1),
    n_proteins_tested = n,
    n_unique_loci = length(unique(x$locus_id)),
    cutoff = cutoff,
    fixed_endpoint = fixed_endpoint,
    n_permutations = n_iterations,
    seed = seed,
    stringsAsFactors = FALSE
  )
  list(result = result, null = null, hit = observed_hit)
}

endpoint_files <- c(
  bmi = "PosCtrl_BMI_all_results.csv",
  t2d = "19_mr_t2d_results.csv",
  hf = "10_mr_hf_results.csv",
  ldl = "30_mr_ldl_results.csv",
  cad = "31_mr_cad_results.csv",
  ckd = "35_mr_ckd_results.csv",
  sbp = "36_mr_sbp_results.csv",
  stroke = "37_mr_stroke_results.csv",
  masld = "39_mr_masld_results.csv"
)
endpoint_data <- lapply(names(endpoint_files), function(endpoint) {
  clean_endpoint(endpoint_files[[endpoint]], endpoint)
})
names(endpoint_data) <- names(endpoint_files)

trajectory_definitions <- list(
  BMI_T2D_CKD = c("bmi", "t2d", "ckd"),
  BMI_CAD_HF = c("bmi", "cad", "hf"),
  LDL_CAD_HF = c("ldl", "cad", "hf"),
  BMI_SBP_Stroke = c("bmi", "sbp", "stroke"),
  SBP_MASLD_HF = c("sbp", "masld", "hf")
)

trajectory_labels <- c(
  BMI_T2D_CKD = "BMI -> T2D -> CKD",
  BMI_CAD_HF = "BMI -> CAD -> HF",
  LDL_CAD_HF = "LDL -> CAD -> HF",
  BMI_SBP_Stroke = "BMI -> SBP -> Stroke",
  SBP_MASLD_HF = "SBP -> MASLD/NAFLD -> HF"
)

calibration_rows <- list()
membership_rows <- list()
null_distributions <- list()

for (trajectory_id in names(trajectory_definitions)) {
  endpoints <- trajectory_definitions[[trajectory_id]]
  x <- merge_trajectory(endpoint_data, endpoints)
  test <- paired_permutation(
    x = x,
    endpoints = endpoints,
    fixed_endpoint = endpoints[2L],
    cutoff = 0.05,
    seed = main_seed,
    n_iterations = n_perm
  )
  row <- test$result
  row$trajectory_id <- trajectory_id
  row$trajectory <- unname(trajectory_labels[[trajectory_id]])
  row$node_1 <- endpoints[1L]
  row$node_2 <- endpoints[2L]
  row$node_3 <- endpoints[3L]
  calibration_rows[[trajectory_id]] <- row
  null_distributions[[trajectory_id]] <- test$null

  if (any(test$hit)) {
    h <- x[test$hit, , drop = FALSE]
    membership_rows[[trajectory_id]] <- data.frame(
      trajectory_id = trajectory_id,
      trajectory = unname(trajectory_labels[[trajectory_id]]),
      id.exposure = h$id.exposure,
      protein_name = h$protein_name,
      locus_id = h$locus_id,
      stringsAsFactors = FALSE
    )
  }
}

calibration <- do.call(rbind, calibration_rows)
calibration <- calibration[, c(
  "trajectory_id", "trajectory", "node_1", "node_2", "node_3",
  setdiff(names(calibration), c("trajectory_id", "trajectory", "node_1", "node_2", "node_3"))
)]
rownames(calibration) <- NULL
write.csv(calibration,
          file.path(output_dir, "multi_trajectory_calibration_corrected.csv"),
          row.names = FALSE)

if (length(membership_rows) > 0L) {
  membership <- do.call(rbind, membership_rows)
  rownames(membership) <- NULL
} else {
  membership <- data.frame(
    trajectory_id = character(), trajectory = character(),
    id.exposure = character(), protein_name = character(),
    locus_id = character(), stringsAsFactors = FALSE
  )
}
write.csv(membership,
          file.path(output_dir, "multi_trajectory_membership_corrected.csv"),
          row.names = FALSE)
saveRDS(null_distributions,
        file.path(output_dir, "multi_trajectory_null_distributions_corrected.rds"))

# Compare against the exact values reported in the submitted manuscript.
submitted <- data.frame(
  trajectory_id = c("BMI_T2D_CKD", "BMI_CAD_HF", "LDL_CAD_HF",
                    "BMI_SBP_Stroke", "SBP_MASLD_HF"),
  submitted_loci = c(2L, 2L, 1L, 1L, 2L),
  submitted_empirical_p = c(0.17, 0.47, 0.74, 0.91, 0.034),
  stringsAsFactors = FALSE
)
comparison <- merge(
  submitted,
  calibration[, c("trajectory_id", "observed_proteins", "observed_loci", "empirical_p")],
  by = "trajectory_id", sort = FALSE
)
comparison$loci_match <- comparison$submitted_loci == comparison$observed_loci
comparison$absolute_p_difference <- abs(
  comparison$submitted_empirical_p - comparison$empirical_p
)
write.csv(comparison,
          file.path(output_dir, "multi_trajectory_submitted_vs_corrected.csv"),
          row.names = FALSE)

# ATTR five-seed stability using the same paired-permutation unit definition.
attr <- read.csv(file.path(input_dir, "11_triangulation_combined.csv"),
                 stringsAsFactors = FALSE, check.names = FALSE)
names(attr)[names(attr) == "beta_ca"] <- "beta_amy"
names(attr)[names(attr) == "se_ca"] <- "se_amy"
names(attr)[names(attr) == "pval_ca"] <- "p_amy"
names(attr)[names(attr) == "pval_cts"] <- "p_cts"
names(attr)[names(attr) == "pval_hf"] <- "p_hf"
attr <- add_locus(attr)

attr_seeds <- c(11L, 42L, 101L, 2024L, 31415L)
attr_rows <- list()
attr_nulls <- list()
for (seed in attr_seeds) {
  test <- paired_permutation(
    x = attr,
    endpoints = c("cts", "amy", "hf"),
    fixed_endpoint = "amy",
    cutoff = 0.05,
    seed = seed,
    n_iterations = n_perm
  )
  attr_rows[[as.character(seed)]] <- test$result
  attr_nulls[[as.character(seed)]] <- test$null
}
attr_stability <- do.call(rbind, attr_rows)
rownames(attr_stability) <- NULL
write.csv(attr_stability,
          file.path(output_dir, "attr_permutation_multiseed_corrected.csv"),
          row.names = FALSE)
saveRDS(attr_nulls,
        file.path(output_dir, "attr_permutation_multiseed_nulls_corrected.rds"))

existing_stability_path <- file.path(project_dir, "source_baseline", "23_permutation_stability.csv")
if (file.exists(existing_stability_path)) {
  existing <- read.csv(existing_stability_path, stringsAsFactors = FALSE)
  attr_compare <- merge(
    existing[, c("seed", "observed_hits", "null_mean", "null_max", "empirical_p")],
    attr_stability[, c("seed", "observed_loci", "null_mean_loci", "null_max_loci", "empirical_p")],
    by = "seed", suffixes = c("_existing", "_corrected")
  )
  write.csv(attr_compare,
            file.path(output_dir, "attr_multiseed_existing_vs_corrected.csv"),
            row.names = FALSE)
}

cat("Corrected multi-trajectory calibration\n")
print(calibration[, c("trajectory", "observed_proteins", "observed_loci", "empirical_p")],
      row.names = FALSE)
cat("\nATTR five-seed empirical P range: ",
    format(min(attr_stability$empirical_p), digits = 6), " to ",
    format(max(attr_stability$empirical_p), digits = 6), "\n", sep = "")
