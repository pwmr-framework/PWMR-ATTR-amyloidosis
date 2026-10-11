#!/usr/bin/env Rscript

# Reviewer-requested sensitivity analyses for JHG-26-335.
#
# Principles:
#   1. Preserve beta/P-value pairing within every endpoint during permutation.
#   2. Collapse proteins sharing a lead SNP to a single locus.
#   3. Keep the original concordance rule explicit and auditable.
#   4. Report effect estimates, uncertainty, and power diagnostics separately
#      because BMI and binary disease endpoints are not on a common beta scale.

# Run from the revision project root. On macOS, Rscript can encode spaces in
# its internal --file argument, so deriving the project path from that value is
# not reliable for this directory name.
project_dir <- normalizePath(getwd())
input_dir <- file.path(project_dir, "inputs")
output_dir <- file.path(project_dir, "rerun_outputs", "calibration")
if (!dir.exists(input_dir)) {
  stop("Run from the revision project root containing the inputs directory.")
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

primary_methods <- c("Inverse variance weighted", "Wald ratio")
cutoffs <- c(0.01, 0.025, 0.05, 0.10)
z_975 <- qnorm(0.975)
z_80 <- qnorm(0.80)

clean_primary <- function(filename, endpoint) {
  x <- read.csv(file.path(input_dir, filename), stringsAsFactors = FALSE,
                check.names = FALSE)
  x <- x[x$method %in% primary_methods, , drop = FALSE]
  if (anyDuplicated(x$id.exposure)) {
    stop(filename, " contains duplicate primary-method rows per protein.")
  }
  protein <- if ("protein_name" %in% names(x)) x$protein_name else x$exposure
  protein <- sub(" \\|\\| id:.*$", "", protein)
  out <- data.frame(
    id.exposure = x$id.exposure,
    protein_name = protein,
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

add_locus <- function(x, positive_control = FALSE) {
  x <- merge(x, lead_snp_map(), by = "id.exposure", all.x = TRUE, sort = FALSE)
  x$locus_id <- ifelse(is.na(x$SNP) | x$SNP == "", x$id.exposure, x$SNP)
  if (positive_control) {
    pon1_ids <- c("prot-a-332", "prot-a-3070", "prot-a-2883")
    x$locus_id[x$id.exposure %in% pon1_ids] <- "PON1_7q21"
  }
  x
}

is_concordant <- function(x, endpoints, cutoff) {
  pmat <- vapply(endpoints, function(e) x[[paste0("p_", e)]], numeric(nrow(x)))
  bmat <- vapply(endpoints, function(e) x[[paste0("beta_", e)]], numeric(nrow(x)))
  complete <- rowSums(!is.finite(pmat) | !is.finite(bmat)) == 0L
  significant <- rowSums(pmat < cutoff) == length(endpoints)
  bsign <- sign(bmat)
  same_sign <- bsign[, 1L] != 0 & bsign[, 1L] == bsign[, 2L] &
    bsign[, 1L] == bsign[, 3L]
  complete & significant & same_sign
}

count_hits <- function(x, endpoints, cutoff) {
  hit <- is_concordant(x, endpoints, cutoff)
  c(proteins = sum(hit), loci = length(unique(x$locus_id[hit])))
}

permutation_test <- function(x, endpoints, cutoff, fixed_endpoint,
                             n_perm = 10000L, seed = 42L) {
  n <- nrow(x)
  observed <- count_hits(x, endpoints, cutoff)
  null <- integer(n_perm)
  p_values <- lapply(endpoints, function(e) x[[paste0("p_", e)]])
  beta_values <- lapply(endpoints, function(e) x[[paste0("beta_", e)]])
  set.seed(seed)
  for (i in seq_len(n_perm)) {
    indices <- lapply(endpoints, function(endpoint) {
      if (identical(endpoint, fixed_endpoint)) seq_len(n) else sample.int(n)
    })
    # The same index vector is deliberately used for beta and P within an
    # endpoint, preserving each MR estimate as one statistical unit.
    p1 <- p_values[[1L]][indices[[1L]]]
    p2 <- p_values[[2L]][indices[[2L]]]
    p3 <- p_values[[3L]][indices[[3L]]]
    b1 <- beta_values[[1L]][indices[[1L]]]
    b2 <- beta_values[[2L]][indices[[2L]]]
    b3 <- beta_values[[3L]][indices[[3L]]]
    hit <- is.finite(p1) & is.finite(p2) & is.finite(p3) &
      is.finite(b1) & is.finite(b2) & is.finite(b3) &
      p1 < cutoff & p2 < cutoff & p3 < cutoff &
      sign(b1) != 0 & sign(b1) == sign(b2) & sign(b1) == sign(b3)
    null[i] <- length(unique(x$locus_id[hit]))
  }
  data.frame(
    cutoff = cutoff,
    fixed_endpoint = ifelse(is.na(fixed_endpoint), "none", fixed_endpoint),
    observed_proteins = unname(observed["proteins"]),
    observed_loci = unname(observed["loci"]),
    null_mean_loci = mean(null),
    null_median_loci = median(null),
    null_max_loci = max(null),
    null_ge_observed = sum(null >= observed["loci"]),
    empirical_p = (sum(null >= observed["loci"]) + 1) / (n_perm + 1),
    n_proteins_tested = n,
    n_unique_loci = length(unique(x$locus_id)),
    n_permutations = n_perm,
    seed = seed,
    stringsAsFactors = FALSE
  )
}

# For the fully shuffled scheme, use a sentinel fixed endpoint not present in
# the endpoint vector. This keeps the same paired-shuffling implementation.
permutation_test_all_shuffled <- function(x, endpoints, cutoff,
                                          n_perm = 10000L, seed = 42L) {
  out <- permutation_test(x, endpoints, cutoff, fixed_endpoint = "__none__",
                          n_perm = n_perm, seed = seed)
  out$fixed_endpoint <- "none"
  out
}

make_positive_control <- function() {
  bmi <- clean_primary("PosCtrl_BMI_all_results.csv", "bmi")
  t2d <- clean_primary("19_mr_t2d_results.csv", "t2d")
  hf <- clean_primary("10_mr_hf_results.csv", "hf")
  t2d <- t2d[, setdiff(names(t2d), "protein_name")]
  hf <- hf[, setdiff(names(hf), "protein_name")]
  x <- merge(merge(bmi, t2d, by = "id.exposure"), hf,
             by = "id.exposure", sort = FALSE)
  add_locus(x, positive_control = TRUE)
}

make_attr <- function() {
  x <- read.csv(file.path(input_dir, "11_triangulation_combined.csv"),
                stringsAsFactors = FALSE, check.names = FALSE)
  names(x)[names(x) == "beta_ca"] <- "beta_amy"
  names(x)[names(x) == "se_ca"] <- "se_amy"
  names(x)[names(x) == "pval_ca"] <- "p_amy"
  names(x)[names(x) == "pval_cts"] <- "p_cts"
  names(x)[names(x) == "pval_hf"] <- "p_hf"
  add_locus(x, positive_control = FALSE)
}

pc <- make_positive_control()
attr <- make_attr()
pc_endpoints <- c("bmi", "t2d", "hf")
attr_endpoints <- c("cts", "amy", "hf")

# 1. Exact invariance of the concordance filter to endpoint ordering.
all_orders <- function(v) {
  rbind(
    v[c(1, 2, 3)], v[c(1, 3, 2)], v[c(2, 1, 3)],
    v[c(2, 3, 1)], v[c(3, 1, 2)], v[c(3, 2, 1)]
  )
}
ordering_rows <- list()
for (label in c("positive_control", "ATTR")) {
  x <- if (label == "positive_control") pc else attr
  endpoints <- if (label == "positive_control") pc_endpoints else attr_endpoints
  for (ord in seq_len(nrow(all_orders(endpoints)))) {
    endpoint_order <- all_orders(endpoints)[ord, ]
    hit <- is_concordant(x, endpoint_order, 0.05)
    ordering_rows[[length(ordering_rows) + 1L]] <- data.frame(
      analysis = label,
      order = paste(endpoint_order, collapse = " -> "),
      protein_hits = sum(hit),
      locus_hits = length(unique(x$locus_id[hit])),
      hit_ids = paste(sort(x$id.exposure[hit]), collapse = ";"),
      stringsAsFactors = FALSE
    )
  }
}
ordering <- do.call(rbind, ordering_rows)
write.csv(ordering, file.path(output_dir, "ordering_invariance.csv"), row.names = FALSE)

# 2. Threshold sensitivity for the full positive-control and ATTR analyses.
threshold_rows <- list()
for (cutoff in cutoffs) {
  a <- permutation_test(pc, pc_endpoints, cutoff, "t2d", 10000L, 42L)
  a$analysis <- "positive_control"
  threshold_rows[[length(threshold_rows) + 1L]] <- a
  b <- permutation_test(attr, attr_endpoints, cutoff, "amy", 10000L, 42L)
  b$analysis <- "ATTR"
  threshold_rows[[length(threshold_rows) + 1L]] <- b
}
threshold_results <- do.call(rbind, threshold_rows)
threshold_results <- threshold_results[, c("analysis", setdiff(names(threshold_results), "analysis"))]
write.csv(threshold_results,
          file.path(output_dir, "threshold_sensitivity_full.csv"), row.names = FALSE)

# 3. Correctly paired anchor-sensitivity permutations at P < 0.05.
anchor_rows <- list()
for (label in c("positive_control", "ATTR")) {
  x <- if (label == "positive_control") pc else attr
  endpoints <- if (label == "positive_control") pc_endpoints else attr_endpoints
  for (anchor in endpoints) {
    z <- permutation_test(x, endpoints, 0.05, anchor, 10000L, 42L)
    z$analysis <- label
    anchor_rows[[length(anchor_rows) + 1L]] <- z
  }
  z <- permutation_test_all_shuffled(x, endpoints, 0.05, 10000L, 42L)
  z$analysis <- label
  anchor_rows[[length(anchor_rows) + 1L]] <- z
}
anchor_results <- do.call(rbind, anchor_rows)
anchor_results <- anchor_results[, c("analysis", setdiff(names(anchor_results), "analysis"))]
write.csv(anchor_results,
          file.path(output_dir, "paired_anchor_sensitivity.csv"), row.names = FALSE)

# 4. Dense power-degradation grid crossed with all requested P thresholds.
# T2D and HF are degraded exactly as in the submitted analysis; BMI is fixed.
fractions <- sort(unique(c(0.005, 0.01, 0.02, 0.05, 0.10, 0.20,
                           seq(0.50, 1.00, by = 0.05))))
hf_full_neff <- 180076
power_rows <- list()
for (fraction in fractions) {
  degraded <- pc
  for (endpoint in c("t2d", "hf")) {
    degraded[[paste0("se_", endpoint)]] <- pc[[paste0("se_", endpoint)]] / sqrt(fraction)
    z <- degraded[[paste0("beta_", endpoint)]] /
      degraded[[paste0("se_", endpoint)]]
    degraded[[paste0("p_", endpoint)]] <- 2 * pnorm(-abs(z))
  }
  for (cutoff in cutoffs) {
    z <- permutation_test(degraded, pc_endpoints, cutoff, "t2d", 1000L, 42L)
    z$sample_fraction <- fraction
    z$weakest_endpoint_neff <- hf_full_neff * fraction
    power_rows[[length(power_rows) + 1L]] <- z
  }
}
power_results <- do.call(rbind, power_rows)
power_results <- power_results[, c("sample_fraction", "weakest_endpoint_neff",
                                   setdiff(names(power_results),
                                           c("sample_fraction", "weakest_endpoint_neff")))]
write.csv(power_results,
          file.path(output_dir, "power_threshold_sensitivity.csv"), row.names = FALSE)

# 5. Effect-size and uncertainty table for concordant hits.
effect_table <- function(x, endpoints, analysis) {
  hit <- is_concordant(x, endpoints, 0.05)
  h <- x[hit, , drop = FALSE]
  rows <- list()
  for (i in seq_len(nrow(h))) {
    for (endpoint in endpoints) {
      beta <- h[[paste0("beta_", endpoint)]][i]
      se <- h[[paste0("se_", endpoint)]][i]
      p <- h[[paste0("p_", endpoint)]][i]
      # Keep the actual source-GWAS scale; a binary phenotype need not use log-odds.
      log_odds_outcome <- endpoint %in% c("t2d", "amy", "hf")
      effect_scale <- if (endpoint == "bmi") {
        "continuous BMI units"
      } else if (endpoint %in% c("cts", "asthma")) {
        "Source linear-regression scale"
      } else if (log_odds_outcome) {
        "log odds ratio"
      } else {
        stop("Unverified outcome scale: ", endpoint)
      }
      rows[[length(rows) + 1L]] <- data.frame(
        analysis = analysis,
        id.exposure = h$id.exposure[i],
        protein_name = h$protein_name[i],
        locus_id = h$locus_id[i],
        endpoint = endpoint,
        beta = beta,
        se = se,
        ci_lower = beta - z_975 * se,
        ci_upper = beta + z_975 * se,
        z = beta / se,
        p = p,
        effect_scale = effect_scale,
        odds_ratio = ifelse(log_odds_outcome, exp(beta), NA_real_),
        or_ci_lower = ifelse(log_odds_outcome, exp(beta - z_975 * se), NA_real_),
        or_ci_upper = ifelse(log_odds_outcome, exp(beta + z_975 * se), NA_real_),
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}
effects <- rbind(effect_table(pc, pc_endpoints, "positive_control"),
                 effect_table(attr, attr_endpoints, "ATTR"))
write.csv(effects, file.path(output_dir, "concordant_effects_with_ci.csv"), row.names = FALSE)

# 6. For each positive-control hit, calculate the retained sample fraction at
# which T2D and HF would both remain nominally significant under fixed betas.
pc_hit <- pc[is_concordant(pc, pc_endpoints, 0.05), , drop = FALSE]
critical_fraction <- function(beta, se, cutoff = 0.05) {
  zcrit <- qnorm(1 - cutoff / 2)
  (zcrit / abs(beta / se))^2
}
emergence <- data.frame(
  id.exposure = pc_hit$id.exposure,
  protein_name = pc_hit$protein_name,
  locus_id = pc_hit$locus_id,
  t2d_critical_fraction = critical_fraction(pc_hit$beta_t2d, pc_hit$se_t2d),
  hf_critical_fraction = critical_fraction(pc_hit$beta_hf, pc_hit$se_hf),
  joint_critical_fraction = pmax(
    critical_fraction(pc_hit$beta_t2d, pc_hit$se_t2d),
    critical_fraction(pc_hit$beta_hf, pc_hit$se_hf)
  ),
  stringsAsFactors = FALSE
)
emergence$joint_critical_neff <- emergence$joint_critical_fraction * hf_full_neff
write.csv(emergence, file.path(output_dir, "positive_control_emergence_fractions.csv"),
          row.names = FALSE)

# 7. Fill the MDE reporting gap in the submitted Supplementary Methods.
amy <- clean_primary("04_mr_results_corrected.csv", "amy")
amy$mde_abs_beta_80pct_power_alpha_0_05 <- (z_975 + z_80) * amy$se_amy
amy$minimum_detectable_or_above_1 <- exp(amy$mde_abs_beta_80pct_power_alpha_0_05)
amy$maximum_detectable_or_below_1 <- exp(-amy$mde_abs_beta_80pct_power_alpha_0_05)
amy$observed_ci_lower <- amy$beta_amy - z_975 * amy$se_amy
amy$observed_ci_upper <- amy$beta_amy + z_975 * amy$se_amy
write.csv(amy, file.path(output_dir, "amyloidosis_mde_80pct.csv"), row.names = FALSE)

# 8. Machine-readable run summary and basic assertions.
stopifnot(all(vapply(split(ordering, ordering$analysis), function(z) {
  length(unique(z$protein_hits)) == 1L &&
    length(unique(z$locus_hits)) == 1L &&
    length(unique(z$hit_ids)) == 1L
}, logical(1))))
stopifnot(nrow(pc) == 1116L)

summary_lines <- c(
  paste("Run time:", format(Sys.time(), tz = "Asia/Shanghai")),
  paste("R version:", R.version.string),
  paste("Positive-control proteins/loci:", nrow(pc), "/", length(unique(pc$locus_id))),
  paste("ATTR proteins/loci:", nrow(attr), "/", length(unique(attr$locus_id))),
  paste("Amyloidosis MDE rows:", nrow(amy)),
  "All permutations preserved beta/P pairing within endpoint; SE was not shuffled.",
  "All six endpoint orderings returned identical concordant memberships."
)
writeLines(summary_lines, file.path(output_dir, "run_summary.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")
