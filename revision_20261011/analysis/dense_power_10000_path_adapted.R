# Path-adapted copy of the recovered August 2026 dense-grid driver.
# Only input/output setup differs. This copy was not executed to generate results.
# Run from the candidate bundle root after input access and provenance review.
project_dir <- normalizePath(getwd())
input_dir <- file.path(project_dir, "inputs")
required_files <- c("01_extraction_progress.rds", "PosCtrl_BMI_all_results.csv",
                    "19_mr_t2d_results.csv", "10_mr_hf_results.csv")
missing_files <- required_files[!file.exists(file.path(input_dir, required_files))]
if (length(missing_files)) stop("Missing inputs: ", paste(missing_files, collapse = ", "))
output_dir <- file.path(project_dir, "rerun_outputs", "dense_grid_10000")
output_file <- file.path(output_dir, "power_degradation_dense_50-100pct_10000perm.csv")
if (file.exists(output_file)) stop("Refusing to overwrite existing output: ", output_file)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

primary_methods <- c("Inverse variance weighted", "Wald ratio")
cutoff <- 0.05
fractions <- seq(0.50, 1.00, by = 0.05)
n_perm <- 10000L
seed <- 42L
hf_full_neff <- 180076

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
    stringsAsFactors = FALSE
  )
  names(out)[3:5] <- paste0(c("beta_", "se_", "p_"), endpoint)
  out
}

iv <- readRDS(file.path(input_dir, "01_extraction_progress.rds"))
iv <- iv[order(iv$id.exposure, iv$pval.exposure, na.last = TRUE),
         c("id.exposure", "SNP", "pval.exposure")]
iv <- iv[!duplicated(iv$id.exposure), c("id.exposure", "SNP")]

bmi <- clean_primary("PosCtrl_BMI_all_results.csv", "bmi")
t2d <- clean_primary("19_mr_t2d_results.csv", "t2d")
hf <- clean_primary("10_mr_hf_results.csv", "hf")
t2d <- t2d[, setdiff(names(t2d), "protein_name")]
hf <- hf[, setdiff(names(hf), "protein_name")]
pc <- merge(merge(bmi, t2d, by = "id.exposure"), hf,
            by = "id.exposure", sort = FALSE)
pc <- merge(pc, iv, by = "id.exposure", all.x = TRUE, sort = FALSE)
pc$locus_id <- ifelse(is.na(pc$SNP) | pc$SNP == "", pc$id.exposure, pc$SNP)
pon1_ids <- c("prot-a-332", "prot-a-3070", "prot-a-2883")
pc$locus_id[pc$id.exposure %in% pon1_ids] <- "PON1_7q21"

endpoints <- c("bmi", "t2d", "hf")

is_concordant <- function(x) {
  pmat <- vapply(endpoints, function(e) x[[paste0("p_", e)]], numeric(nrow(x)))
  bmat <- vapply(endpoints, function(e) x[[paste0("beta_", e)]], numeric(nrow(x)))
  rowSums(!is.finite(pmat) | !is.finite(bmat)) == 0L &
    rowSums(pmat < cutoff) == length(endpoints) &
    sign(bmat[, 1L]) != 0 &
    sign(bmat[, 1L]) == sign(bmat[, 2L]) &
    sign(bmat[, 1L]) == sign(bmat[, 3L])
}

permutation_row <- function(x, fraction) {
  hit_observed <- is_concordant(x)
  observed_proteins <- sum(hit_observed)
  observed_loci <- length(unique(x$locus_id[hit_observed]))

  p_values <- lapply(endpoints, function(e) x[[paste0("p_", e)]])
  beta_values <- lapply(endpoints, function(e) x[[paste0("beta_", e)]])
  n <- nrow(x)
  null <- integer(n_perm)
  set.seed(seed)

  for (i in seq_len(n_perm)) {
    indices <- list(sample.int(n), seq_len(n), sample.int(n))
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

  hit_table <- x[hit_observed,
                 c("id.exposure", "protein_name", "locus_id"), drop = FALSE]
  ge <- sum(null >= observed_loci)
  data.frame(
    sample_fraction = fraction,
    weakest_endpoint_neff = hf_full_neff * fraction,
    cutoff = cutoff,
    fixed_endpoint = "t2d",
    observed_proteins = observed_proteins,
    observed_loci = observed_loci,
    hit_ids = paste(sort(hit_table$id.exposure), collapse = ";"),
    hit_protein_names = paste(hit_table$protein_name[order(hit_table$id.exposure)],
                              collapse = ";"),
    hit_locus_ids = paste(sort(unique(hit_table$locus_id)), collapse = ";"),
    null_mean_loci = mean(null),
    null_median_loci = median(null),
    null_max_loci = max(null),
    null_ge_observed = ge,
    empirical_p = (ge + 1) / (n_perm + 1),
    n_proteins_tested = n,
    n_unique_loci = length(unique(x$locus_id)),
    n_permutations = n_perm,
    seed = seed,
    stringsAsFactors = FALSE
  )
}

start_time <- proc.time()[["elapsed"]]
rows <- vector("list", length(fractions))
for (i in seq_along(fractions)) {
  fraction <- fractions[i]
  degraded <- pc
  for (endpoint in c("t2d", "hf")) {
    degraded[[paste0("se_", endpoint)]] <-
      pc[[paste0("se_", endpoint)]] / sqrt(fraction)
    z <- degraded[[paste0("beta_", endpoint)]] /
      degraded[[paste0("se_", endpoint)]]
    degraded[[paste0("p_", endpoint)]] <- 2 * pnorm(-abs(z))
  }
  rows[[i]] <- permutation_row(degraded, fraction)
  cat(sprintf("completed fraction %.2f (%d/%d)\n",
              fraction, i, length(fractions)))
}

result <- do.call(rbind, rows)
stopifnot(nrow(result) == 11L,
          identical(result$n_permutations, rep(10000L, 11L)),
          all(diff(result$observed_proteins) >= 0),
          all(diff(result$observed_loci) >= 0))
write.csv(result, output_file, row.names = FALSE)
cat("output:", output_file, "\n")
cat("elapsed_seconds:", proc.time()[["elapsed"]] - start_time, "\n")
