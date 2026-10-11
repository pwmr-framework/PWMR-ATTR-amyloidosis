# Diagnostic reproduction ONLY: saved Sun prot-a-1586, amyloidosis, three SNPs.
# This is NOT the old four-cis-SNP analysis and MUST NOT replace that table.
# No network, no harmonisation rerun, no new instrument selection, no sourcing.

options(stringsAsFactors = FALSE, digits = 17)
outdir <- file.path(normalizePath(getwd()), "rerun_outputs", "itih1")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
sources <- c(
  exposure = file.path(getwd(), "upstream_inputs", "itih1", "01_extraction_progress.rds"),
  harmonised = file.path(getwd(), "upstream_inputs", "itih1", "02_harmonised_data.rds"),
  saved_mr = file.path(getwd(), "upstream_inputs", "itih1", "03_mr_results_raw.rds"),
  original_pipeline_read_only = file.path(getwd(), "upstream_inputs", "itih1", "01_PWMR_pipeline.R"),
  diagnostic_script = file.path(getwd(), "analysis", "recompute_saved_sun_itih1_diagnostic.R")
)
stopifnot(all(file.exists(sources)))
sha256 <- function(path) {
  ans <- system2("/usr/bin/shasum", c("-a", "256", shQuote(path)), stdout = TRUE)
  stopifnot(length(ans) == 1L, is.null(attr(ans, "status")))
  sub(" .*", "", ans)
}
before_hash <- vapply(sources, sha256, character(1))
manifest <- data.frame(role = names(sources), path = unname(sources),
                       sha256 = unname(before_hash),
                       bytes = unname(file.info(sources)$size))
write.csv(manifest, file.path(outdir, "input_source_manifest.csv"), row.names = FALSE)

exposure_id <- "prot-a-1586"
outcome_id <- "finn-b-E4_AMYLOIDOSIS"
expected_snps <- sort(c("rs1042779", "rs1056524", "rs4298806"))
extract <- readRDS(sources[["exposure"]])
dat <- readRDS(sources[["harmonised"]])
saved <- readRDS(sources[["saved_mr"]])
e <- extract[extract$id.exposure == exposure_id, , drop = FALSE]
h <- dat[dat$id.exposure == exposure_id & dat$id.outcome == outcome_id, , drop = FALSE]
r <- saved[saved$id.exposure == exposure_id & saved$id.outcome == outcome_id &
             saved$method == "Inverse variance weighted", , drop = FALSE]
e <- e[order(e$SNP), , drop = FALSE]
h <- h[order(h$SNP), , drop = FALSE]

checks <- character()
check <- function(label, condition) {
  stopifnot(length(condition) == 1L, !is.na(condition), condition)
  checks <<- c(checks, label)
}
check("exact exposure ID and three expected extracted SNPs",
      nrow(e) == 3L && identical(e$SNP, expected_snps))
check("exact outcome ID and three expected harmonised SNPs",
      nrow(h) == 3L && identical(h$SNP, expected_snps))
check("no duplicate SNPs", !anyDuplicated(h$SNP))
check("all harmonised SNPs eligible", all(h$mr_keep & h$mr_keep.exposure & h$mr_keep.outcome))
check("no remove, ambiguous or palindromic flags", all(!h$remove & !h$ambiguous & !h$palindromic))
check("same effect and other alleles in aligned exposure/outcome",
      all(h$effect_allele.exposure == h$effect_allele.outcome &
            h$other_allele.exposure == h$other_allele.outcome))
check("harmonisation action recorded as 2", all(h$action == 2L))
check("no rs2286795 in this branch", !("rs2286795" %in% c(e$SNP, h$SNP)))
num_fields <- c("beta.exposure", "se.exposure", "beta.outcome", "se.outcome")
check("finite effect and standard error inputs", all(is.finite(as.matrix(h[, num_fields]))))
check("strictly positive standard errors", all(h$se.exposure > 0 & h$se.outcome > 0))
check("extraction and harmonised exposure records agree",
      identical(e$SNP, h$SNP) && all(e$beta.exposure == h$beta.exposure) &&
        all(e$se.exposure == h$se.exposure) &&
        all(e$effect_allele.exposure == h$effect_allele.exposure) &&
        all(e$other_allele.exposure == h$other_allele.exposure))
check("saved primary MR record is unique and has three SNPs", nrow(r) == 1L && r$nsnp == 3L)
write.csv(e, file.path(outdir, "saved_sun_itih1_extracted_inputs.csv"), row.names = FALSE)
write.csv(h, file.path(outdir, "saved_sun_itih1_amyloidosis_harmonised_inputs.csv"), row.names = FALSE)

# Reproduce the existing mr_ivw method only, without invoking data-access code.
calc <- TwoSampleMR::mr_ivw(h$beta.exposure, h$beta.outcome,
                          h$se.exposure, h$se.outcome)
# Independent explicit IVW calculation, with multiplicative random effects and
# the same underdispersion correction as TwoSampleMR::mr_ivw 0.7.4.
bx <- h$beta.exposure
by <- h$beta.outcome
w <- 1 / h$se.outcome^2
b_manual <- sum(w * bx * by) / sum(w * bx^2)
q_manual <- sum(w * (by - b_manual * bx)^2)
se_manual <- sqrt(max(1, q_manual / (length(bx) - 1)) / sum(w * bx^2))
manual <- c(b = b_manual, se = se_manual,
            pval = 2 * pnorm(abs(b_manual / se_manual), lower.tail = FALSE))
metrics <- c("b", "se", "pval")
comparison <- data.frame(
  metric = metrics,
  saved_value = unname(unlist(r[1, metrics])),
  recomputed_TwoSampleMR = unname(unlist(calc[metrics])),
  independent_manual = unname(manual[metrics])
)
comparison$absolute_difference_vs_saved <- abs(comparison$recomputed_TwoSampleMR - comparison$saved_value)
comparison$absolute_difference_vs_manual <- abs(comparison$recomputed_TwoSampleMR - comparison$independent_manual)
check("primary beta, SE and P match saved MR within absolute 1e-12",
      all(comparison$absolute_difference_vs_saved < 1e-12))
check("package and independent manual beta, SE and P agree within absolute 1e-12",
      all(comparison$absolute_difference_vs_manual < 1e-12))
write.csv(comparison, file.path(outdir, "diagnostic_recomputation_comparison.csv"), row.names = FALSE)
result <- data.frame(
  analysis_status = "DIAGNOSTIC_ONLY_NOT_FOUR_CIS_REPRODUCTION",
  exposure_id = exposure_id, outcome_id = outcome_id,
  method = "Inverse variance weighted", nsnp = calc$nsnp,
  snps = paste(h$SNP, collapse = ";"),
  beta = calc$b, se = calc$se, p = calc$pval,
  diagnostic_Q = calc$Q, diagnostic_Q_df = calc$Q_df,
  diagnostic_Q_p = calc$Q_pval
)
write.csv(result, file.path(outdir, "diagnostic_saved_sun_itih1_result.csv"), row.names = FALSE)
check("input file SHA-256 hashes unchanged after diagnostic computation",
      identical(before_hash, vapply(sources, sha256, character(1))))
write.csv(data.frame(check = checks, status = "PASS"),
          file.path(outdir, "diagnostic_assertions.csv"), row.names = FALSE)
writeLines(c(
  "Diagnostic reproduction only. This does not reproduce or replace the old four-cis-SNP analysis.",
  paste("Run UTC:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  paste("Current R:", R.version.string),
  paste("Current TwoSampleMR:", as.character(packageVersion("TwoSampleMR"))),
  "Original calculation runtime version is not established from the saved RDS alone.",
  "No API calls, new instruments, re-harmonisation, or alternative MR methods were used.",
  "The Q output is recorded only as an IVW implementation diagnostic, not a new sensitivity claim.",
  "",
  capture.output(sessionInfo())
), file.path(outdir, "diagnostic_runtime.txt"))
print(result)
print(comparison)
cat("PASS:", length(checks), "assertions. Old four-cis-SNP branch remains BLOCKED_INPUT_UNRESOLVED.\n")
