# Required inputs and access limits

No upstream source data, tokens or account files are included. The verbatim dense-grid script remains in the private working archive. This candidate includes only a separately labelled path-adapted copy. Exclusion pending terms review is **not** a claim that redistribution is prohibited.
`source_inputs_manifest.json` lists checksums for the locally retained calibration inputs;
those checksums identify snapshots but do not establish their scientific provenance.

## Calibration scripts

Both scripts expect `inputs/01_extraction_progress.rds` (saved instrument records).
The sensitivity script also expects:

- `inputs/PosCtrl_BMI_all_results.csv`
- `inputs/19_mr_t2d_results.csv`
- `inputs/10_mr_hf_results.csv`
- `inputs/11_triangulation_combined.csv`
- `inputs/04_mr_results_corrected.csv`

The five-configuration script uses BMI, T2D and HF above, plus:

- `inputs/30_mr_ldl_results.csv`
- `inputs/31_mr_cad_results.csv`
- `inputs/35_mr_ckd_results.csv`
- `inputs/36_mr_sbp_results.csv`
- `inputs/37_mr_stroke_results.csv`
- `inputs/39_mr_masld_results.csv`
- `inputs/11_triangulation_combined.csv` for ATTR seed sensitivity.

An optional `source_baseline/23_permutation_stability.csv` produces a historical
comparison only. It is not required for the corrected calibration outputs.
Saved inputs are endpoint-level MR records, not raw participant-level GWAS data.
Exact original CTS/HF harmonisation and all primary-run versions remain incompletely
traced; these later scripts do not prove the earlier workflow.

## PON1 exploratory analysis

Place the exact saved snapshots (after authorised access) under `upstream_inputs/pon1/`:

- `olink_PON1_chr7.gz`: UKB-PPP European discovery Olink PON1 chromosome 7 file;
  upstream dataset Synapse `syn51364943`. This is a local alias, not an upstream filename.
- `finngen_R12_E4_AMYLNAS.gz`: FinnGen R12 E4_AMYLNAS summary data.
- `associations_direct.json`, `associations_direct.json.provenance.json`,
  `gwas_metadata.json`: saved OpenGWAS rs662 records and query provenance for
  `ebi-a-GCST006867` (T2D), `ebi-a-GCST009541` (HF), `ieu-b-40` (BMI).

The script validates alleles, IDs and the association-snapshot checksum. A fresh
API retrieval may differ and must not be silently substituted for the saved snapshot.
The deCODE PON1 file and GWAS Catalog `GCST90085767` support separate supplementary
checks; they are not inputs to this Olink Wald-ratio script and are not included here.

## ITIH1 diagnostic

Under `upstream_inputs/itih1/`, the script requires saved
`01_extraction_progress.rds`, `02_harmonised_data.rds`, `03_mr_results_raw.rds`
and `01_PWMR_pipeline.R`. The last file is hashed only, never sourced.
The branch is Sun `prot-a-1586` with FinnGen `finn-b-E4_AMYLOIDOSIS`, using
rs1042779, rs1056524 and rs4298806. Historical software versions cannot be inferred
from the saved RDS objects. No four-SNP dataset is inferred or supplied.

## Specific remaining gaps

1. The main 10,000-permutation dense-grid driver and execution record have been recovered, and its historical CSV matches the saved file byte for byte. Historical input hashes and full runtime information are not present in that record.
2. Historical input identity and complete runtime links remain limited to recovered records. R5/R12 endpoint metadata and the local rs1042779 source record have been checked; do not treat that as a replay of all primary analyses.
3. Author-approved MIT applies to author code; source-data compatibility and result-sharing scope remain separately unresolved.
4. Cite and verify the deposited Git commit for this directory; no complete data archive or DOI is implied.

The recovered original expects the saved instrument RDS plus BMI, T2D and HF result CSVs listed above. Its statistical-code suffix is preserved verbatim in the path-adapted copy, which has not been run. Finding the driver does not establish historical input identity or close the other gaps.
