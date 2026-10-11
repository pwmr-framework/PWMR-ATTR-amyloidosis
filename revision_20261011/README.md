# JHG-26-335: revision analysis code

Code-only revision snapshot prepared 11 October 2026 for https://github.com/pwmr-framework/PWMR-ATTR-amyloidosis. This directory is an additive update; earlier repository files and history are retained. Cite the exact Git commit containing this directory. No DOI is assigned.

## Scope and licence

Five author analysis scripts are supplied with the author-approved MIT licence (Copyright 2026 Jincheng Xing). LICENSE covers author-owned code and original accompanying documentation, not third-party data, source identifiers, citations, dependencies, or other people's code. Software dependencies retain their own licences.

This directory contains **no saved result datasets and no upstream data**. The author approved publication of this code-only scope. Existing manuscript tables and privately retained results have not been removed or replaced.

## Files and dependencies

- analysis/dense_power_10000_path_adapted.R: path-adapted recovered 10,000-permutation dense-grid driver; base R.
- analysis/reviewer_sensitivity_analysis.R: threshold/anchor/order sensitivity; its power branch uses 1,000 permutations and must not replace the dense-grid driver; base R.
- analysis/unified_trajectory_calibration.R: additional-configuration and ATTR seed calibration; base R.
- analysis/run_pon1.R: exploratory Olink PON1 single-variant analysis; TwoSampleMR, digest and jsonlite.
- analysis/recompute_saved_sun_itih1_diagnostic.R: saved Sun three-SNP ITIH1 diagnostic; TwoSampleMR and /usr/bin/shasum.
- INPUTS_AND_LIMITATIONS.md: required saved snapshots and missing provenance.
- source_inputs_manifest.json: twelve calibration-input hashes, not a complete PON1/ITIH1 input manifest; inputs are not included.
- manifest.json: included-file hashes and lengths. It excludes its own hash.

## Running and limitations

Read INPUTS_AND_LIMITATIONS.md first and obtain permitted inputs from their providers. Source GWAS records are not interchangeable with the saved harmonised or MR intermediate tables expected by these scripts. This is a **partial code archive**, not a complete reconstruction of the primary MR pipeline or a certified one-command reproduction package.

Run a selected script from this directory only after its exact required inputs and dependencies are available, for example:

    Rscript analysis/dense_power_10000_path_adapted.R

Do not run the calibration scripts concurrently. Outputs go to rerun_outputs/, which is excluded from version control; some outputs contain upstream rows and must not be published automatically. The dense-grid driver refuses to overwrite its existing output; other scripts may overwrite their corresponding rerun outputs.

The scripts were copied byte-for-byte from the reviewed local candidate and were not executed in preparing this package. Parse-only checks do not verify results. Earlier diagnostic session records used R 4.5.3 and TwoSampleMR 0.7.4; these do not establish every historical primary-run environment.

## Evidence boundaries

PON1 is exploratory single-variant analysis, not certified independent disease-outcome replication. ITIH1 is the saved three-SNP diagnostic, not a reinstated four-SNP validation. Locus grouping is the implemented analysis rule, not new LD verification. The source record for the historical dense grid was recovered, but full historical input/runtime identity remains incompletely documented. The package does not include complete source data for every manuscript figure.

## Source attribution

Source conditions and citations remain applicable: OpenGWAS https://opengwas.io/policies/eula ; UKB-PPP Olink syn51364943 https://doi.org/10.1038/s41586-023-06592-6 ; FinnGen https://www.finngen.fi/en/access_results . The supplementary manuscript Methods and Table 1 contain the study-specific source descriptions.

For OpenGWAS access, cite the infrastructure (doi:10.1101/2020.08.10.244293) and MR-Base (doi:10.7554/eLife.34408), as well as the original dataset publications. Source access or a source citation does not transfer its data licence to this code package.
