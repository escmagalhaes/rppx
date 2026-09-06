Changelog

All changes to rppx will be documented below

## [Unreleased]

### Added
- plot_heatmap(): new optional col_id argument to specify the column with sample identifiers. Used when heatmap colnames need to show sample identifiers.

- run_QC_analysis(): new logical strict argument (default=TRUE). This is used when calling clean_protein_names() in this QC wrapper. Previously, strict=TRUE was hardcoded, making run_DE_analysis() to fail.

- run_DE_analysis(): new logical strict argument (default=TRUE) added to this wrapper. his is used when calling clean_protein_names() in this QC wrapper. Previously, strict=TRUE was hardcoded, making run_DE_analysis() to fail.

### Changed
- plot_heatmap(): heatmap colnames now use sample identifiers (from col_id, or rownames(df) if col_id is NULL) instead of being coerced to sequential numbers

- run_QC_analysis(): now strict argument (default=TRUE) allows for function to run without stop if there is an issue with feature names (non-HGCN approved). Use strict=FALSE to prevent function to stop. 

### Fixed

- plot_heatmap(): column names in col_ann and ht_mtx were overwritten with seq_len(). This caused heatmaps to display numeric indexes instead of sample identifiers.

- run_DE_analysis(): new strict = strict argument added to qc_results<-run_QC_analysis(). Now  strict argument can be controlled inside run_QC_analysis() from the run_DE_analysis() wrapper.  

[0.1.3-beta] - 2026-08-26

Added

Initial public beta release.
