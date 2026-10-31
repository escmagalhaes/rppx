NEWS/Changelog

All changes to rppx will be documented below

## [Unreleased]

### Changed
- get_mycolors(): changed function name to get_colors() and added roxygen documentation

- clean_protein_names(): added roxygen documentation

- get_rawtable_values(): added roxygen documentation

- safe_clean_raw_data(): added roxygen documentation

- get_feature_dataset(): added roxygen documentation

- get_feature_dataset(): added roxygen documentation

- get_factor_levels(): added roxygen documentation and adjusted arguments

- remove_extra_terms(): added roxygen documentation

- adjust_datatable(): added roxygen documentation and added `date_format` argument

- determine_optimal_quantile(): added roxygen documentation and p-values are now used as second level sorting parameter to determine the best quantile split in `dplyr::arrange(q_value,p_value,n_quantiles)` and performed other small fixes

- get_quantile_groups(): added roxygen documentation and performed small fixes

- collapse_quantile_groups(): added roxygen documentation and performed small fixes

- get_var_groups(): added roxygen documentation and altered behavior to fit into get_subsets() and create_all_subsets() pipeline and performed small fixes

- get_subsets(): added roxygen documentation and altered behavior to fit into get_var_groups() and create_all_subsets() pipeline and performed small fixes

- create_all_subsets(): added roxygen documentation and altered behavior to fit into get_var_groups() and get_subsets() pipeline and performed small fixes

### Fixed
- manuscript_export(): small adjustments

- render_manuscript(): small adjustments

- filter_lookup_table(): small adjustments

- make_float_chunks(): small adjustments

- plot_density(): adjusted get_mycolors() calls to get_colors()

- get_stringdb_data(): adjusted arguments

- create_base_network(): adjusted arguments

- apply_base_visual_style(): small adjustments

- highlight_edge_sets(): small adjustments

- apply_visual_contrast(): small adjustments

- get_single_ptn_network(): adjusted arguments

- compare_ref_network(): adjusted arguments

- build_string_network(): adjusted arguments

- set_contrast(): small adjustments

- get_cbio_data(): adjusted arguments

- plot_heatmap(): adjusted get_mycolors() calls to get_colors() and adjusted `colors_anno_col` and `colors_anno_row` arguments by calling ggsci::pal_jco()

- create_bxplot(): added package calling ggpubr::stat_compare_means()

- plot_delta_volcano(): added package calling ggrepel::geom_text_repel()

- apply_visual_contrast(): corrected argument call of setEdgeLineWidthBypass() from `new.width` to `new.widths`

- plot_bars(): adjusted get_mycolors() calls to get_colors()

- run_corr_pipeline(): adapted function to accept changes in create_all_subsets(), added an option to use pre-defined quantile groups (e.g. High vs Low), using `discretize_proteins = FALSE`

### Removed
- get_ptn_groups(): removed completely in accordance to create_all_subsets() pipeline change


## [0.1.4-beta] - 2026-10-18

### Added
- plot_heatmap(): new optional col_id argument to specify the column with sample identifiers. Used when heatmap colnames need to show sample identifiers.

- plot_heatmap(): safety check to guarantee that grouping_var is a factor with pre-defined levels

- plot_heatmap(): added optional argument row_anno_colors_override to specify the colors used for
row annotation.

- run_QC_analysis(): new logical strict argument (default=TRUE). This is used when calling clean_protein_names() in this QC wrapper. Previously, strict=TRUE was hardcoded, making run_DE_analysis() to fail.

- run_DE_analysis(): new logical strict argument (default=TRUE) added to this wrapper. his is used when calling clean_protein_names() in this QC wrapper. Previously, strict=TRUE was hardcoded, making run_DE_analysis() to fail.

- plot_bars(): function to create barplots with ggplot2 according to package style

- plot_corr(): function to create correlation plots with ggplot2 according to package style

### Changed
- plot_heatmap(): heatmap colnames now use sample identifiers (from col_id, or rownames(df) if col_id is NULL) instead of being coerced to sequential numbers

- plot_heatmap(): added a fallback to colors_anno_row which passes more colors if the user pass
less colors than the levels of the row annotation variable

- run_QC_analysis(): now strict argument (default=TRUE) allows for function to run without stop if there is an issue with feature names (non-HGCN approved). Use strict=FALSE to prevent function to stop. 

- get_mycolors(): changed color order in "vector" mode ('darkturquoise' is now at the 6th position and 'darkcyan' at 10th)

### Fixed

- plot_heatmap(): column names in col_ann and ht_mtx were overwritten with seq_len(). This caused heatmaps to display numeric indexes instead of sample identifiers.

- plot_heatmap(): adjusted grouping_var wrangling to not drop unused levels if one is empty

- plot_heatmap(): removed sort() calls when adjusting the order of rows according to row annotation.
Now, the user can pass the desired row order and it will not be coerced to alphabetical anymore.

- run_DE_analysis(): new strict = strict argument added to qc_results<-run_QC_analysis(). Now  strict argument can be controlled inside run_QC_analysis() from the run_DE_analysis() wrapper.  

- determine_optimal_quantile(): q- and p-values being rounded only in the final results table

- plot_km(): added "pkg::" before function calls and adjusted theme, including scales, adjusted order of colors in "colors" argument

- plot_density(): adjusted ggplot theme

## [0.1.3-beta] - 2026-08-26

### Added
Initial public beta release.
