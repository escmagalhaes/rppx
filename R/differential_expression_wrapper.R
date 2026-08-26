#' Function to wrap DE analysis pipeline
#' @export
#' @examples
#' ## Placeholder Example ##
run_DE_analysis<-function(df,proteins,grouping_variable=NULL,verbose=TRUE,
                          run_sample_qc=TRUE,sample_outlier_policy="keep_all",
                          sample_rm_thr=30,sample_diff_thr=20,
                          run_feature_qc=TRUE,feature_outlier_policy="keep_all",
                          check_cor=TRUE,check_mad=FALSE,cor_cutoff=0.95,
                          ptn_z_thrs=3,freq_thrs=0.10,de_mode=c("binary","ordered"),
                          pval_cutoff=0.05,LFC_cutoff=0.25,ref_level=NULL,strict=FALSE,
                          treat=FALSE,covars=NULL,vplot_legend_position=c(0.8,0.9)){

  de_mode<-match.arg(de_mode)

  #QC analysis
  qc_results<-run_QC_analysis(
    df                     = df,
    proteins               = proteins,
    grouping_variable      = grouping_variable,
    verbose                = verbose,

    run_sample_qc          = run_sample_qc,
    sample_outlier_policy  = sample_outlier_policy,
    sample_rm_thr          = sample_rm_thr,
    sample_diff_thr        = sample_diff_thr,

    run_feature_qc         = run_feature_qc,
    feature_outlier_policy = feature_outlier_policy,
    check_cor              = check_cor,
    check_mad              = check_mad,
    cor_cutoff             = cor_cutoff,
    ptn_z_thrs             = ptn_z_thrs,
    freq_thrs              = freq_thrs
  )

  #Run DE model
  if(verbose) cat("Running differential expression analysis model...\n")
  df_qc<-qc_results$df_final
  proteins_qc<-qc_results$proteins_final

  if (de_mode=='binary'){
    #Function to run binary DE limma model
    mod<-binary_de_model(
      df           = df_qc,
      grouping_var = grouping_variable,
      proteins     = proteins_qc,
      covars       = covars,
      LFC_cutoff   = LFC_cutoff,
      pval_cutoff  = pval_cutoff,
      treat        = treat,
      ref_level    = ref_level,
      verbose      = verbose
    )
  } else {
    #Function to run ordered DE limma model with more than 2 levels
    mod<-ordered_de_model(
      df           = df_qc,
      grouping_var = grouping_variable,
      proteins     = proteins_qc,
      covars       = covars,
      LFC_cutoff   = LFC_cutoff,
      pval_cutoff  = pval_cutoff,
      verbose      = verbose
    )
  }

  res_data<-mod$results
  res_data<-tibble::rownames_to_column(res_data,'original')

  #Merge DE results with name mapping table for downstream analysis
  name_map_tab_final<-res_data |>
    dplyr::left_join(qc_results$ptn_name_mapping,by="original") |>
    dplyr::mutate(map_id=row_number()) |>
    dplyr::select(map_id,everything())
  name_map_tab_final<-name_map_tab_final[order(name_map_tab_final$final_names),]

  #Extract cleaned names of DE proteins
  #Create table subset for protein networks and extract cleaned names
  name_map_tab_final_subset<-name_map_tab_final[name_map_tab_final$is_signif==TRUE
                                                & !is.na(name_map_tab_final$final_names),]
  name_map_tab_final_subset<-name_map_tab_final_subset[order(name_map_tab_final_subset$final_names),]
  cleaned_sig_names<-unique(na.omit(name_map_tab_final_subset[,'final_names']))

  direction_col<- if (de_mode=="binary") "class_high" else "direction"
  cleaned_up_names<-unique(na.omit(
    name_map_tab_final_subset[name_map_tab_final_subset[[direction_col]]=='up','final_names']))
  cleaned_down_names<-unique(na.omit(
    name_map_tab_final_subset[name_map_tab_final_subset[[direction_col]]=='down','final_names']))

  return(list(
    expr_data=df_qc,
    ptn=proteins_qc,
    qc_results=qc_results,
    summary=list(
      fit=mod$fit,
      vplot=mod$vplot,
      ptn=list(original_names=mod$all_sig_ptn
               ,repository_names=cleaned_sig_names),
      ptn_class=list(original=mod$sig_ptn_class
                     ,repository=list(all=cleaned_sig_names
                                      ,up=cleaned_up_names
                                      ,down=cleaned_down_names)),
      ptn_size=mod$sig_ptn_size,
      mapping_table_full=name_map_tab_final,
      mapping_table_subset=name_map_tab_final_subset
    )
  ))

}
