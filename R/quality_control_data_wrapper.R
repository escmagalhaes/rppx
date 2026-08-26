#' Function to wrap QC pipeline
#' @export
#' @examples
#' ## Placeholder Example ##
run_QC_analysis<-function(df,proteins,verbose=TRUE,grouping_variable=NULL
                          ,run_sample_qc=TRUE,sample_outlier_policy="keep_all"
                          ,sample_rm_thr=30,sample_diff_thr=20
                          ,run_feature_qc=TRUE,feature_outlier_policy="keep_all"
                          ,check_cor=TRUE,check_mad=FALSE,cor_cutoff=0.95
                          ,ptn_z_thrs=3,freq_thrs=0.10){

  #Match arguments
  sample_outlier_policy<-match.arg(sample_outlier_policy,
                                   c("keep_all","remove_strong","remove_all","interactive"))

  feature_outlier_policy<-match.arg(feature_outlier_policy,
                                    c("keep_all","remove","interactive"))

  #Ensure protein vector is correct
  proteins_filtered<-unique(intersect(proteins,colnames(df)))

  #Create empty objects in case both QC are FALSE
  df_filtered<-df
  removed_samples<-character(0)
  removed_features<-character(0)
  imbalance_flag<-FALSE
  sample_detect<-NULL
  feature_detect<-NULL

  ##Sample QC##
  if(run_sample_qc){

    if(verbose) cat("Running sample QC...\n")

    sample_detect<-detect_sample_outliers(df=df,proteins=proteins_filtered)

    #show plots
    print(sample_detect$plots$plot1_sample)
    print(sample_detect$plots$plot2_sample)

    sample_policy<-apply_sample_policy(
      df=df,strong=sample_detect$strong,
      moderate=sample_detect$moderate,
      policy=sample_outlier_policy,
      grouping_variable=grouping_variable,
      sample_rm_thr=sample_rm_thr,
      sample_diff_thr=sample_diff_thr,
      verbose=verbose)

    df_filtered<-sample_policy$df_filtered
    removed_samples<-sample_policy$removed
    imbalance_flag<-sample_policy$imbalance_flag
  }

  ##Feature QC##
  if(run_feature_qc){

    if(verbose) cat("Running feature QC...\n")
    proteins_feature<-intersect(proteins_filtered,colnames(df_filtered))
    feature_detect<-detect_feature_outliers(df=df_filtered,proteins=proteins_feature
                                            ,check_cor=check_cor,check_mad=check_mad
                                            ,cor_cutoff=cor_cutoff,ptn_z_thrs=ptn_z_thrs
                                            ,freq_thrs=freq_thrs,feature_outlier_policy=feature_outlier_policy)

    feature_policy<-apply_feature_policy(df=df_filtered,proteins=proteins_feature
                                         ,nzv=feature_detect$nzv,cor=feature_detect$cor
                                         ,mad=feature_detect$mad,policy=feature_outlier_policy
                                         ,freq_thrs=freq_thrs,verbose=verbose)

    df_final<-feature_policy$df
    ptn_final<-intersect(feature_policy$proteins, colnames(df_final))
    removed_features<-feature_policy$removed

  } else {
    df_final<-df_filtered
    ptn_final<-proteins_filtered
  }

  #Make a QC report
  qc_report<-build_qc_report(
    n_samples_original=nrow(df),
    n_samples_final=nrow(df_final),
    n_proteins_original=length(proteins_filtered),
    n_proteins_final=length(ptn_final),
    strong_outliers=if(!is.null(sample_detect)) sample_detect$strong else character(0),
    moderate_outliers=if(!is.null(sample_detect)) sample_detect$moderate else character(0),
    samples_removed=removed_samples,
    proteins_removed=removed_features,
    freq_threshold=freq_thrs,
    imbalance_flag=imbalance_flag,
    sample_plots=if(!is.null(sample_detect)) sample_detect$plots else NULL,
    feature_plots=if(!is.null(feature_detect)) feature_detect$raw else NULL
  )

  #Create a table with HGCN compatible protein names for downstream analysis
  clean_names<-clean_protein_names(protein_names=ptn_final,strict=TRUE)
  name_map_tab<-clean_names$name_mapping[
    ,c('original',"base_name","ptm",'ptm_type','mapping_level','final_names')]

  #Print partial QC report
  if(verbose) {
    cat(paste0(
      "\nQC completed successfully. You are all set!\n",
      "Samples analyzed: ", nrow(df), "\n",
      "Sample outliers removed: ",length(removed_samples),
      " (", if(length(removed_samples)==0) "none" else paste(removed_samples,collapse=", "), ")\n",
      "Proteins analyzed: ", length(proteins_filtered), "\n",
      "Protein outliers removed: ", length(removed_features),
      " (", if(length(removed_features)==0) "none" else paste(removed_features,collapse=", "), ")\n"
    ))
  }

  invisible(list(
    df_final=df_final,
    proteins_final=ptn_final,
    qc_report=qc_report,
    qc_results=list(removed_samples=removed_samples,removed_proteins=removed_features),
    diagnostics=list(sample_qc=sample_detect,feature_qc=feature_detect
                     ,strong_outliers=if(!is.null(sample_detect)) sample_detect$strong else character(0)
                     ,moderate_outliers=if(!is.null(sample_detect)) sample_detect$moderate else character(0)
                     ,imbalance_flag=imbalance_flag,freq_threshold=freq_thrs),
    ptn_name_mapping=name_map_tab
  ))
}
