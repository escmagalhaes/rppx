#' Function to wrap correlations pipeline for multiple datasets
#' @export
#' @examples
#' ## Placeholder Example ##
run_corr_pipeline<-function(df_input,proteins_of_interest,all_proteins,
                            grouping_variable=NULL,var_labs=NULL,df_normal=NULL,
                            df_subsetting=FALSE,min_samples=10,ncol_grid=2,verbose=TRUE,
                            method='spearman',conf_level=0.95,
                            pval_cutoff=0.05,abs_cutoff=0.20,
                            run_sample_qc=TRUE,sample_outlier_policy="keep_all",
                            sample_rm_thr=30,sample_diff_thr=20,
                            run_feature_qc=TRUE,feature_outlier_policy="keep_all",
                            check_cor=TRUE,check_mad=FALSE,cor_cutoff=0.95,
                            ptn_z_thrs=3,freq_thrs=0.10,p_adjust_method="BH",
                            single_dataset_label="all") {

  #Helper for blank WF plot
  blank_wf_plot<-function(df_name,protein,n) {
    ggplot()+geom_blank()+
      annotate("text",x=0.5,y=0.5,label=paste0("Analysis not done\nDataset: ",df_name,
                                               "\nProtein: ",protein, "\nN = ",n)
               ,size=5,hjust=0.5,vjust=0.5)+theme_void()
  }

  #Normalize grouping_variable
  if (!is.null(grouping_variable)) {
    if (is.factor(grouping_variable)) grouping_variable<-as.character(grouping_variable)
    if (!isTRUE(df_subsetting) && !is.null(grouping_variable) && length(grouping_variable) > 1) {
      stop("When df_subsetting=FALSE, grouping_variable must be NULL or length=1")
    }
  }

  #Ensure proteins are character
  proteins_of_interest<-as.character(proteins_of_interest)

  #QC analysis
  qc_results<-run_QC_analysis(
    df=df_input,
    grouping_variable=grouping_variable,
    proteins=all_proteins,
    verbose=verbose,
    run_sample_qc=run_sample_qc,
    sample_outlier_policy=sample_outlier_policy,
    sample_rm_thr=sample_rm_thr,
    sample_diff_thr=sample_diff_thr,
    run_feature_qc=run_feature_qc,
    feature_outlier_policy=feature_outlier_policy,
    check_cor=check_cor,
    check_mad=check_mad,
    cor_cutoff=cor_cutoff,
    ptn_z_thrs=ptn_z_thrs,
    freq_thrs=freq_thrs
  )
  df_qc<-qc_results$df_final
  proteins_qc<-qc_results$proteins_final

  #Create subsets if requested
  df_list<-list()
  subset_summary<-data.frame()
  if (isTRUE(df_subsetting)) {
    if (is.null(proteins_of_interest) && is.null(grouping_variable)) {
      stop("Subsetting requires at least one of protein(s) or grouping_variable.")
    }
    subset_results<-create_all_subsets(
      df_normal=df_normal,
      df_disease=df_input,
      df_apply=df_qc,
      protein_list=proteins_of_interest,
      variable_list=grouping_variable,
      var_labs=var_labs,
      print_summary=verbose
    )
    #Flatten nested list to df_list and create names
    df_list<-list()
    subset_summary<-data.frame()
    subset_counter<-1
    for(subset_key in names(subset_results)) {
      for(subname in names(subset_results[[subset_key]]$subsets)) {
        dataset_id<-paste0("dataset", sprintf("%03d", subset_counter))
        dataset_label<-subname
        df_list[[dataset_id]]<-subset_results[[subset_key]]$subsets[[subname]]
        subset_info<-subset_results[[subset_key]]$summary[
          subset_results[[subset_key]]$summary$subset_name==subname,]
        subset_info$dataset<-dataset_id
        subset_info$dataset_label<-dataset_label
        subset_summary<-rbind(subset_summary,subset_info)
        subset_counter<-subset_counter+1
      }
    }
    if(length(df_list)==0){
      stop("No subsets created. Check protein thresholds or grouping variables.")
    }
    if(verbose){
      cat("\nTotal subsets created:",length(df_list),"\n")
    }
  } else {
    df_list<-list(dataset1=df_qc)
    if(verbose) cat("Single dataframe detected, wrapping into list 'dataset1'\n")
  }

  #Initialize results
  results_list<-vector("list",length(df_list))
  names(results_list)<-names(df_list)

  summary_table<-data.frame(
    dataset_name=names(df_list),
    n_samples=NA_integer_,
    analyzed=FALSE,
    n_signif_correlations=NA_integer_,
    note="",
    row.names=NULL
  )

  #Loop over datasets / subsets
  for (dataset_idx in seq_along(df_list)) {
    df_name<-names(df_list)[dataset_idx]
    df_sub<-df_list[[dataset_idx]]
    n<-nrow(df_sub)
    summary_table$n_samples[dataset_idx]<-n

    if (n < min_samples) {
      #Skip analysis, create blank WF plots
      wf_plots<-lapply(proteins_of_interest,function(prot) blank_wf_plot(df_name,prot,n))
      names(wf_plots)<-proteins_of_interest
      wf_grid<-ggpubr::ggarrange(plotlist=wf_plots,ncol=min(ncol_grid,length(wf_plots)),
                                 nrow=ceiling(length(wf_plots) / ncol_grid))
      wf_grid<-ggpubr::annotate_figure(
        wf_grid,top=ggpubr::text_grob(paste0(df_name," (N=",n, ")"),face="bold",size=12))

      results_list[[dataset_idx]]<-list(
        df_filtered=df_sub,
        all_proteins_filtered=all_proteins,
        qc_results=NULL,
        summary=list(
          all_correlations=NULL,
          ptn_class=NULL,
          ptn_size=0,
          mapping_table_full=NULL,
          mapping_table_subset=NULL,
          wfplots=wf_plots,
          wf_plot_grid=wf_grid
        )
      )
      summary_table$analyzed[dataset_idx]<-FALSE
      summary_table$n_signif_correlations[dataset_idx]<-0
      summary_table$note[dataset_idx]<-"Skipped: < min_samples"
      next
    }

    #Run correlation analysis
    if (verbose) cat("Running correlation analysis for dataset:",df_name,"\n")
    cor_res<-ptn_corr(df=df_sub,
                      all_proteins=proteins_qc,
                      proteins_of_interest=proteins_of_interest,
                      method=method,
                      conf_level=conf_level,
                      pval_cutoff=pval_cutoff,
                      abs_cutoff=abs_cutoff,
                      p_adjust_method=p_adjust_method)

    #Add Mean expression value per protein
    mean_expr<-df_sub |>
      dplyr::select(dplyr::all_of(proteins_qc)) |>
      colMeans(na.rm=TRUE) |>
      tibble::enframe(name="protein",value="mean_expr")
    cor_res<-cor_res |>
      dplyr::left_join(mean_expr,by="protein") |>
      dplyr::rename(mean_expr_protein=mean_expr) |>
      dplyr::left_join(mean_expr,by=c("query_protein"="protein")) |>
      dplyr::rename(mean_expr_query_protein=mean_expr) |>
      dplyr::select(query_protein,protein,mean_expr_query_protein,mean_expr_protein,everything())

    #Export dataframe with cleaned protein names for downstream analysis
    ptn_name_mapping<-qc_results$ptn_name_mapping
    name_map_tab_final<-cor_res |>
      dplyr::left_join(ptn_name_mapping,by=c("protein"="original"),relationship="many-to-many")|>
      mutate(map_id=row_number()) |> select(map_id,everything()) |> arrange(final_names)

    #Significant correlations and names
    sig<-name_map_tab_final[name_map_tab_final$is_signif
                            & !is.na(name_map_tab_final$final_names), ]
    original_sig_names<-unique(sig$protein)
    original_split_names<-split(sig$protein,sig$sign)
    original_pos_names<-unique(original_split_names$pos)
    original_neg_names<-unique(original_split_names$neg)
    cleaned_sig_names<-unique(sig$final_names)
    cleaned_split_names<-split(sig$final_names,sig$sign)
    cleaned_pos_names<-unique(cleaned_split_names$pos)
    cleaned_neg_names<-unique(cleaned_split_names$neg)

    #Waterfall plots
    wf_plot_list<-list()
    cor_res_by_query<-split(cor_res,cor_res$query_protein)
    for(prot in proteins_of_interest) {
      cor_res_sub<-cor_res_by_query[[prot]]
      if(is.null(cor_res_sub) || nrow(cor_res_sub)==0) next
      neg_plot<-wf_plot(cor_res_sub,type='Negative')
      pos_plot<-wf_plot(cor_res_sub,type='Positive')
      wf_combo<-ggpubr::ggarrange(neg_plot, pos_plot, ncol=2, labels=NULL)
      wf_plot_list[[prot]]<-ggpubr::annotate_figure(
        wf_combo,top=ggpubr::text_grob(paste0("Significant correlations "
                                              ,prot," (N=",n," patients)")
                                       ,face="bold", size=8))
    }
    if(length(wf_plot_list)==0) wf_plot_list<-lapply(proteins_of_interest,function(prot) blank_wf_plot(df_name,prot,n))
    wf_grid<-ggpubr::ggarrange(plotlist=wf_plot_list,ncol=min(ncol_grid,length(wf_plot_list)),
                               nrow=ceiling(length(wf_plot_list)/ncol_grid))

    #Compile results
    results_list[[dataset_idx]]<-list(
      summary=list(
        all_correlations=cor_res,
        ptn_class=list(
          original=list(
            all=original_sig_names,
            pos=original_pos_names,
            neg=original_neg_names),
          repository=list(
            all=cleaned_sig_names,
            pos=cleaned_pos_names,
            neg=cleaned_neg_names)),
        ptn_size=length(original_sig_names),
        mapping_table_full=name_map_tab_final,
        mapping_table_subset=sig,
        wfplots=wf_plot_list,
        wf_plot_grid=wf_grid
      )
    )
    summary_table$analyzed[dataset_idx]<-TRUE
    summary_table$n_signif_correlations[dataset_idx]<-nrow(sig)
    summary_table$note[dataset_idx]<-"Analyzed"
  }

  #Combine correlations across datasets
  combined_correlations<-dplyr::bind_rows(lapply(seq_along(results_list),function(i) {
    cor_tab<-results_list[[ i ]]$summary$all_correlations
    if (is.null(cor_tab)) return(NULL)
    cor_tab$dataset<-names(results_list)[ i ]
    cor_tab$pair_id<-paste(cor_tab$query_protein,cor_tab$protein,sep="__")
    cor_tab
  }))
  if (!is.null(combined_correlations) && nrow(subset_summary) > 0) {
    combined_correlations<-dplyr::left_join(combined_correlations,
                                            subset_summary[,c("dataset","dataset_label")]
                                            ,by="dataset")
  } else {
    combined_correlations$dataset_label<-single_dataset_label
  }
  combined_correlations<-dplyr::left_join(combined_correlations,
                                          summary_table[,c("dataset_name","n_samples")]
                                          ,by=c("dataset"="dataset_name"))
  combined_correlations<-combined_correlations |>
    dplyr::select(dataset,query_protein,dataset_label,n_samples,everything())

  #Create an object with combined correlations with ptn names cleaned and mapped
  combined_correlations_mapped<-combined_correlations |>
    dplyr::left_join(ptn_name_mapping,by=c("protein"="original"),relationship="many-to-many")

  return(list(
    results=results_list,
    summary_table=summary_table,
    combined_correlations=combined_correlations,
    combined_correlations_mapped=combined_correlations_mapped,
    ptn_name_mapping=ptn_name_mapping,
    subset_summary=subset_summary  #Optional: only if df_subsetting=TRUE
  ))
}
