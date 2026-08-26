#' Function to wrap Pathway analysis pipeline
#' @export
#' @examples
#' ## Placeholder Example ##
run_enriched_paths<-function(proteins,analysis=c("all","PFG","EnrichR")
                             ,original_name_key="original",pfg_list=NULL
                             ,repository_name_key="repository"
                             ,libs='default',enrichr_filt_ptn=FALSE
                             ,min_ptn_enrichr_path=3,enrichr_filt_path=FALSE
                             ,min_path_enrichr_ptn=3,out_rm=TRUE,drop_paths=TRUE
                             ,path_number=20,max_str_width=40){

  analysis<-match.arg(analysis)

  if (!is.list(proteins))
    stop("proteins (input) must be a list of proteins.")

  # Extract keys only if needed
  if (analysis %in% c("all", "PFG")) {
    if (!original_name_key %in% names(proteins))
      stop("original_name_key not found in proteins list.")

    list1<-proteins[[original_name_key]]

    if (!is.list(list1))
      stop("original_name_key must refer to a list of protein vectors.")
  }

  if (analysis %in% c("all","EnrichR")) {
    if (!repository_name_key %in% names(proteins))
      stop("repository_name_key not found in proteins list.")

    list2<-proteins[[repository_name_key]]

    if (!is.list(list2))
      stop("repository_name_key must refer to a list of protein vectors.")
  }

  if (analysis %in% c("all","PFG")) {
    subset_names<-names(list1)
  } else {
    subset_names<-names(list2)
  }

  if (analysis=="all") {
    if (!identical(names(list1),names(list2)))
      stop("original and repository lists must have identical subset names.")
  }

  if (is.null(subset_names)) subset_names<-"all"  #fallback if the input has just one class

  #Inner function to process each subset
  process_subset<-function(original_proteins,repository_proteins) {

    results<-list()

    #If PFG analysis requested
    if (analysis %in% c("all","PFG")) {
      original_names<-unique(as.character(unlist(original_proteins))) #Fetch protein for PFG analysis

      #Safety checks
      if (length(original_names)==0) stop("No valid names for PFG")
      if (is.null(pfg_list)) stop("Reference PFG should be passed as pfg_list")

      pfg_res_df<-pfg_path(proteins=original_names,pfg_list=pfg_list) #function for PFG analysis
      results$pfg_tab<-pfg_res_df

      #Create plot for PFG and store other results
      results$pfg_plot<-plot_enrichment_paths(df=pfg_res_df
                                              ,path_number=path_number
                                              ,type='PFG')
      results$input_pfg<-original_names
      results$size_input_pfg<-length(original_names)
    }

    #If EnrichR analysis requested
    if (analysis %in% c("all","EnrichR")) {
      repository_names<-unique(as.character(unlist(repository_proteins))) #Fetch protein for EnrichR analysis
      if (length(repository_names)==0)
        stop("No valid names for EnrichR")

      enrichr_res_df<-enrichr_path(
        proteins=repository_names,
        libs=libs,
        enrichr_filt_ptn=enrichr_filt_ptn,
        min_ptn_enrichr_path=min_ptn_enrichr_path,
        enrichr_filt_path=enrichr_filt_path,
        min_path_enrichr_ptn=min_path_enrichr_ptn,
        out_rm=out_rm,
        drop_paths=drop_paths
      )

      #Store Results
      results$enrichr_tab<-enrichr_res_df
      results$input_enrichr<-repository_names
      results$size_input_enrichr<-length(repository_names)

      #Adjust EnrichR df for plotting and store results
      enrichr_res_df_adj<-enrichr_res_df
      enrichr_res_df_adj$Library<-gsub("_"," ",enrichr_res_df_adj$Library)
      enrichr_res_df_adj$Library<-gsub("[0-9]+", "", enrichr_res_df_adj$Library)
      enrichr_res_df_adj$Library<-gsub("Human", "", enrichr_res_df_adj$Library)
      enrichr_res_df_adj$Library<-gsub("Pathways", "", enrichr_res_df_adj$Library)
      enrichr_res_df_adj$Pathway<-gsub("\\s*\\([^\\)]+\\)","",enrichr_res_df_adj$Pathway)
      enrichr_res_df_adj$Pathway<-gsub("R-HSA-.*","",enrichr_res_df_adj$Pathway)
      enrichr_res_df_adj$Pathway<-gsub(" WP.*","",enrichr_res_df_adj$Pathway)
      enrichr_res_df_adj$Pathway<-gsub(" WP.*","",enrichr_res_df_adj$Pathway)
      enrichr_res_df_adj$Pathway<-stringr::str_trunc(enrichr_res_df_adj$Pathway
                                                     ,width=max_str_width,side="right")
      results$enrichr_plot<-plot_enrichment_paths(df=enrichr_res_df_adj
                                                  ,path_number=path_number
                                                  ,type='EnrichR')

    }

    #If all analysis were requested
    if (analysis=="all") {
      #Create final dataframes
      pfg_df_final<-results$pfg_tab |> select('Pathway','Library','# Proteins in Pathway','# Sig. Correlated Proteins'
                                               ,'% Correlated','pvalue','qvalue','Significance'
                                               ,'Odds_Ratio','Combined_Score','Correlated_Proteins')
      enrichr_df_final<-results$enrichr_tab
      combined_df<-rbind(pfg_df_final,enrichr_df_final)
      rownames(combined_df)<-NULL
      results$enriched_paths_tab<-combined_df
    }
    return(results)
  }

  #Loop through all subsets
  results<-lapply(subset_names,function(sub) {
    ori<- if (analysis %in% c("all","PFG")) { list1[[sub]] } else NULL
    repo<- if (analysis %in% c("all","EnrichR")) { list2[[sub]] } else NULL
    process_subset(ori,repo)
  })
  names(results)<-subset_names
  return(results)
}
