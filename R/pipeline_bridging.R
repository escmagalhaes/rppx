### Bridge data between pipelines ###

#' Function to pre-process Correlations output for pathways and protein networks input
#' @export
#' @examples
#' ## Placeholder Example ##
preprocess_corr_sig<-function(cor_df,datasets,df_ptn_map,query_ptn=NULL) {

  cor_df<-cor_df |> dplyr::filter(dataset_label %in% datasets)

  #Optional filter by query protein
  if (!is.null(query_ptn)) cor_df<-cor_df |> dplyr::filter(query_protein %in% query_ptn)

  #Convert mean_expr_query to meaningful row
  query_rows<-cor_df |> dplyr::group_by(query_protein,dataset_label) |>
    dplyr::slice_head(n=1) |>       #one row per query protein, all columns intact
    dplyr::ungroup() |>
    dplyr::mutate(
      dplyr::across(!c(dataset,query_protein,dataset_label,n_samples
                       ,mean_expr_query_protein),~NA),  #set all non-meaningful to NA
      protein=query_protein,corr=1,abs_corr=1,is_signif=TRUE,sign="pos"
      ,mean_expr_protein=mean_expr_query_protein)
  cor_df<-dplyr::bind_rows(cor_df,query_rows)

  #Pivot table
  cor_df<-cor_df |> tidyr::pivot_wider(id_cols=c(query_protein,protein)
                                        ,names_from=dataset_label
                                        ,values_from=c(is_signif,abs_corr,sign,mean_expr_protein)) |>
    dplyr::filter(dplyr::if_all(dplyr::starts_with("is_signif"), ~.x %in% TRUE)) |>
    dplyr::left_join(df_ptn_map,by=c("protein"="original"),relationship="many-to-many")
  #Supress dplyr warning: many-to-many is expected because may have multiple PTM and total with the same base name

  #Extract Significant correlations and names
  sig<-cor_df[!is.na(cor_df$final_names),]
  original_sig_names<-unique(sig$protein)
  original_pos_names<-unique(sig |>
                               dplyr::filter(dplyr::if_all(dplyr::starts_with("sign_"), ~.x %in% 'pos')) |> pull(protein))
  original_neg_names<-unique(sig |>
                               dplyr::filter(dplyr::if_all(dplyr::starts_with("sign_"), ~.x %in% 'neg')) |> pull(protein))

  cleaned_sig_names<-unique(sig$final_names)
  cleaned_pos_names<-unique(sig |>
                              dplyr::filter(dplyr::if_all(dplyr::starts_with("sign_"), ~.x %in% 'pos')) |> pull(final_names))
  cleaned_neg_names<-unique(sig |>
                              dplyr::filter(dplyr::if_all(dplyr::starts_with("sign_"), ~.x %in% 'neg')) |> pull(final_names))

  return(list(mapping_df=cor_df,
              ptn_class=list(
                original=list(
                  all=original_sig_names,
                  pos=original_pos_names,
                  neg=original_neg_names),
                repository=list(
                  all=cleaned_sig_names,
                  pos=cleaned_pos_names,
                  neg=cleaned_neg_names))
  ))
}
