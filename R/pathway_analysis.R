### Pathway analysis ###

#' Function for PFG analysis
#' @export
#' @examples
#' ## Placeholder Example ##
pfg_path<-function(proteins,pfg_list){

  if (length(proteins) == 0)
    stop("proteins cannot be empty")

  if (!is.list(pfg_list))
    stop("pfg_list must be a list of pathways and their respective proteins")

  if (is.null(names(pfg_list)))
    stop("pfg_list must have pathway names")

  all_ptns<-unique(unlist(pfg_list)) #Define the total number of proteins in the dataset

  #Create list of pathways where selected proteins overlap
  pfg_gp_ptn<-list()
  for (i in seq_along(pfg_list)) {
    pfg_gp_ptn[[ i ]]<-pfg_list[[ i ]][pfg_list[[ i ]] %in% proteins]
  }
  names(pfg_gp_ptn)<-names(pfg_list)

  #Table pathway results
  tab_list<-list()
  for (i in seq_along(pfg_list)) {
    tab_temp<-data.frame(list=c("selected_ptn","not_selected_ptn"),
                         pathway=c( length(pfg_gp_ptn[[ i ]])
                                    ,( length(pfg_list[[ i ]]) - length(pfg_gp_ptn[[ i ]]) )
                         ),
                         not_pathway=c( ( length(proteins)-length(pfg_gp_ptn[[ i ]]) )
                                        ,( length(all_ptns) - length(proteins) - length(pfg_list[[ i ]]) + length(pfg_gp_ptn[[ i ]]) )
                         )
    )
    tab<-as.matrix(tab_temp|> select(pathway,not_pathway))
    rownames(tab)<-tab_temp$list
    tab_list[[ i ]]<-tab
  }
  names(tab_list)<-names(pfg_list)

  #Calculate statistics and report to PFG table
  colnames_pfg_df<-c('Pathway','Library','# Proteins in Pathway','# Sig. Correlated Proteins'
                     ,'% Correlated','pvalue','lowerCI','upperCI'
                     ,'qvalue','Significance','Odds_Ratio','Combined_Score'
                     ,'Correlated_Proteins','All PFG Proteins'
  )
  pfg_df<-data.frame(matrix(NA,nrow=length(pfg_list),ncol=length(colnames_pfg_df)))
  colnames(pfg_df)<-colnames_pfg_df

  for (i in seq_along(pfg_list)) {
    ftest<-fisher.test(tab_list[[ i ]],alternative="two.sided",conf.int=T) #fisher.test on the 2x2 contingency table
    pfg_df[[ i , '# Proteins in Pathway']]<-length(pfg_list[[ i ]])
    pfg_df[[ i , '# Sig. Correlated Proteins']]<-length(pfg_gp_ptn[[ i ]])
    pfg_df[[ i , '% Correlated']]<-round((length(pfg_gp_ptn[[ i ]])/length(pfg_list[[ i ]]))*100,digits=1)
    pfg_df[[ i , 'pvalue']]<-ftest$p
    pfg_df[[ i , 'lowerCI']]<-ftest$conf.int[1]
    pfg_df[[ i , 'upperCI']]<-ftest$conf.int[2]
    pfg_df[[ i , 'Odds_Ratio']]<-(1.0 * tab_list[[ i ]][1,1] * tab_list[[ i ]][2,2]) / max(1.0 * tab_list[[ i ]][2,1] * tab_list[[ i ]][1,2], 1) #manual OR to avoid Inf when zero cells occur
    pfg_df[[ i , 'Combined_Score']]<-(-log(pfg_df[[ i , 'pvalue']]))*(pfg_df[[ i ,'Odds_Ratio']])
    pfg_df[[ i , 'Correlated_Proteins']]<-paste(pfg_gp_ptn[[ i ]],collapse= "; ")
    pfg_df[[ i , 'All PFG Proteins']]<-paste(pfg_list[[ i ]],collapse= ", ")
  }
  pfg_df$Pathway<-names(pfg_list)
  pfg_df$Library<-rep('PFG',nrow(pfg_df))
  pfg_df$lowerCI[!is.finite(pfg_df$lowerCI)]<-NA #turns Inf values into NA for excel handling
  pfg_df$upperCI[!is.finite(pfg_df$upperCI)]<-NA #turns Inf values into NA for excel handling
  pfg_df$qvalue<-p.adjust(pfg_df$pvalue,"BH",n=length(pfg_df$pvalue))
  pfg_df$Significance<-factor(ifelse(pfg_df$qvalue<0.05,'q-value<0.05','q-value>=0.05'))
  pfg_df$Pathway<-factor(pfg_df$Pathway,levels=unique(pfg_df$Pathway[order(pfg_df$Combined_Score,decreasing=T)]))
  pfg_df<-pfg_df[order(pfg_df$pvalue,decreasing=F),]
  pfg_df<-pfg_df[order(pfg_df$qvalue,decreasing=F),]
  pfg_df<-pfg_df[order(pfg_df$Odds_Ratio,decreasing=T),]
  pfg_df<-pfg_df[order(pfg_df$Combined_Score,decreasing=T),]
  rownames(pfg_df)<-NULL

  return(pfg_df)
}

#' Function to get library sets from EnrichR
#' @export
#' @examples
#' ## Placeholder Example ##
get_libs_sets<-function() {
  list(default=c("Reactome_Pathways_2024","GO_Molecular_Function_2025"
                 ,"WikiPathways_2024_Human","KEGG_2026"),
       exploratory=c("Reactome_Pathways_2024","GO_Molecular_Function_2025"
                     ,'GO_Cellular_Component_2025','GO_Biological_Process_2025'
                     ,"WikiPathways_2024_Human","KEGG_2026","MSigDB_Hallmark_2024"),
       react="Reactome_Pathways_2024",
       wiki="WikiPathways_2024_Human",
       kegg="KEGG_2026",
       gobp="GO_Biological_Process_2025",
       gomf="GO_Molecular_Function_2025",
       gocc="GO_Cellular_Component_2025"
  )
}

#' Function wrapper to adjust library argument for EnrichR
#' @export
#' @examples
#' ## Placeholder Example ##
get_enrichr_libs<-function(names="default") {

  libs_sets<-get_libs_sets() #defined library sets
  all_enrichr_libs<-enrichR::listEnrichrDbs()$libraryName #all library sets from EnrichR

  #Convert keyword(s) to actual sets
  libs<-unlist(lapply(names,function(n) {
    if ( n %in% names(libs_sets)) {
      libs_sets[[ n ]]            #predefined set
    } else {
      n                          #keep as it is (user-specified library)
    }
  }))

  #Remove duplicates
  libs<-unique(libs)

  #Validate against actual EnrichR DBs
  invalid_libs<-setdiff(libs,all_enrichr_libs)
  if (length(invalid_libs) > 0) {
    warning("The following libraries are not in EnrichR and will be ignored: ",
            paste(invalid_libs,collapse=", "))
    libs<-setdiff(libs,invalid_libs)
  }
  return(libs)
}

#' Function for EnrichR analysis
#' @export
#' @examples
#' ## Placeholder Example ##
enrichr_path<-function(proteins,libs,enrichr_filt_ptn=FALSE,min_ptn_enrichr_path=3
                       ,enrichr_filt_path=FALSE,min_path_enrichr_ptn=3,out_rm=TRUE
                       ,drop_paths=TRUE,path_number=20){

  #EnrichR libraries
  libs<-get_enrichr_libs(libs)

  #Upload protein names to EnrichR and adjust dataset
  enrichr_df<-enrichr(proteins,libs)
  for (i in 1:length(enrichr_df)){
    enrichr_df[[ i ]]$Library<-rep(names(enrichr_df)[ i ],nrow(enrichr_df[[ i ]]))
  }
  enrichr_df_full<-do.call(rbind,enrichr_df)
  rownames(enrichr_df_full)<-NULL
  enrichr_df_full<-enrichr_df_full[,!(names(enrichr_df_full) %in% c("Old.P.value","Old.Adjusted.P.value"))]
  colnames(enrichr_df_full)<-c('Pathway','Overlap_ptn','pvalue','qvalue','Odds_Ratio'
                               ,'Combined_Score','Correlated_Proteins','Library')

  #(optional) Filter Minimum protein per pathway (minimum number of protein in each protein path)
  if(isTRUE(enrichr_filt_ptn)){
    enr_filt_ptn<-enrichr_df_full
    enr_filt_ptn$N_overlaps<-as.numeric(sub("/.*","",enr_filt_ptn$Overlap_ptn))
    enr_filt_ptn<-enr_filt_ptn[enr_filt_ptn$N_overlaps>=min_ptn_enrichr_path,]
  } else {
    enr_filt_ptn<-enrichr_df_full
  }

  #(optional) Filter Minimum links per protein (minimum number of pathways each protein appears in)
  if(isTRUE(enrichr_filt_path)){
    ptn_path<-enr_filt_ptn[,c('Pathway','Correlated_Proteins')]|> tidyr::separate_rows(Correlated_Proteins)
    ptn_path_cts<-ptn_path|> dplyr::distinct(Correlated_Proteins,Pathway)|> dplyr::count(Correlated_Proteins,name="ptn_count")
    filt_ptn<-ptn_path_cts$Correlated_Proteins[ptn_path_cts$ptn_count>=min_path_enrichr_ptn]
    enr_filt_ptn_path<-enr_filt_ptn[sapply(strsplit(enr_filt_ptn$Correlated_Proteins,";"),function(x) any(x %in% filt_ptn)),]
  } else {
    enr_filt_ptn_path<-enr_filt_ptn
  }

  #(optional) Clean unecessary Library pathways
  default_drop<-c("syndrome","deletion","duplication","copy\\s*number","\\bcnv\\b",
                  "trisomy","microdeletion","microduplication","chromosome|chromosomal",
                  "\\bdisease\\b","infection","viral|virus","bacterial|bacteria",
                  "cancer|carcinoma|leukemia|melanoma|mesothelioma|tumou?r",
                  "sars","covid","coronavirus","hiv","covs",
                  "influenza","hepatitis","ebola","zika",
                  "ectoderm|endoderm|mesoderm","mesonephric|metanephric|nephric",
                  "embryo|embryonic|embryogenesis","fetal|foetal|fetus|foetus",
                  "gastrulat","trophoblast","placenta","morphogenesis",
                  "somite|somitogenesis","organogenesis")
  if (is.logical(drop_paths)) {
    #TRUE=use default drop list, FALSE=drop nothing
    paths_to_drop<-if (drop_paths) default_drop else character(0)
  } else if (is.character(drop_paths)) {
    #user provides custom list of patterns
    paths_to_drop<-drop_paths
  } else {
    stop("\ndrop_paths argument must be TRUE, FALSE, or a character vector.\n")
  }
  if(length(paths_to_drop) > 0){
    enr_filt_ptn_path_clean<-enr_filt_ptn_path[!grepl(paste(paths_to_drop,collapse="|"),tolower(enr_filt_ptn_path$Pathway)),]
  } else {
    enr_filt_ptn_path_clean<-enr_filt_ptn_path
  }

  #(optional) Remove outlier paths
  if(isTRUE(out_rm)){
    #Calculate Q1, Q3, and IQR
    enrich_q1<-quantile(enr_filt_ptn_path_clean$Combined_Score,0.25)
    enrich_q3<-quantile(enr_filt_ptn_path_clean$Combined_Score,0.75)
    enrich_iqr<-IQR(enr_filt_ptn_path_clean$Combined_Score)
    #Define bounds
    lower_bound_enrichr<-enrich_q1 - 1.5 * enrich_iqr
    upper_bound_enrichr<-enrich_q3 + 1.5 * enrich_iqr
    #Identify and remove outliers
    outliers_enrichr<-subset(enr_filt_ptn_path_clean$Pathway,enr_filt_ptn_path_clean$Combined_Score < lower_bound_enrichr
                             | enr_filt_ptn_path_clean$Combined_Score > upper_bound_enrichr)
    enr_filt_ptn_path_clean_out_rm<-enr_filt_ptn_path_clean[!grepl(paste(outliers_enrichr,collapse="|"),enr_filt_ptn_path_clean$Pathway),]
  } else {
    enr_filt_ptn_path_clean_out_rm<-enr_filt_ptn_path_clean
  }

  enrichr_df_final<-enr_filt_ptn_path_clean_out_rm|>
    separate(Overlap_ptn,into=c("# Sig. Correlated Proteins","# Proteins in Pathway"),sep="/")
  enrichr_df_final$"# Sig. Correlated Proteins"<-as.numeric(enrichr_df_final$"# Sig. Correlated Proteins")
  enrichr_df_final$"# Proteins in Pathway"<-as.numeric(enrichr_df_final$"# Proteins in Pathway")
  enrichr_df_final$'% Correlated'<-round((enrichr_df_final$"# Sig. Correlated Proteins"/enrichr_df_final$"# Proteins in Pathway")*100,digits=1)
  enrichr_df_final$Significance<-factor(ifelse(enrichr_df_final$qvalue<0.05,'q-value<0.05','q-value>=0.05'))
  enrichr_df_final<-enrichr_df_final|> select(c('Pathway','Library','# Proteins in Pathway'
                                                  ,'# Sig. Correlated Proteins','% Correlated'
                                                  ,'pvalue','qvalue','Significance','Odds_Ratio'
                                                  ,'Combined_Score','Correlated_Proteins'))
  enrichr_df_final<-enrichr_df_final[order(enrichr_df_final$pvalue,decreasing=F),]
  enrichr_df_final<-enrichr_df_final[order(enrichr_df_final$qvalue,decreasing=F),]
  enrichr_df_final<-enrichr_df_final[order(enrichr_df_final$Odds_Ratio,decreasing=T),]
  enrichr_df_final<-enrichr_df_final[order(enrichr_df_final$Combined_Score,decreasing=T),]
  rownames(enrichr_df_final)<-NULL
  return(enrichr_df_final)
}

#' Function to generate Enrichment Barplots
#' @export
#' @examples
#' ## Placeholder Example ##
plot_enrichment_paths<-function(df,type=c('PFG','EnrichR'),path_number=20){

  type<-match.arg(type)

  if (!"Combined_Score" %in% names(df))
    stop("Combined_Score column missing from dataframe. Please check input")

  #Define minimal cutoff for combined Score (has to be >0 otherwise should not plot)
  df_valid<-df[df$Combined_Score > 0, ]
  df_valid<-df_valid[order(df_valid$Combined_Score,decreasing=TRUE), ]

  if (nrow(df_valid)==0)
    stop("No positive Combined_Score values to plot.") #If no Combined_Score is > 0 then there is an issue. This warns that something is not right

  #Define number of paths to show
  path<-min(nrow(df_valid),path_number) #Extract what is smaller: number of paths with valid scores of the user-defined path number

  if (type=='PFG'){
    barcolors<-brewer.pal(4,'Pastel1')[4]
    legend_position<-'inside'
  } else {
    barcolors<-brewer.pal(5,'Pastel1')[c(3,2,1,5)]
    legend_position<-'bottom'
    }

  plot_paths<-ggplot(df_valid[1:path,],aes(x=reorder(fct_rev(Pathway),Combined_Score),y=Combined_Score,fill=Library))+
    ggtitle(paste('Top',type,'Pathways'))+
    geom_bar(color='black',position=position_dodge(width=0.7),stat="identity",linewidth=0.1)+
    geom_text(aes(label=Pathway),vjust=0.5,y=0.25,hjust=0,color="black",size=3,fontface='bold')+
    #geom_text(aes(label=Correlated_Proteins),vjust=1.5,y=0.25,hjust=0,color="black",size=2)+
    labs(y='Combined Score',x=NULL)+scale_fill_manual(name="Library",values=barcolors)+
    theme(plot.title=element_text(size=12,hjust=0,face='bold',color='black',margin=margin(t=1,r=0,b=5,l=0))
          ,axis.text.x=element_text(size=12,face="bold",color='black',margin=margin(t=1,r=0,b=0,l=0))
          ,axis.title.x=element_text(size=12,face="bold",margin=margin(t=2,r=0,b=0,l=0))
          ,axis.line=element_line(linewidth=0.25,linetype="solid",colour="black")
          ,axis.ticks=element_line(linewidth=0.25)
          ,axis.ticks.length=unit(0.025,"cm")
          ,aspect.ratio=1/1
          ,legend.position=legend_position
          ,legend.title=element_text(face="bold")
          ,legend.position.inside=c(0.8,0.25)
          ,legend.text=element_text(face="bold")
          ,legend.key.spacing.y=unit(0.1,"lines")
          ,legend.margin=margin(t=0,r=0,b=0,l=0,unit="cm")
          ,legend.key.size=unit(0.5,"lines")
          ,legend.box.margin=margin(t=-0.25,r=0,b=0,l=0,unit="cm")
          ,axis.text.y=element_blank()
          ,axis.ticks.y=element_blank()
          ,axis.title.y=element_blank()
          ,panel.background=element_blank()
          ,panel.grid.major=element_blank()
          ,panel.grid.minor=element_blank()
          ,panel.border=element_blank()
          )+
    guides(fill=guide_legend(byrow=TRUE))+coord_flip()
  return(plot_paths)
}

#' Function to plot Heatmap of pathway Combined scores
#' @export
#' @examples
#' ## Placeholder Example ##
path_heatmap<-function(df_list,libs='default',analysis=c("PFG","EnrichR")
                       ,legend_title="Combined\nScore",legend_direction='vertical'
                       ,cluster_columns=TRUE,anno_row=TRUE,colors_ht=NULL,plot_title=NULL
                       ,row_fontsize=10,col_fontsize=10,title_fontsize=12,max_str_width=40
                       ,n_legend_breaks=5,n_color_breaks=11,score_cutoff=1,top_n_paths=20
                       ,colors_anno_row=c(pal_jco()(2)[2:1],'red2',brewer.pal(7,"Dark2")[c(7,1,2:6)])){

  analysis<-match.arg(analysis)

  #Validate inputs
  if (length(df_list) < 2)
    stop("df_list must contain at least 2 datasets to compare.")

  if (analysis %in% "PFG") {
    libs<-"PFG"
  } else {
    libs<-get_enrichr_libs(libs)
  }

  #Combine datasets
  df_merged<-df_list |>
    purrr::imap(~ .x |> dplyr::filter(Library %in% libs) |>
                  dplyr::select(Pathway,Combined_Score) |>
                  dplyr::rename(!!paste0("score_", .y) := Combined_Score)
    ) |> purrr::reduce(full_join,by="Pathway") |>
    dplyr::mutate(across(starts_with("score_"), ~replace_na(.x, 0)))

  #Helper to add overlap label as row annotation
  make_overlap_label<-function(row,col_names) {
    present<-col_names[row > 0]
    present<-sub("^score_","",present)

    if (length(present)==0) {
      return("None")
    } else if (length(present)==1) {
      return(paste0(present," only"))
    } else {
      return(paste(sort(present),collapse=" & "))
    }
  }

  #Add overlap label and filter out paths with no overlaps
  score_cols<-grep("^score_",names(df_merged),value=TRUE)
  df_merged<-df_merged |> dplyr::rowwise() |> dplyr::mutate(
    overlaps = make_overlap_label(dplyr::c_across(dplyr::all_of(score_cols)),score_cols)) |>
    dplyr::ungroup() |> dplyr::mutate(
      n_sets   = rowSums(dplyr::across(starts_with("score_")) > 0),
      overlaps = factor(overlaps)) |> dplyr::arrange(dplyr::desc(n_sets)) |>
    dplyr::mutate(overlaps = factor(overlaps,levels=unique(overlaps))) |>
    dplyr::select(-n_sets)

  #Filter by Combined Score Cutoff and the ones without any overlap
  df_filtered<-df_merged |> dplyr::filter(overlaps !='None') |>
    dplyr::filter(apply(dplyr::pick(starts_with("score_")),1,max) >= score_cutoff) |>
    droplevels()

  if (nrow(df_filtered)==0)
    stop ("No pathways remain after filtering. Consider lowering score_cutoff.")

  #Order df by overlaps and then highest score
  mat<-as.matrix(df_filtered |> dplyr::select(starts_with("score_")))
  rownames(mat)<-df_filtered$Pathway
  mat[is.na(mat)]<-0
  row_max<-apply(mat,1,max)
  n_present<-rowSums(mat > 0)
  order_idx<-names(sort(row_max,decreasing=TRUE))
  order_idx<-order_idx[order(n_present[order_idx],decreasing=TRUE)]
  df_top<-df_filtered[match(order_idx,df_filtered$Pathway),] |>
    dplyr::slice_head(n=top_n_paths) |> droplevels()

  #Adjust row names in case of EnrichR
  if (analysis %in% "EnrichR") {
    df_top$Pathway<-gsub("\\s*\\([^\\)]+\\)","",df_top$Pathway)
    df_top$Pathway<-gsub("R-HSA-.*","",df_top$Pathway)
    df_top$Pathway<-gsub(" WP.*","",df_top$Pathway)
    df_top$Pathway<-stringr::str_trunc(df_top$Pathway,width=max_str_width,side="right")
    df_top<-dplyr::distinct(df_top,Pathway,.keep_all=TRUE)
  }

  #Add row annotation if called
  if(isTRUE(anno_row)){
    if (nlevels(df_top$overlaps) > length(colors_anno_row)) {
      colors_anno_row<-grDevices::colorRampPalette(colors_anno_row)(nlevels(df_top$overlaps))
    }
    cols_row_anno<-list(Overlaps=stats::setNames(
      colors_anno_row[1:nlevels(df_top$overlaps)],levels(df_top$overlaps)))
    row_anno_ht<-ComplexHeatmap::rowAnnotation(
      Overlaps             = df_top$overlaps,
      col                  = cols_row_anno,
      annotation_name_rot  = 45,
      annotation_name_side = "bottom")
  } else {
    row_anno_ht<-NULL
  }

  #Heatmap matrix
  mtx<-as.matrix(df_top |> dplyr::select(starts_with("score_")))
  rownames(mtx)<-df_top$Pathway
  mtx_plot<-log1p(mtx) #always log transform

  #Adjust heatmap colors
  if (is.null(colors_ht)) {
    breaks_log<-seq(min(mtx_plot),max(mtx_plot),length.out=n_color_breaks)
    cols_ht<-circlize::colorRamp2(breaks_log,matlab::jet.colors(n_color_breaks))
  } else {
    breaks_log<-seq(min(mtx_plot),max(mtx_plot),length.out=n_color_breaks)
    cols_ht<-colors_ht
  }

  #Adjust legend breaks and labels (heatmap is in log, but legend show linear vals)
  legend_idx<-round(seq(1,n_color_breaks,length.out=n_legend_breaks))
  legend_breaks<-breaks_log[legend_idx]
  legend_labels<-round(expm1(legend_breaks))

  #Plot heatmap
  path_ht<-ComplexHeatmap::Heatmap(mtx_plot,
                   row_names_gp=gpar(fontsize=row_fontsize),
                   column_names_gp=gpar(fontsize=col_fontsize),
                   column_title_gp=gpar(fontsize=title_fontsize,fontface="bold"),
                   name=legend_title,
                   column_labels=names(df_list),
                   show_row_dend=FALSE,
                   show_column_dend=FALSE,
                   column_title=plot_title,
                   column_names_rot=45,
                   left_annotation=row_anno_ht,
                   cluster_rows=FALSE,
                   cluster_columns=cluster_columns,
                   col=cols_ht
                   ,heatmap_legend_param=list(
                     at=legend_breaks
                     ,labels=legend_labels
                     ,direction=legend_direction
                   )
  )
  return(list(plot=path_ht,table=df_merged))
}

