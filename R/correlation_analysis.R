### Correlation analysis core functions ###

#' Function to perform correlation analysis
#' @export
#' @examples
#' ## Placeholder Example ##
ptn_corr<-function(df,all_proteins,proteins_of_interest,method='spearman'
                   ,conf_level=0.95,pval_cutoff=0.05,abs_cutoff=0.20
                   ,p_adjust_method="BH"){

  #Subset dataframe for all proteins
  df<-df[,all_proteins,drop=F]
  df_mtx<-as.matrix(df)

  #Stop in case proteins of interest are not in main data matrix
  missing<-setdiff(proteins_of_interest,colnames(df_mtx))
  if(length(missing) > 0)
    stop("Proteins not found: ",paste(missing,collapse=", "))

  #Make sure proteins of interest are in main data matrix
  proteins_of_interest<-intersect(proteins_of_interest,colnames(df_mtx))

  #Rank in case of spearman corr
  if(method=="spearman")  df_mtx<-apply(df_mtx,2,function( x ) rank( x, na.last="keep"))
  df_mtx<-as.matrix(df_mtx)
  colnames(df_mtx)<-colnames(df)

  #Full correlation matrix
  cor_full<-cor(df_mtx,use="pairwise.complete.obs")

  #Subset to proteins of interest
  cor_mtx<-cor_full[proteins_of_interest,,drop=FALSE]

  #Pairwise sample size matrix
  valid<-!is.na(df_mtx)
  n_matrix<-t(valid) %*% valid
  n_mtx<-n_matrix[proteins_of_interest,,drop=FALSE]

  #Sanity check
  if(!identical(dim(cor_mtx),dim(n_mtx)))
    stop("Correlation and sample size matrices have different dimensions")

  #Compute t statistic
  tval<- cor_mtx * sqrt(( n_mtx -2) / pmax(1 - cor_mtx ^2,.Machine$double.eps))

  #Compute pvalues
  pval<- 2 * pt( -abs( tval ),df= n_mtx -2)

  #Confidence intervals
  z_alpha<-qnorm(1-(1- conf_level )/2)
  z<-atanh( cor_mtx )
  se<-1/sqrt(pmax( n_mtx -3, 1))
  ci_lower<-tanh(z - z_alpha * se)
  ci_upper<-tanh(z + z_alpha * se)

  #Pivot to long format
  grid<-expand.grid(query_protein=rownames(cor_mtx),protein=colnames(cor_mtx)
                    ,KEEP.OUT.ATTRS=FALSE,stringsAsFactors=FALSE)

  cor_df<-data.frame(grid,corr=as.vector(cor_mtx),tval=as.vector(tval),
                     pval=as.vector(pval),lowerCI=as.vector(ci_lower),upperCI=as.vector(ci_upper))

  #Remove self-correlations
  cor_df<-cor_df[cor_df$query_protein !=cor_df$protein,]

  #Adjust pvalues
  cor_df$qval<-p.adjust(cor_df$pval,method=p_adjust_method)

  #Add extra columns
  cor_df$abs_corr<-abs(cor_df$corr)
  cor_df$sign<-ifelse(cor_df$corr>0,"pos","neg")
  cor_df$is_signif<-cor_df$abs_corr >= abs_cutoff & cor_df$qval <= pval_cutoff

  #Sort results
  cor_df<-cor_df[order(cor_df$query_protein,-cor_df$abs_corr),]
  rownames(cor_df)<-NULL

  return(cor_df)
}

#' Function to generate Waterfall plots
#' @export
#' @examples
#' ## Placeholder Example ##
wf_plot<-function(df,type=c('Positive','Negative'),abs_cutoff=0.20){

  type<-match.arg(type)

  if (!"corr" %in% names(df))
    stop("corr column missing from dataframe. Please check input")

  if (type=='Positive'){
    barcol<-'red'
    position<-'bottom'
    sign<-'pos'
    df_valid<-df[df$abs_corr > abs_cutoff & df$is_signif & df$sign==sign,]
    df_valid<-df_valid[order(df_valid$corr,decreasing=FALSE),]
    df_valid$protein<-factor(df_valid$protein,levels=df_valid$protein)
  }

  if (type=='Negative'){
    barcol<-'blue'
    position<-'top'
    sign<-'neg'
    df_valid<-df[df$abs_corr > abs_cutoff & df$is_signif & df$sign==sign,]
    df_valid<-df_valid[order(df_valid$corr,decreasing=TRUE),]
    df_valid$protein<-factor(df_valid$protein,levels=df_valid$protein)
  }

  if (nrow(df_valid)==0) {
    return(ggplot()+theme_void()+ggtitle(paste(type,"correlations: none")))
  }

  wfplot<-ggplot(aes(x=protein,y=corr,fill=sign),data=df_valid)+
    geom_bar(position='dodge',stat="identity")+scale_x_discrete(position=position)+
    theme_classic(base_size=11,base_line_size=0.05,base_rect_size=0.05)+
    scale_fill_manual(values=barcol)+ylab(paste(type,"Correlation"))+xlab(NULL)+
    theme(text=element_text(size=4),legend.position="none"
          ,axis.text.x=element_text(size=6,angle=0,hjust=1)
          ,axis.title.x=element_text(size=7,angle=0,hjust=1)
          ,plot.title=element_text(size=8,margin=margin(t=0,r=0,b=0,l=0,unit="cm"))
          ,plot.margin=margin(t=0.5,r=0.5,b=0.5,l=0.5,unit="cm"))+coord_flip()

  return(wfplot)
}

#' Function to generate scatter plots
#' @export
#' @examples
#' ## Placeholder Example ##
plot_corr_scatter<-function(comb_corr,datasetA,datasetB,datasetA_name,datasetB_name
                            ,query_protein,highlight_n=20
                            ,colors=c('red2','blue2','purple2','grey50')){

  #Helper functions for Z transformation
  fisher_z<-function(r) 0.5 * log((1 + r) / (1 - r))
  se_z<-function(n) 1 / sqrt(n - 3)

  df<-comb_corr |>
    select(query_protein,dataset_label,protein,corr,sign,is_signif) |>
    pivot_wider(names_from=dataset_label,values_from=c(corr,sign,is_signif))

  #Adjust names
  corA<-paste0("corr_",datasetA)
  corB<-paste0("corr_",datasetB)
  is_signif_A<-paste0("is_signif_",datasetA)
  is_signif_B<-paste0("is_signif_",datasetB)
  sign_A<-paste0("sign_",datasetA)
  sign_B<-paste0("sign_",datasetB)

  #Remove NAs
  df<-df |> dplyr::filter(!is.na(.data[[corA]]),!is.na(.data[[corB]]))

  #Z transform data before calculating global similarity
  df<-df |> dplyr::mutate(z_A=fisher_z(.data[[corA]]),z_B=fisher_z(.data[[corB]]))

  #Calculate global similarity between networks
  global_cor<-cor(df$z_A,df$z_B,use="complete.obs",method="pearson")

  #Test if global_cor is significantly different from 0
  n_proteins<-nrow(df)
  z_global<-fisher_z(global_cor)
  se_global<- 1 / sqrt(n_proteins - 3)
  p_global<- 2 * pnorm(-abs(z_global / se_global))

  sim_flag<-dplyr::case_when(
    p_global >= 0.05             ~ " [No significant similarity]",
    global_cor < 0               ~ " [anti-correlation]",
    global_cor < 0.3             ~ " [weak similarity]",
    global_cor < 0.6             ~ " [moderate similarity]",
    global_cor < 0.8             ~ " [strong similarity]",
    TRUE                         ~ " [very strong similarity]"
  )

  #Compute Stats for reports with relatioship flag
  df_lm<-df |> dplyr::rename(x=all_of(corA),y=all_of(corB))
  lm_fit<-lm(y~ x,data=df_lm)
  lm_summary<-summary(lm_fit)
  slope<-coef(lm_fit)[2]
  r2<-lm_summary$r.squared
  pval<-lm_summary$coefficients[2,'Pr(>|t|)']

  #Relationship flag
  attenuation_flag<-dplyr::case_when(
    slope < -0.05                 ~ " [anti-correlation]",
    slope > -0.05 & slope < 0.2   ~ " [no linear agreement]",
    slope >= 0.2  & slope < 0.9   ~ " [attenuation]",
    slope >= 0.9  & slope < 1.0   ~ " [near-perfect agreement]",
    slope == 1.0                  ~ " [perfect agreement]",
    slope >= 1.1                  ~ " [amplification]"
  )

  #Identify rate of concordance between datasets for all proteins (overall concordance)
  df<-df |> dplyr::mutate(concordant=.data[[ sign_A ]]==.data[[ sign_B ]]
                    ,weight=sqrt(abs(.data[[corA]]) * abs(.data[[corB]]))
                    ,strength=abs(.data[[ corA ]]) + abs(.data[[ corB ]]))
  concordance_rate<-weighted.mean(df$concordant,df$weight,na.rm=TRUE)  #mean of a logical var= TRUE %

  #Identify significant corr in each datatset
  df<-df |> dplyr::mutate(sig_A=.data[[ is_signif_A ]]
                    ,sig_B=.data[[ is_signif_B ]]
                    ,Significant=dplyr::case_when(sig_A & sig_B~"Both"
                                           ,sig_A~datasetA_name
                                           ,sig_B~datasetB_name
                                           ,TRUE~"NS"))
  df$Significant<-factor(df$Significant,levels=c(datasetA_name,datasetB_name,"Both","NS"))

  #Correlation between concordant correlations only (concordance among signif in both)
  df_sig<-df |> dplyr::filter(Significant=="Both")                        #subset signif in both
  df_sig_cor<-if (nrow(df_sig) > 2) {
    cor(df_sig$z_A,df_sig$z_B,use="complete.obs",method="pearson")  #cor of signif in both
  } else { NA }
  sig_concordance<-weighted.mean(df_sig$concordant,df_sig$weight,na.rm=TRUE)
  #concordant=logical var (TRUE if sign in A is the same as in B, FALSE otherwise)
  #rate=mean of a logical var (TRUE %), or the % of concordant cases

  #Highlight top concordant proteins
  top<-df |> dplyr::filter(sig_A & sig_B & concordant) |> arrange(desc(strength)) |> slice_head(n=highlight_n)

  #Adjust colors for labelling
  color_vals<-setNames(colors,c(datasetA_name,datasetB_name,"Both","NS"))

  #Plot
  scplot<-ggplot(df,aes(x=.data[[corA]],y=.data[[corB]],color=Significant))+
    geom_point(alpha=0.6)+scale_color_manual(values=color_vals)+
    #geom_text_repel(data=top,aes(label=protein),size=2,color='black',max.overlaps=30)+
    #geom_vline(xintercept=0,linetype="dashed")+geom_hline(yintercept=0,linetype="dashed")+
    geom_smooth(method="lm",se=FALSE,color="darkgreen",linetype="solid",linewidth=0.5)+
    geom_abline(slope=1,intercept=0,linetype="dashed",color="orange3",linewidth=0.5)+
    theme_classic()+theme(plot.title=element_text(size=10)
                          ,plot.subtitle=element_text(size=7.5,face="italic")
                          ,plot.margin=margin(t=0.25,r=0.25,b=0.25,l=0.25,unit="cm")
                          ,legend.title=element_text(size=8)
                          ,legend.text=element_text(size=6.5)
                          ,legend.position="bottom"
                          ,legend.margin=margin(t=0,r=0,b=0,l=0,unit="cm")
                          ,legend.box.margin=margin(t=-0.25,r=0,b=-0.25,l=0,unit="cm")
                          ,legend.key.width=unit(0.25,"cm")
                          ,legend.key.height=unit(0.25,"cm")
                          ,legend.spacing.x=unit(0.25,"cm")
                          )+
    labs(title=paste(datasetA_name,"vs",datasetB_name,"correlations")
         ,x=paste("Correlations",datasetA_name),y=paste("Correlations",datasetB_name)
         ,subtitle=paste0("Global Correlation: R = ",round(global_cor,2)
                          ,", p = ", signif(p_global, 2)
           ,"\nSignificant overlaps N = ",nrow(df_sig)
           )
         )

  return(list(scplot=scplot
              ,top_sig_overlaps=top$protein
              ,stats=list(global_cor=global_cor,p_global=p_global
                          ,slope=slope,r2=r2,pval=pval,relationship=attenuation_flag)))
}

#' Function togenerate delta-delta scatter plots
#' @export
#' @examples
#' ## Placeholder Example ##
plot_delta_scatter<-function(comb_corr,datasetA,datasetB
                               ,baseline_dataset,query_protein
                               ,datasetA_name,datasetB_name
                               ,colors=c('red2','blue2','purple2','grey50')
                               ,baseline_name,highlight_n=20) {

  #Helper functions for Z transformation
  fisher_z<-function(r) 0.5 * log((1 + r) / (1 - r))
  se_z<-function(n) 1 / sqrt(n - 3)

  #Pivot dataset
  df<-comb_corr |> select(query_protein,dataset_label,protein,corr,is_signif,n_samples) |>
    pivot_wider(names_from=dataset_label,values_from=c(corr,is_signif,n_samples),names_sep="_")

  #Adjust names
  corA<-paste0("corr_",datasetA)
  corB<-paste0("corr_",datasetB)
  sigA<-paste0("is_signif_",datasetA)
  sigB<-paste0("is_signif_",datasetB)
  corbase<-paste0("corr_",baseline_dataset)

  #Compute deltas (corr - baseline corr)
  df<-df |> dplyr::mutate(delta_A= .data[[corA]] - .data[[corbase]]
                   ,delta_B= .data[[corB]] - .data[[corbase]])

  #Remove NAs
  df<-df |> dplyr::filter(!is.na(delta_A),!is.na(delta_B))

  #Global delta-delta similarity
  delta_cor<-cor(df$delta_A,df$delta_B,use="complete.obs")

  #Add stats
  lm_fit<-lm(delta_B~delta_A,data=df)
  lm_summary<-summary(lm_fit)
  slope<-coef(lm_fit)[2]
  r2<-lm_summary$r.squared
  pval<-lm_summary$coefficients[2,'Pr(>|t|)']

  #Identify significant corr in each datatset
  df<-df |> dplyr::mutate(sig_A= .data[[sigA]],sig_B= .data[[sigB]]
                    ,Significant=dplyr::case_when(sig_A & sig_B~"Both",sig_A~datasetA_name,
                                           sig_B~datasetB_name,TRUE~"NS"))
  df$Significant<-factor(df$Significant,levels=c(datasetA_name,datasetB_name,"Both","NS"))

  #Identify rate of concordance between datasets for all proteins (overall concordance)
  df<-df |> dplyr::mutate(concordant=sign(delta_A)==sign(delta_B),strength=abs(delta_A)+abs(delta_B))
  concordance_rate<-mean(df$concordant,na.rm=TRUE) #mean of a logical var= TRUE %

  #Correlation between concordant correlations only (concordance among signif in both)
  sig_both_df<-df |> dplyr::filter(Significant=="Both")                    #subset signif in both
  sig_both_cor<-if (nrow(sig_both_df) > 2) {
    cor(sig_both_df$delta_A,sig_both_df$delta_B,use="complete.obs")  #cor of signif in both
  } else { NA }
  sig_concordance<-sig_both_df |> summarise(rate=mean(concordant,na.rm=TRUE)) |> pull(rate)
  #concordant=logical var (TRUE if sign in A is the same as in B, FALSE otherwise)
  #rate=mean of a logical var (TRUE %), or the % of concordant cases

  #Highlight top concordant proteins
  top<-df |> dplyr::filter(sig_A & sig_B & concordant) |> arrange(desc(strength)) |> slice_head(n=highlight_n)

  #Adjust colors
  color_vals<-setNames(colors,c(datasetA_name,datasetB_name,"Both","NS"))

  #Plot
  ddscplot<-ggplot(df,aes(x=delta_A,y=delta_B,color=Significant))+
    geom_point(alpha=0.6)+scale_color_manual(values=color_vals)+
    geom_abline(slope=1,intercept=0,linetype="dashed",color="orange3",size=0.5)+
    #geom_text_repel(data=top,aes(label=protein),size=2,max.overlaps=30,color='black')+
    geom_vline(xintercept=0,linetype="dashed")+geom_hline(yintercept=0,linetype="dashed")+
    geom_smooth(method="lm",se=FALSE,color="darkgreen",linetype="solid",size=0.5)+theme_classic()+
    theme(plot.title=element_text(size=10)
          ,plot.margin=margin(t=0.25,r=0.25,b=0.25,l=0.25,unit="cm")
          ,plot.subtitle=element_text(size=7.5,face="italic")
          ,legend.title=element_text(size=8)
          ,legend.text=element_text(size=6.5)
          ,legend.position="bottom"
          ,legend.margin=margin(t=0,r=0,b=0,l=0,unit="cm")
          ,legend.box.margin=margin(t=-0.25,r=0,b=-0.25,l=0,unit="cm")
          ,legend.key.width=unit(0.25,"cm")
          ,legend.key.height=unit(0.25,"cm")
          ,legend.spacing.x=unit(0.25,"cm")
          )+
    labs(title=paste(datasetA_name,"vs",datasetB_name,"(Delta vs Delta)")
         ,x=bquote(Delta~Corr~"(" * .(datasetA_name) * " - " * .(baseline_name) * ")")
         ,y=bquote(Delta~Corr~"(" * .(datasetB_name) * " - " * .(baseline_name) * ")")
      ,subtitle=paste0("Delta-Delta similarity: R = ",round(delta_cor,2)
                       ," | Significant overlaps N =",nrow(sig_both_df)
                       ,"\nConcordance in overlaps: R = ",round(sig_both_cor,2)
                       ," | Rate: ",round(sig_concordance*100,1),"%"
                       ,"\n(R-square = ",round(r2,2)
                       ,", slope = ",round(slope,2)
                       ,", p = ",signif(pval,2),")"
                       )
      )

  return(list(scplot=ddscplot,top_sig_overlaps=top$protein))
}

#' Function to generate correlation change heatmap by protein
#' @export
#' @examples
#' ## Placeholder Example ##
plot_corr_heatmap<-function(comb_corr,subsets,dataset_names,title_suffix=NULL
                            ,top_n=50,cluster_rows=TRUE,signif_mode=c("any","all")
                            ,annotate_signif=TRUE,signif_symbol="\u2731",subtitle=TRUE
                            ,switch_mode=c("none","dissimilar","prioritize_dissimilar"
                                           ,"similar","prioritize_similar")) {
  #Match arguments
  signif_mode<-match.arg(signif_mode)
  switch_mode<-match.arg(switch_mode)

  #Summarise per protein
  top_df<-comb_corr |> dplyr::filter(dataset_label %in% subsets) |>
    dplyr::group_by(protein) |> dplyr::summarise(
      signif_count=sum(is_signif,na.rm=TRUE), #count datasets that protein is significant
      signif_flag=if (signif_mode=="any") {
        any(is_signif,na.rm=TRUE) #TRUE if a protein is signif in at least one dataset
      } else {
        all(is_signif,na.rm=TRUE) #TRUE if a protein is signif in all datasets
        },change=max(corr,na.rm=TRUE) - min(corr,na.rm=TRUE)
      ,sign_switch=dplyr::n_distinct(sign[!is.na(sign)]) > 1, .groups="drop") |>
    dplyr::filter(signif_flag)

  #Apply switch_mode
  if (switch_mode=="dissimilar") {
    #Keep only proteins that change sign (dissimilar) across datasets
    top_df<-top_df |> dplyr::filter(sign_switch)

  } else if (switch_mode=="prioritize_dissimilar") {
    #Keep all proteins, but rank sign-changing proteins first, then by magnitude of change
    top_df<-top_df |> dplyr::arrange(dplyr::desc(signif_count),dplyr::desc(sign_switch), dplyr::desc(change))

  } else if (switch_mode=="similar") {
    #Only similar proteins
    top_df<-top_df |> dplyr::filter(!sign_switch)

  } else if (switch_mode=="prioritize_similar") {
    #Similar first
    top_df<-top_df |>dplyr::arrange(dplyr::desc(signif_count),sign_switch,change)

  } else {
    #Default
    top_df<-top_df |> dplyr::arrange(dplyr::desc(change))
  }

  #Select top proteins
  top_ptn<-top_df |> dplyr::slice_head(n=top_n) |> dplyr::pull(protein)

  #Filter original data
  top_comb_corr<-comb_corr |> dplyr::filter(protein %in% top_ptn)

  #Annotate significant proteins in each dataset
  if (annotate_signif) {
    signif_long<-top_comb_corr |>
      dplyr::filter(dataset_label %in% subsets) |>
      dplyr::select(protein,dataset_label,is_signif) |>
      dplyr::group_by(protein,dataset_label) |>
      dplyr::summarise(is_signif=any(is_signif,na.rm=TRUE),.groups="drop")

    #Map dataset_label: display name using subsets/dataset_names vectors
    label_map<-stats::setNames(dataset_names,subsets)
    signif_long<-signif_long |> dplyr::mutate(dataset=dplyr::recode(dataset_label,!!!label_map))
  }

  #Prepare heatmap data
  heat_df<-top_comb_corr |> dplyr::filter(dataset_label %in% subsets) |>
    dplyr::select(protein,dataset_label,corr) |>
    tidyr::pivot_wider(names_from=dataset_label,values_from=corr) |>
    dplyr::select(protein,dplyr::all_of(subsets))
  colnames(heat_df)<-c("protein",dataset_names)

  #(optional) clustering
  if (cluster_rows) {
    mat<-as.data.frame(heat_df)
    rownames(mat)<-mat$protein
    mat<-as.matrix(mat[,-1])
    mat[is.na(mat)]<-0  #NA handling

    row_clust<-stats::hclust(stats::dist(mat))
    protein_order<-rownames(mat)[row_clust$order]

    #Force order direction: high to low
    avg_corr<-rowMeans(mat,na.rm=TRUE)
    if (avg_corr[protein_order[1]] < avg_corr[protein_order[length(protein_order)]]) {
      protein_order<-rev(protein_order)
    }

  } else {
    protein_order<-top_df |> dplyr::filter(protein %in% top_ptn) |> dplyr::pull(protein)
  }

  #Back to long format
  heat_long<-heat_df |> tidyr::pivot_longer(-protein,names_to="dataset",values_to="corr")

  #Apply ordering (reverse for ggplot top-down)
  heat_long$protein<-factor(heat_long$protein,levels=rev(protein_order))
  heat_long$dataset<-factor(heat_long$dataset,levels=dataset_names)

  #Join significance on plotting data
  if (annotate_signif) {
    signif_long$dataset<-factor(signif_long$dataset,levels=dataset_names)
    signif_long$protein<-factor(signif_long$protein,levels=rev(protein_order))

    heat_long<-heat_long |>
      dplyr::left_join(
        signif_long |> dplyr::select(protein,dataset,is_signif),by=c("protein", "dataset")) |>
      dplyr::mutate(signif_label=dplyr::if_else(is_signif,signif_symbol,""))
  }

  #Plotting
  ht<-ggplot(heat_long,aes(dataset,protein,fill=corr))+geom_tile()+
    scale_fill_gradient2(low="blue3",mid="white",high="red3",midpoint=0,limits=c(-1,1)
                         ,breaks=c(-1,0,1),labels=c("-1","0","1"),na.value="grey50")+
    theme_classic()+theme(axis.text.y=element_text(size=9)
                          ,axis.text.x=element_text(angle=45,vjust=1,hjust=0.9,size=7)
                          ,plot.title=element_text(size=11)
                          ,legend.title=element_text(size=7)
                          ,legend.text=element_text(size=7)
                          ,plot.subtitle=element_text(size=8,face='italic'))+
    labs(x="",y="",fill="Correlation",title=paste("Correlations across",title_suffix))

  #Add subtitle
  if (subtitle) {
    ht<-ht+labs(subtitle=paste("Top",length(unique(heat_long$protein)),"proteins"))
  }

  #Add significance to the plot
  if (annotate_signif) {
    ht<-ht+ggplot2::geom_point(data=heat_long |> dplyr::filter(is_signif)
                               ,color="black",size=1,shape=8)
  }
  return(ht)
}

#' Function to generate delta-volcano plots
#' @export
#' @examples
#' ## Placeholder Example ##
plot_delta_volcano<-function(comb_corr,datasetA,datasetB
                             ,datasetA_name,datasetB_name,highlight_n=15){
  corr_mat<-comb_corr |>
    select(dataset_label,protein,corr,qval) |>
    pivot_wider(names_from=dataset_label,values_from=c(corr,qval))
  xcor<-paste0("corr_",datasetA)
  ycor<-paste0("corr_",datasetB)
  xp<-paste0("qval_",datasetA)
  yp<-paste0("qval_",datasetB)
  df<-corr_mat |> mutate(delta_corr=.data[[xcor]] - .data[[ycor]],
                          p_change=pmax(.data[[xp]],.data[[yp]],na.rm=TRUE))
  top<-df |> arrange(desc(abs(delta_corr))) |> slice_head(n=highlight_n)
  ggplot(df,aes(delta_corr, -log10(p_change)))+
    geom_point(alpha=0.6,color="grey50")+
    geom_point(data=top,color="red",size=2)+
    geom_text_repel(data=top,aes(label=protein),size=3)+
    geom_vline(xintercept=0,linetype="dashed",color="red")+theme_classic()+
    labs(title=paste("Delta volcano:",datasetA_name,"vs",datasetB_name)
         ,x="Delta correlation",y="-log10(p)")
}

#' Function to generate correlation plot of similarity matrix datasets
#' @export
#' @examples
#' ## Placeholder Example ##
corr_sim_mtx_plot<-function(comb_corr,subsets,dataset_names,title_suffix){

  #Prepare correlation matrix
  df<-comb_corr |>
    dplyr::filter(dataset_label %in% subsets) |>
    dplyr::select(dataset_label,protein,corr) |>
    tidyr::pivot_wider(names_from=dataset_label,values_from=corr)
  colnames(df)<-c("protein",dataset_names)
  similarity_matrix<-cor(df[-1],use="pairwise.complete.obs")

  #Convert to long format
  sim_df<-as.data.frame(similarity_matrix) |>
    tibble::rownames_to_column("dataset1") |>
    tidyr::pivot_longer(-dataset1,names_to="dataset2",values_to="similarity")

  #Keep only lower triangle
  sim_df<-sim_df |>
    dplyr::mutate( i =match(dataset1,dataset_names),j =match(dataset2,dataset_names)) |>
    dplyr::filter(i > j)

  #Plot
  ggplot(sim_df,aes(dataset2,dataset1,fill=similarity))+
    geom_tile(color="white")+
    geom_text(aes(label=round(similarity,2)), size=4)+
    scale_fill_gradient2(low="#4575B4",mid="white",high="#D73027",midpoint=0.5,limits=c(0,1)) +
    coord_fixed()+theme_classic()+
    labs(title=paste("Similarity matrix of correlations between",title_suffix)
         ,x="",y="",fill="Network\nsimilarity")
}
