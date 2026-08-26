### Differential Expression Analysis with limma

#' Function to compute DE analysis with limma with two variables
#' @export
#' @examples
#' ## Placeholder Example ##
binary_de_model<-function(df,grouping_var,proteins,covars=NULL
                          ,pval_cutoff=0.05,LFC_cutoff=0.25
                          ,treat=FALSE,ref_level=NULL,verbose=TRUE){

  #Make sure df is a dataframe
  df<-as.data.frame(df)
  rownames(df)<-seq_len(nrow(df))

  #Grouping variable (should be a factor with 2 levels)
  if (!is.null(ref_level)) {
    group<-relevel(factor(df[[grouping_var]]),ref_level)
  } else {
    #default=biggest group as reference
    ref_level<-names(table(df[[grouping_var]]))[which.max(table(df[[grouping_var]]))]
    group<-relevel(factor(df[[grouping_var]]),ref_level)
  }
  names(group)<-rownames(df)

  #Safety check
  if (nlevels(group) != 2) {
    stop("\ngroup must have exactly 2 levels.\n")
  }
  if (verbose) {
    message("Reference level set to: ",levels(group)[1])
    message("Comparing ",levels(group)[2]," vs ",levels(group)[1])
  }

  #Adjust dataset for limma (should be a initial df of samples x variables)
  #Safety check
  if(!all(proteins %in% colnames(df)))
    stop("Some proteins not found in dataframe.")

  expr<-t(df[,proteins]) #input here the protein names to analyse

  #Generate metadata df from main dataframe
  meta_df<-data.frame(.grp_var=group)

  if (!is.null(covars)) {
    if (!all(covars %in% colnames(df))) {
      stop("\nSome covariates not in main dataframe.\n")
    }
    meta_df<-cbind(meta_df,df[,covars,drop=FALSE])

    #Adjust covars
    meta_df[covars]<-lapply(meta_df[covars], function(x) {
      if (is.character(x)) factor(x) else x
    })
  }

  #Generate formula for model matrix
  if (!is.null(covars)) {
    formula_mod<-paste("~ .grp_var +",paste(covars,collapse=" + "))
  } else {
    formula_mod<-"~ .grp_var"
  }

  #Design matrix
  design<-model.matrix(as.formula(formula_mod),data=meta_df)
  colnames(design)<-make.names(colnames(design))

  #Safety checks
  if (!all(colnames(expr)==rownames(design))) {
    stop("\nSample order mismatch between expression matrix and design matrix. Please, revise.\n")
  }

  if (qr(design)$rank < ncol(design)) {
    stop("\nDesign matrix is not full rank. Possible collinearity between covariates.\n")
  }

  #Fit model
  fit<-limma::lmFit(expr,design)

  #Define workflow based on optional treat (default=F)
  #Safety check for group variable that matters for extracting coefs
  group_coef<-colnames(design)[grep("^\\.grp_var",colnames(design))]
  if (length(group_coef) != 1) {
    stop("\nCould not uniquely identify group coefficient. Revise variable names for multiple variables with 'group' in the name\n")
  }

  if(treat){
    fit<-limma::treat(fit,lfc=LFC_cutoff)
    results<-limma::topTreat(fit,coef=group_coef,number=Inf,sort.by="P")
    results$is_signif<-results$adj.P.Val < pval_cutoff
    results$class_high<-factor(ifelse(results$is_signif==TRUE & results$logFC > 0,'up'
                                      ,ifelse(results$is_signif==TRUE & results$logFC < 0,'down','ns'))
                               ,levels=c('up','down','ns'))
    results$class_low<-factor(ifelse(results$class_high=='up','down'
                                     ,ifelse(results$class_high=='down','up','ns'))
                              ,levels=c('up','down','ns'))
    results$abs_diff_high<-abs(results$logFC)
    results$abs_diff_low<-abs(results$logFC)
    signif_ptn<-rownames(results[results$is_signif,])
    up_ptn<-rownames(results[results$class_high=="up",])
    down_ptn<-rownames(results[results$class_high=="down",])
  } else {
    fit<-limma::eBayes(fit,robust=TRUE)
    results<-limma::topTable(fit,coef=group_coef,number=Inf,sort.by="P")
    results$is_signif<-results$adj.P.Val < pval_cutoff & abs(results$logFC) > LFC_cutoff
    results$class_high<-factor(ifelse(results$is_signif==TRUE & results$logFC > LFC_cutoff,'up'
                                      ,ifelse(results$is_signif==TRUE & results$logFC < -LFC_cutoff,'down','ns'))
                               ,levels=c('up','down','ns'))
    results$class_low<-factor(ifelse(results$class_high=='up','down'
                                     ,ifelse(results$class_high=='down','up','ns'))
                              ,levels=c('up','down','ns'))
    results$abs_diff_high<-abs(results$logFC)
    results$abs_diff_low<-abs(results$logFC)
    signif_ptn<-rownames(results[results$is_signif,])
    up_ptn<-rownames(results[results$class_high=="up",])
    down_ptn<-rownames(results[results$class_high=="down",])
  }

  ##Generate volcano plot
  plot_vol<-plot_vol_de(df=results,LFC='logFC',pval='adj.P.Val',groups=levels(group)
                        ,down_ptn=down_ptn,up_ptn=up_ptn,pval_cutoff=pval_cutoff
                        ,LFC_cutoff=LFC_cutoff)

  #Calculate mean expression of each ALL protein by group and attach to results (matching by rownames)
  mean_level1<-colMeans(df[group==levels(group)[1],proteins,drop=FALSE],na.rm=TRUE)
  mean_level2<-colMeans(df[group==levels(group)[2],proteins,drop=FALSE],na.rm=TRUE)
  results[[paste0("mean_expr_", tolower(levels(group)[1]))]]<-mean_level1[rownames(results)]
  results[[paste0("mean_expr_", tolower(levels(group)[2]))]]<-mean_level2[rownames(results)]

  if (verbose) message(length(signif_ptn)," significant proteins found.")

  return(list(
    fit=fit,
    design=design,
    results=results,
    vplot=plot_vol,
    all_sig_ptn=signif_ptn,
    sig_ptn_size=length(signif_ptn),
    sig_ptn_class=list(all=signif_ptn,up=up_ptn,down=down_ptn)
  ))
}

#' Function to compute DE analysis with limma with three or more variables (ordered model)
#' @export
#' @examples
#' ## Placeholder Example ##
ordered_de_model<-function(df,grouping_var,proteins,covars=NULL
                           ,pval_cutoff=0.05,LFC_cutoff=0.25,verbose=TRUE){

  #Make sure df is a dataframe
  df<-as.data.frame(df)
  rownames(df)<-seq_len(nrow(df))

  #Treat grouping_var as numeric
  group<-as.numeric(df[[grouping_var]])
  if (verbose) message("Running ordered limma model on: ",grouping_var)

  #Expression matrix
  if (!all(proteins %in% colnames(df)))
    stop("Some proteins not found in dataframe.")
  expr<-t(df[,proteins,drop=FALSE])

  #Metadata
  meta_df<-data.frame(.grp_var=group)
  if (!is.null(covars)) {
    if (!all(covars %in% colnames(df)))
      stop("Some covariates not in main dataframe.")
    meta_df<-cbind(meta_df,df[,covars,drop=FALSE])
    meta_df[covars]<-lapply(meta_df[covars],function(x) {
      if (is.character(x)) factor(x) else x
    })
  }

  #Formula and design
  formula_mod<-if (!is.null(covars)) {
    paste("~ .grp_var +",paste(covars,collapse=" + "))
  } else "~ .grp_var"

  design<-model.matrix(as.formula(formula_mod),data=meta_df)
  colnames(design)<-make.names(colnames(design))

  #Safety checks
  if (!all(colnames(expr)==rownames(design)))
    stop("Sample order mismatch.")
  if (qr(design)$rank < ncol(design))
    stop("Design matrix is not full rank.")

  #Model Fit
  fit<-limma::lmFit(expr,design)
  fit<-limma::eBayes(fit,robust=TRUE)
  group_coef<-colnames(design)[grep("^\\.grp_var", colnames(design))]
  results<-limma::topTable(fit,coef=group_coef,number=Inf,sort.by="P")
  results$is_signif<-results$adj.P.Val < pval_cutoff & abs(results$logFC) > LFC_cutoff
  results$direction<-ifelse(results$is_signif & results$logFC > 0,"up",
                            ifelse(results$is_signif & results$logFC < 0,"down","ns"))

  signif_ptn<-rownames(results[results$is_signif,])
  up_ptn<-rownames(results[results$direction=="up",])
  down_ptn<-rownames(results[results$direction=="down",])

  if (verbose) message(length(signif_ptn)," significant proteins found.")

  return(list(
    fit           = fit,
    design        = design,
    results       = results,
    vplot         = NULL,
    all_sig_ptn   = signif_ptn,
    sig_ptn_size  = length(signif_ptn),
    sig_ptn_class = list(all=signif_ptn,up=up_ptn,down=down_ptn)
  ))

}

#' Function to generate Volcano plot from DE analysis with limma
#' @export
#' @examples
#' ## Placeholder Example ##
plot_vol_de<-function(df,LFC,pval,groups,down_ptn,up_ptn,pval_cutoff=0.05
                      ,LFC_cutoff=0.25){

  if(length(groups) != 2)
    stop("groups must contain exactly two group names")

  if(!LFC %in% colnames(df))
    stop("LFC column not found in dataframe.")

  if(!pval %in% colnames(df))
    stop("pval column not found in dataframe.")

  if(!all(c(up_ptn,down_ptn) %in% rownames(df)))
    warning("Some proteins in up_ptn/down_ptn not present in dataframe.")

  #Create key color variable
  down_res<-rownames(df) %in% down_ptn
  up_res<-rownames(df) %in% up_ptn

  if(any(up_res & down_res))
    stop("A protein cannot be both up- and down-regulated at the same time.")

  key_col<-ifelse(down_res,'deepskyblue2',ifelse(up_res,'firebrick1','gray50'))
  names(key_col)[key_col=='firebrick1']<-'Upregulated'
  names(key_col)[key_col=='gray50']<-'NS'
  names(key_col)[key_col=='deepskyblue2']<-'Downregulated'

  #Adjust ylim of plot
  yvals<-(-log10(pmax(df[[pval]],0.000001)))
  ylim_max<-ceiling(max(yvals[is.finite(yvals)]))

  #Adjust xlim of the plot
  x_limits<-c(min(df[[LFC]],na.rm=TRUE) - 0.5,max(df[[LFC]],na.rm=TRUE) + 0.5)

  vplot<-EnhancedVolcano(df,lab=row.names(df),cutoffLineType='twodash'
                         ,cutoffLineWidth=0.25,widthConnectors=0.15,pointSize=1.5,labSize=3
                         ,colAlpha=1,colCustom=key_col,caption=NULL,drawConnectors=T
                         ,legendLabSize=10,legendIconSize=3,FCcutoff=LFC_cutoff
                         ,pCutoff=pval_cutoff,subtitle=NULL,arrowheads=F,border='partial'
                         ,borderColour='black',borderWidth=0.25
                         ,title=paste('Up and Down-regulated Proteins in',groups[2],'vs',groups[1])
                         ,x=LFC,y=pval,boxedLabels=T,xlim=x_limits#,ylim=c(0,ylim_max+0.1)
                         ,axisLabSize=10,titleLabSize=10,gridlines.major=F,gridlines.minor=F
  )+theme(axis.ticks=element_line(linewidth=0.25),axis.ticks.length=unit(0.1,"cm")
          ,legend.position="bottom",legend.margin=margin(t=0,r=0,b=0.25,l=0,unit="cm")
          ,plot.margin=margin(t=0.25,r=0.25,b=0,l=0.25,unit="cm")
          ,legend.box.margin=margin(t=-0.8,r=-1,b=0.25,l=-1,unit="cm"))#+coord_flip()
  return(vplot)
}
