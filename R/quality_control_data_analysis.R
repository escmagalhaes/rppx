### QC analysis functions ###

#' Function to check for sample outliers
#' @export
#' @examples
#' ## Placeholder Example ##
qc_sample<-function(df,proteins){ #df is samples x proteins

  #Sanity check
  if(sum(!proteins %in% colnames(df))>0){
    stop("Some protein names are not in dataframe. Please re-check.\n")
  }

  ##Identify Sample outliers##
  #Plot1
  sample_dist_mat<-as.matrix(dist(df[,proteins]))
  median_dist<-apply(sample_dist_mat,1,median)
  z_med_dist<-scale(median_dist)
  z_med_df<-data.frame(values=z_med_dist[,1,drop=F])
  z_med_df<-tibble::rownames_to_column(z_med_df,'id')
  z_med_df$idx<-1:nrow(z_med_df)
  z_med_df$outliers<-factor(ifelse(abs(z_med_df$values) > 3,'strong_outlier'
                                   ,ifelse(abs(z_med_df$values) > 2,'moderate_outlier'
                                           ,'not_outlier')))
  z_med_df$labs<-factor(ifelse(z_med_df$outliers=='strong_outlier',z_med_df$id
                               ,ifelse(z_med_df$outliers=='moderate_outlier'
                                       ,z_med_df$id,'not_outlier')))
  plot_out<-ggplot(z_med_df,aes(x=idx,y=values,color=outliers,label=labs))+geom_point()+
    geom_text(data=subset(z_med_df,labs !='not_outlier'),aes(label=labs),vjust=1.5)+ggtitle("Sample Outlier Score")+
    ylim((min(z_med_df$values)-5),(max(z_med_df$values)+5))+ylab('Z-score')+
    geom_hline(yintercept=c(-3,3),linetype="dashed",color="red",linewidth=1)+
    geom_hline(yintercept=c(-2,2),linetype="dashed",color="blue",linewidth=1)+
    scale_color_manual(values=c("not_outlier"="gray","strong_outlier"="red","moderate_outlier"="blue"))+
    theme(axis.text.y=element_text(size=10)
          ,axis.text.x=element_blank(),axis.title.x=element_blank(),legend.position='none'
          ,plot.title=element_text(size=12,hjust=0,face='bold',color='black',margin=margin(t=0,r=0,b=5,l=0))
          ,panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank()
          ,axis.line=element_line(linewidth=0.25,linetype="solid",colour="black"),panel.border=element_blank()
          ,axis.ticks=element_line(linewidth=0.25),axis.ticks.length=unit(0.025,"cm"))

  strong_outliers<-z_med_df[z_med_df$outliers=='strong_outlier','id']
  moderate_outliers<-z_med_df[z_med_df$outliers=='moderate_outlier','id']

  #Plot2
  pca_raw<-prcomp(df[,proteins],center=T,scale.=T)
  pca_df<-data.frame(PC1=pca_raw$x[,1],PC2=pca_raw$x[,2])
  pca_df<-tibble::rownames_to_column(pca_df,'id')
  pca_df<-merge(pca_df,z_med_df[,c('id','labs','outliers')],by='id',all=T)
  plot_pca_out<-ggplot(pca_df,aes(x=PC1,y=PC2,color=outliers))+geom_point()+
    geom_text(data=subset(pca_df,labs !='not_outlier'),aes(label=labs),vjust=1.5)+
    ggtitle("PCA samples")+
    scale_color_manual(values=c("not_outlier"="gray","strong_outlier"="red","moderate_outlier"="blue"))+
    theme(axis.text.x=element_blank(),legend.position='none'
          ,plot.title=element_text(size=12,hjust=0,face='bold',color='black',margin=margin(t=0,r=0,b=5,l=0))
          ,panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank()
          ,axis.line=element_line(linewidth=0.25,linetype="solid",colour="black"),panel.border=element_blank()
          ,axis.ticks=element_line(linewidth=0.25),axis.ticks.length=unit(0.025,"cm"))

  return(list(
    plot1_sample=plot_out,
    plot2_sample=plot_pca_out,
    strong_out_sample=strong_outliers,
    size_strong_out_sample=length(strong_outliers),
    moderate_out_sample=moderate_outliers,
    size_moderate_out_sample=length(moderate_outliers)
  ))
}

#' Function for checking group imbalance (var=variable of choice)
#' @export
#' @examples
#' ## Placeholder Example ##
chk_gp_imbalance<-function(df_unfilt,df_filt,var,sample_rm_thr=30,sample_diff_thr=20) {

  # Accept either a column name string, vector, or factor
  if (!is.factor(df_unfilt[[var]])) df_unfilt[[var]]<-factor(df_unfilt[[var]])
  if (!is.factor(df_filt[[var]])) df_filt[[var]]<-factor(df_filt[[var]])

  # Compute group sizes
  original_counts<-table(df_unfilt[[var]])
  filtered_counts<-table(df_filt[[var]])

  # Align levels in case a group was completely removed
  all_levels<-union(names(original_counts),names(filtered_counts))
  original_counts<-original_counts[all_levels]
  filtered_counts<-filtered_counts[all_levels]
  filtered_counts[is.na(filtered_counts)]<-0

  # Compute percentage removed per group
  sample_removed<-ifelse(original_counts==0,0
                         ,(original_counts - filtered_counts) / original_counts * 100)
  names(original_counts)<-all_levels
  names(filtered_counts)<-all_levels
  names(sample_removed)<-all_levels

  # Report
  cat("Group balance check:\n")
  for (g in all_levels) {
    cat(sprintf("  %s: %d to %d (%.1f%% removed)\n",
                g, original_counts[g], filtered_counts[g], sample_removed[g]))
  }

  # Checks
  warning_flag<-FALSE
  if (any(sample_removed > sample_rm_thr)) {
    warning(sprintf("QC removed more than %d%% of samples in at least one group: %s",
                    sample_rm_thr,
                    paste(names(sample_removed),round(sample_removed,1),sep="=",collapse=", ")))
    warning_flag<-TRUE
  }

  if (max(sample_removed) - min(sample_removed) > sample_diff_thr) {
    warning(sprintf("QC removed samples disproportionately across groups (>%d%% difference): %s",
                    sample_diff_thr,
                    paste(names(sample_removed),round(sample_removed,1),sep="=",collapse=", ")))
    warning_flag<-TRUE
  }

  if (!warning_flag) {
    cat("No imbalances detected. You are good to go!\n")
  }

  #Return filtered df only if warnings are FALSE, else skip QC
  df_final<- if (warning_flag) df_unfilt else df_filt

  invisible(list(
    df_filtered=df_final,
    original_counts=as.vector(original_counts),
    filtered_counts=as.vector(filtered_counts),
    sample_removed=as.vector(sample_removed),
    warning_flag=warning_flag,
    groups=all_levels
  ))
}

#' Function for sample outlier detection
#' @export
#' @examples
#' ## Placeholder Example ##
detect_sample_outliers<-function(df,proteins) {

  sample_out<-qc_sample(df,proteins)
  strong_out<-sample_out$strong_out_sample
  moderate_out<-sample_out$moderate_out_sample

  list(strong=strong_out,moderate=moderate_out,plots=sample_out)
}

#' Function for interaction after sample outlier detection
#' @export
#' @examples
#' ## Placeholder Example ##
handle_all_sample_outliers<-function(strong_out_sample=character(0),
                                     moderate_out_sample=character(0),
                                     remove_strong_only=FALSE,
                                     remove_all=FALSE,
                                     keep_all=FALSE,
                                     allow_interactive_override=TRUE,
                                     verbose=TRUE) {

  # Combine outlier candidates
  outlier_candidates <- c(strong_out_sample, moderate_out_sample)
  if (length(outlier_candidates) == 0) {
    if (verbose) cat("\nNo outliers detected.\n")
    return(character(0))
  }

  # Default removal based on flags
  samples_to_remove <- character(0)
  if (keep_all) {
    samples_to_remove <- character(0)
  } else if (remove_all) {
    samples_to_remove <- outlier_candidates
  } else if (remove_strong_only) {
    samples_to_remove <- strong_out_sample
  }

  # Interactive override
  if (allow_interactive_override && interactive() && length(outlier_candidates) > 0) {
    cat("\nDetected sample outliers:\n")

    # Create table with sample IDs
    all_out_df <- data.frame(
      sample_id = outlier_candidates,
      type = ifelse(outlier_candidates %in% strong_out_sample, "strong", "moderate"),
      stringsAsFactors = FALSE
    )

    # Print with color hints
    for (i in seq_len(nrow(all_out_df))) {
      color <- ifelse(all_out_df$type[i] == "strong", "\033[31m", "\033[34m")
      reset <- "\033[39m"
      cat(i, ":", color, all_out_df$sample_id[i], reset, "(", all_out_df$type[i], ")\n")
    }

    cat("\nInteractive option: enter sample IDs to KEEP, space separated;\n",
        "or press ENTER to keep all;\n",
        "or write 'del' to remove all detected outliers:\n")
    input <- readline(prompt = "")
    input <- trimws(input)

    if (tolower(input) == "del") {
      samples_to_remove <- outlier_candidates
      if (verbose) cat("\nAll outliers removed.\n")
    } else if (input == "") {
      samples_to_remove <- character(0)
      if (verbose) cat("\nAll outliers retained.\n")
    } else {
      keep_ids <- strsplit(input, "\\s+")[[1]]
      keep_ids <- keep_ids[keep_ids %in% outlier_candidates]  # validate IDs
      if (length(keep_ids) == 0) {
        warning("No valid sample IDs entered. No samples removed.")
        samples_to_remove <- character(0)
      } else {
        samples_to_remove <- setdiff(outlier_candidates, keep_ids)
        if (verbose) cat("\nRemoving:", paste(samples_to_remove, collapse = ", "), "\n")
      }
    }
  } else {
    if (verbose) {
      cat("\nNon-interactive or auto-keep mode: removing", length(samples_to_remove),
          "sample(s) by default.\n")
    }
  }

  return(samples_to_remove)
}

#' Function for determining how to deal with sample outliers
#' @export
#' @examples
#' ## Placeholder Example ##
apply_sample_policy<-function(df,strong,moderate,policy,grouping_variable
                              ,sample_rm_thr,sample_diff_thr,verbose=TRUE) {

  #Decide which samples should be removed
  samples_to_remove<-switch(policy,keep_all=character(0),remove_strong=strong,remove_all=c(strong,moderate)
                            ,interactive=handle_all_sample_outliers(strong,moderate,verbose=verbose))

  #If nothing to remove, return original data
  if(length(samples_to_remove)==0){
    if(verbose) cat("No samples removed.\n")
    return(list(df_filtered=df,removed=character(0),imbalance_flag=FALSE))
  }

  #Temporary filtered data
  df_temp<-df[!rownames(df) %in% samples_to_remove,,drop=FALSE]

  #Check group imbalance only if grouping_variable is provided
  if (!is.null(grouping_variable) && grouping_variable %in% colnames(df)) {
    imbalance_res<-chk_gp_imbalance(df_unfilt=df,df_filt=df_temp,var=grouping_variable,
                                    sample_rm_thr=sample_rm_thr,sample_diff_thr=sample_diff_thr)
    imbalance_flag<-imbalance_res$warning_flag
    df_final<-imbalance_res$df_filtered
  } else {
    if (verbose && is.null(grouping_variable)) cat("No grouping variable supplied. Skipping imbalance check.\n")
    imbalance_flag<-FALSE
    df_final<-df_temp
  }
  if(imbalance_res$warning_flag){

    if(verbose) { warning("Group imbalance detected. Skipping sample removal.") }

    return(list(df_filtered=df,removed=character(0),imbalance_flag=TRUE))
  }

  #Return filtered data
  list(df_filtered=imbalance_res$df_filtered,removed=samples_to_remove,imbalance_flag=FALSE)
}

#' Function to check for protein (feature) outliers
#' @export
#' @examples
#' ## Placeholder Example ##
qc_feature<-function(df_sample_filtered,proteins,check_cor=TRUE,check_mad=FALSE
                     ,cor_cutoff=0.95,ptn_z_thrs=3,freq_thrs=0.10){

  #Sanity check
  if (sum(!proteins %in% colnames(df_sample_filtered)) > 0) {
    stop("Some protein names are not in dataframe. Please re-check.\n")
  }

  #Define inputs
  df<-df_sample_filtered[,proteins]
  n_samples<-nrow(df)

  #Scale dataset (column-wise)
  df_scaled<-scale(df,center=T,scale=T)

  ##Identify and delete near-zero-variance outliers##
  nzv<-caret::nearZeroVar(df_scaled,saveMetrics=TRUE)
  if (any(nzv$nzv)){
    nzv_ptn<-rownames(nzv[nzv$nzv==TRUE,])
    df_nzv_filt<-df[,!colnames(df) %in% nzv_ptn]
    cat("Near-zero variance proteins detected.\n")
  } else {
    nzv_ptn<-character(0)
    df_nzv_filt<-df
    cat("No near-zero variance proteins detected.\n")
  }

  #(optional) Flag and delete near-duplicates by correlation
  high_cor_all<-character(0)
  high_cor_remove<-character(0)
  high_cor_keep<-character(0)
  if(isTRUE(check_cor)){
    df_nzv_filt_scaled<-scale(df_nzv_filt,center=T,scale=T)
    df_cor_full<-cor(df_nzv_filt_scaled,use="pairwise.complete.obs")
    diag(df_cor_full)<-NA  #exclude self-correlation for mean absolute correlation computation downstream
    df_cor<-df_cor_full
    df_cor[lower.tri(df_cor,diag=FALSE)]<-NA #turn lower triangle of the matrix into NA to extract high_cor without duplicates for pair detection
    high_cor<-which(abs(df_cor) > cor_cutoff,arr.ind=TRUE)
    #Identify pairs that are extremely correlated
    if (nrow(high_cor) > 0) {
      high_cor_pairs<-data.frame(
        ptn1=rownames(df_cor)[high_cor[,1]],
        ptn2=colnames(df_cor)[high_cor[,2]],
        corr=df_cor[high_cor]
      )
      #Calculate mean absolute correlation and remove the proteins that have lower value
      mean_abs_cor<-rowMeans(abs(df_cor_full),na.rm=TRUE)
      remove_ptn<-apply(high_cor_pairs,1, function(pair) {
        ptn1<-pair["ptn1"]
        ptn2<-pair["ptn2"]
        if (mean_abs_cor[ptn1] >= mean_abs_cor[ptn2]) ptn1 else ptn2
      })
      high_cor_remove<-unique(remove_ptn)
      high_cor_keep<-unique(c(high_cor_pairs$ptn1,high_cor_pairs$ptn2))
      high_cor_keep<-high_cor_keep[!high_cor_keep %in% high_cor_remove]
      high_cor_all<-unique(c(high_cor_remove,high_cor_keep))
      df_nzv_cor_filt<-df_nzv_filt[,!colnames(df_nzv_filt) %in% high_cor_remove]
      cat("Near-duplicate proteins detected.\n")

    } else {
      high_cor_all<-character(0)
      high_cor_remove<-character(0)
      high_cor_keep<-character(0)
      df_nzv_cor_filt<-df_nzv_filt
      cat("No near-duplicate proteins detected.\n")
    }
  } else {
    df_nzv_cor_filt<-df_nzv_filt
    high_cor_all<-character(0)
    high_cor_remove<-character(0)
    high_cor_keep<-character(0)
  }

  #Check for outliers in the remaining data using Mean absolute deviation(MAD)-based SD
  #This is better than standard SD because it adressess skewness better
  #But be careful: if MAD=0, then scores will be Inf
  #Use conditional statement to avoid MAD=Inf
  out_ptn<-character(0)
  out_plot<-NULL
  if(isTRUE(check_mad)){
    robust_scale<-function(x) {
      if(all(is.na(x))) return(rep(NA,length(x)))
      m<-mad(x,constant=1.4826,na.rm=TRUE)
      if (m==0) {
        s<-sd(x,na.rm=TRUE)
        if (s==0) return(rep(0,length(x)))  #NZV will caught this, so it is just safeguard
        return((x - mean(x,na.rm=TRUE)) / s)
      }
      (x - median(x,na.rm=TRUE)) / m
    }
    df_scaled_robust<-apply(df_nzv_cor_filt,2,robust_scale)

    #Check for outliers based on 3 SD
    ptn_outliers<-apply(df_scaled_robust,2, function(x) mean(abs(x) > ptn_z_thrs,na.rm=TRUE))
    out_ptn<-names(ptn_outliers[ptn_outliers > max(freq_thrs, 2/n_samples)])
    out_plot<-ggplot(data.frame(value=ptn_outliers),aes(x=value))+geom_vline(xintercept=freq_thrs,color="red",linetype="dashed")+
      geom_histogram(binwidth=0.005,fill="gray40")+ggtitle("Histogram of proportion of samples with outliers (>3 MAD-SD) per protein")+
      theme(axis.text.x=element_text(size=8,face="bold",color='black')
            ,axis.title.x=element_blank(),axis.title.y=element_blank(),legend.position='none'
            ,plot.title=element_text(size=10,hjust=0.5,face='bold',color='black',margin=margin(t=0,r=0,b=5,l=0))
            ,panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank()
            ,axis.line=element_line(linewidth=0.25,linetype="solid",colour="black"),panel.border=element_blank()
            ,axis.ticks=element_line(linewidth=0.25),axis.ticks.length=unit(0.025,"cm"))#+coord_flip()
    if (length(out_ptn)==0){
      cat("No protein outliers detected by SD threshold.\n")
    }
  } else {
    out_ptn<-character(0)
  }

  #Create object with all feature outliers
  feature_outliers<-unique(c(nzv_ptn,high_cor_remove,out_ptn))

  return(list(
    df_after_nzv_filt=df_nzv_filt,
    near_zero_var_proteins=nzv_ptn,
    size_near_zero_var_proteins=length(nzv_ptn),
    df_after_cor_filt=df_nzv_cor_filt,
    high_cor_proteins=high_cor_all,
    size_high_cor_proteins=length(high_cor_all),
    high_cor_removed=high_cor_remove,
    size_high_cor_removed=length(high_cor_remove),
    high_cor_kept=high_cor_keep,
    size_high_cor_kept=length(high_cor_keep),
    outliers_proteins=out_ptn,
    size_outliers_proteins=length(out_ptn),
    feature_out=feature_outliers,
    size_feature_out=length(feature_outliers),
    histogram_outliers=out_plot

  ))
}

#' Function to detect feature outliers
#' @export
#' @examples
#' ## Placeholder Example ##
detect_feature_outliers<-function(df,proteins,check_cor,check_mad,cor_cutoff,ptn_z_thrs
                                  ,freq_thrs,feature_outlier_policy="keep_all") {

  feature_out<-qc_feature(df_sample_filtered=df,proteins=proteins,check_cor=check_cor,
                          check_mad=check_mad,cor_cutoff=cor_cutoff,ptn_z_thrs=ptn_z_thrs,freq_thrs=freq_thrs)

  mad_outliers<-character(0)
  if (check_mad) {

    print(feature_out$histogram_outliers) #always print plot

    if(feature_outlier_policy=="interactive"){

      cat("\nCurrent MAD frequency threshold:",freq_thrs,"\n")

      response<-trimws(readline("Press ENTER to accept threshold or type a new value (0-1): "))

      if(response !=""){
        new_thr<-as.numeric(response)

        if(!is.na(new_thr) && new_thr >= 0 && new_thr <= 1){
          freq_thrs<-new_thr
        } else {
          warning("Invalid threshold. Using default.")
        }
      }

      cat("Using MAD threshold:",freq_thrs,"\n")
    }

    robust_scale<-function( x ) {
      if(all(is.na( x ))) return(rep(NA,length( x )))
      m<-mad( x, constant=1.4826,na.rm =TRUE)
      if (m==0) {
        s<-sd( x, na.rm=TRUE)
        if (s==0) return(rep(0,length(x)))
        return(( x - mean( x, na.rm=TRUE)) / s)
      }
      (x - median(x, na.rm = TRUE)) / m
    }
    df_scaled<-apply(feature_out$df_after_cor_filt,2,robust_scale)
    df_scaled<-as.matrix(df_scaled)
    ptn_outliers<-apply(df_scaled,2,function( x ) mean(abs(x) > ptn_z_thrs,na.rm=TRUE))
    mad_outliers<-names(ptn_outliers[ptn_outliers > freq_thrs])
  }
  list(nzv=feature_out$near_zero_var_proteins,cor=feature_out$high_cor_removed,
       mad=mad_outliers,raw=feature_out)
}

#' Function for determining how to deal with feature outliers
#' @export
#' @examples
#' ## Placeholder Example ##
apply_feature_policy<-function(df,proteins,nzv,cor,mad,policy,freq_thrs,verbose=TRUE) {

  #decide what to do with MAD outliers
  mad_remove<-switch(policy,keep_all=character(0),remove=mad,interactive=mad)

  #proteins to remove (NZV + correlation always removed)
  proteins_to_remove<-unique(c(nzv,cor,mad_remove))

  if (length(proteins_to_remove)==0) {
    if (verbose) cat("No proteins removed.\n")
    return(list(df=df,proteins=proteins,removed=character(0)))
  }

  #filter dataset
  df_final<-df[,!colnames(df) %in% proteins_to_remove,drop=FALSE]

  #update protein vector
  ptn_final<-setdiff(proteins,proteins_to_remove)

  list(df=df_final,proteins=ptn_final,removed=proteins_to_remove)
}

#' Function to get QC summary
#' @export
#' @examples
#' ## Placeholder Example ##
build_qc_report<-function(n_samples_original,n_samples_final,n_proteins_original
                          ,n_proteins_final,strong_outliers,moderate_outliers
                          ,samples_removed,proteins_removed,freq_threshold
                          ,imbalance_flag,sample_plots=NULL,feature_plots=NULL) {

  report <- list(

    summary = data.frame(
      metric = c(
        "Samples (original)",
        "Samples (final)",
        "Proteins (original)",
        "Proteins (final)",
        "Strong sample outliers",
        "Moderate sample outliers",
        "Samples removed",
        "Proteins removed",
        "MAD frequency threshold",
        "Group imbalance detected"
      ),

      value = c(
        n_samples_original,
        n_samples_final,
        n_proteins_original,
        n_proteins_final,
        length(strong_outliers),
        length(moderate_outliers),
        length(samples_removed),
        length(proteins_removed),
        freq_threshold,
        imbalance_flag
      ),

      stringsAsFactors = FALSE
    ),

    decisions = list(
      strong_outliers = strong_outliers,
      moderate_outliers = moderate_outliers,
      samples_removed = samples_removed,
      proteins_removed = proteins_removed,
      freq_threshold = freq_threshold,
      imbalance_flag = imbalance_flag
    ),

    plots = list(
      sample_qc = sample_plots,
      feature_qc = feature_plots
    )
  )

  class(report) <- "qc_report"
  return(report)
}
