### Quantile splitting of data ###

#' Function to get several quantile groups for many proteins
#' @export
#' @examples
#' ## Placeholder Example ##
get_quantile_groups<-function(data,features,n_quantiles=NULL,data_normal=NULL,
                              mode=c("quantile","normal_based","both"),
                              normal_labels=c("Below Normal","Lower Normal"
                                              ,"Upper Normal","Above Normal")){
  mode<-match.arg(mode)

  #Safety checks
  if (mode %in% c("quantile","both")) {
    if (is.null(n_quantiles))
      stop("n_quantiles must be provided for mode '",mode,"'")

    if (!all(n_quantiles==as.integer(n_quantiles)) || any(n_quantiles < 2)) {
      stop("n_quantiles must be integers >= 2")
    }

    if (any(duplicated(n_quantiles))) {
      stop("n_quantiles contains duplicates: ",
           paste(n_quantiles[duplicated(n_quantiles)],collapse=", "))
    }
  }

  if (mode %in% c("normal_based","both")) {
    if (is.null(data_normal))
      stop(paste0("data_normal must be provided for '",mode,"' mode"))

    if (length(normal_labels) != 4) {
      stop("normal_labels must have exactly 4 labels: normal_based mode always produces 4 intervals")
    }

    missing_normal<-setdiff(features,names(data_normal))
    if (length(missing_normal) > 0)
      stop("These features are missing from data_normal: ",
           paste(missing_normal,collapse=", "))
  }

  if (!all(features %in% names(data))) {
    stop("These features are missing from datatable: ",
         paste(setdiff(features,names(data)),collapse=", "))
  }

  existing_gp_cols<-grep(paste0("^(",paste(features,collapse="|"),")_gp"),names(data),value=TRUE)
  if (length(existing_gp_cols) > 0) {
    warning("These columns will be overwritten: ",paste(existing_gp_cols,collapse=", "))
  }

  #Helper make bins based on features
  make_quantile_bins<-function(data,feature,n_quantiles=NULL,data_normal=NULL,
                               mode=c("quantile","normal_based","both")) {
    mode<-match.arg(mode)
    bins<-list()
    if (mode %in% c("quantile","both")) {
      qt<-unname(stats::quantile(data[[feature]],probs=seq(0,1,length.out=n_quantiles + 1),na.rm=TRUE))
      qt<-round(qt,digits=4)
      qt[1]<- -Inf           #adjust the lower bound
      qt[length(qt)]<-Inf    #adjust the upper bound
      bins$quantile<-qt
    }

    if (mode %in% c("normal_based","both")) {
      bins$normal<-c(
        -Inf,  #adjust the lower bound
        round(min(data_normal[[feature]]),digits=4),
        round(stats::median(data_normal[[feature]]),digits=4),
        round(max(data_normal[[feature]]),digits=4),
        Inf    #adjust the upper bound
      )

      if (is.unsorted(bins$normal)) {
        warning("For '",feature,"':
                either max or min expression from data did not
                surpass or are very close to data_normal levels.")
      }
    }
    return(bins)
  }

  #Helper to cut data based on bins
  cut_feature<-function(data,feature,bins) {
    cut(data[[feature]],breaks=bins,include.lowest=TRUE)
  }

  #Helper to adjust labels
  format_cut_labels<-function(x) {
    levels(x)<-gsub(",","  ,  ",levels(x))
    return(x)
  }

  #Loop through features
  for (feature in features) {

    #Quantile-based columns
    if (mode %in% c("quantile","both")) {
      for ( q in n_quantiles) {
        bins_q<-make_quantile_bins(data,feature,n_quantiles= q ,
                                   data_normal=data_normal,mode="quantile")
        col_name_q<-paste0(feature,"_gp", q )
        data[[col_name_q]]<-format_cut_labels(cut_feature(data,feature,bins_q$quantile))
      }
    }

    #Normal-based column (single binning)
    if (mode %in% c("normal_based","both")) {
      bins_n<-make_quantile_bins(data,feature,data_normal=data_normal,mode="normal_based")
      col_name_n<-paste0(feature,"_gp_n")
      data[[col_name_n]]<-cut_feature(data,feature,bins_n$normal)
      levels(data[[col_name_n]])<-normal_labels
    }
  }
  return(data)
}

#' Function to collapse quantile groups into other groups (e.g. 'high' vs 'low')
#' @export
#' @examples
#' ## Placeholder Example ##
collapse_quantile_groups<-function(data,feature,groups,new_col=NULL,verbose=TRUE) {

  #Safety checks
  if (!is.list(groups) || is.null(names(groups))) {
    stop("'groups' must be a named list with levels (e.g. list(low=c(1,2,3),high=c(4,5)))")
  }
  if (!feature %in% names(data)) {
    stop("Variable '",feature,"' not found in data")
  }
  if (!is.factor(data[[feature]])) {
    stop("Variable '",feature,"' must be a factor. Use get_quantile_groups() first")
  }

  #Get number of levels/quantiles
  n_levels<-nlevels(data[[feature]])
  all_indices<-unlist(groups)

  #More Safety checks
  if (!all(all_indices %in% seq_len(n_levels))) {
    stop("groups contains invalid indices. Variable '",feature,
         "' has ",n_levels," levels")
  }

  if (any(duplicated(all_indices))) {
    stop("'groups' contains duplicated indices: ",
         paste(all_indices[duplicated(all_indices)],collapse=", "))
  }

  unassigned_indices<-setdiff(seq_len(n_levels),all_indices)
  if (length(unassigned_indices) > 0) {
    warning("These level indices are unassigned and will become NA: ",
            paste(unassigned_indices,collapse=", "))
  }

  #Recode variables
  level_names<-levels(data[[feature]])
  new_levels<-names(groups)
  recoded<-rep(NA_character_,nrow(data))

  for ( gp in new_levels) {
    recoded[as.integer(data[[feature]]) %in% groups[[ gp ]]]<- gp
  }

  if (is.null(new_col)) new_col<-paste0(feature,"_gp_n",length(groups))
  data[[new_col]]<-factor(recoded,levels=new_levels)

  #Inform which levels were grouped
  if(verbose) cat("\nVariable",feature,"changed to",new_col,"according to the following:\n")
  for ( gp in new_levels) {
    if(verbose) cat(paste(level_names[groups[[ gp ]]],collapse=" & ")
                    ,"recoded to", gp ,"\n")
  }
  return(data)
}

#' Function to automate best quantile split determination
#' @export
#' @examples
#' ## Placeholder Example ##
determine_optimal_quantile<-function(data,feature,time_var,event_var,
                                     q_range=2:6,padj_method="BH") {

  #Safety check
  padjust_methods<-c("holm","hochberg","hommel","bonferroni","BH","BY","fdr","none")
  if (!padj_method %in% padjust_methods) {
    stop("Invalid p-value adjustment method. Choose one of: ",
         paste(padjust_methods,collapse=", "))
  }

  results<-list()
  for ( q in q_range ) {

    #Create quantiles
    data_qt<-get_quantile_groups(
      data        = data,
      features    = feature,
      n_quantiles = q ,
      mode        = "quantile"
    )
    var<-paste0(feature,"_gp",q)

    for ( split in 1:( q - 1)) {

      #Collapse quantiles
      data_qt_binary<-collapse_quantile_groups(
        data    = data_qt,
        feature = var,
        groups  = list(Low = 1: split , High = (split + 1): q ),
        new_col = "group",
        verbose = FALSE
      )

      #Make sure quantile groups have enough cases, otherwise stop
      grp_counts<-table(data_qt_binary$group)
      if (length(grp_counts) < 2 || any(grp_counts==0)) { next }

      fit<-survival::survdiff(survival::Surv(
        data_qt_binary[[time_var]],data_qt_binary[[event_var]])~group,
        data=data_qt_binary)

      p_value<- 1 - stats::pchisq(fit$chisq,df=1)

      results[[paste( q , split , sep="_")]]<-data.frame(
        n_quantiles   =  q ,
        split_point   = split,
        low_group     = paste0(split, "/", q ),
        high_group    = paste0( q - split, "/", q ),
        n_low_group   = sum(data_qt_binary$group=="Low"),
        n_high_group  = sum(data_qt_binary$group=="High"),
        p_value       = round(p_value,4)
      )
    }
  }
  results_data<-dplyr::bind_rows(results) |>
    dplyr::mutate(q_value = stats::p.adjust(p_value,method = padj_method)) |>
    dplyr::arrange(q_value,p_value,n_quantiles)

  #Safety check
  if (nrow(results_data) == 0) {
    stop("No valid quantile split could be generated.",
         "Feature may contain too few unique values.")
  }

  #Apply best split to the actual data
  best<-results_data[1,]
  cat("\nBest split: ",best$low_group,
      " Low vs ",best$high_group," High",
      " (p = ",signif(best$p_value,3),
      ", q = ",signif(best$q_value,3),")\n")

  best_var<-paste0(feature,"_gp",best$n_quantiles)

  #Call functions one last time to attach best quantile to the data
  data_final<-get_quantile_groups(
    data        = data,
    features    = feature,
    n_quantiles = best$n_quantiles,
    mode        = "quantile"
  ) |>
    collapse_quantile_groups(
      feature = best_var,
      groups  = list(Low  = 1:best$split_point,
                     High = (best$split_point + 1):best$n_quantiles),
      new_col = paste0(feature,"_group"),
      verbose = FALSE
    )

  attr(data_final,"optimal_split_results")<-results_data
  return(data_final)
}
