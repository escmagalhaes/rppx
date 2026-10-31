### Quantile splitting of samples ###

#' Transform continuous features into categorical groups via quantile-based or
#' reference-based splitting
#'
#' For each feature, creates one or more grouping columns in the input dataframe
#' by cutting values into quantile-based intervals or reference-based cutpoints,
#' which uses the values from normal samples (controls) to create the intervals.
#' This is the primary tool for continuous-to-categorical transformation in the
#' rppx subsetting workflow. Designed for protein expression data but works on
#' any numeric feature.
#'
#' @param data A data frame or tibble containing the features to split.
#' @param features A character vector of column names to split into groups.
#'   All must be present in `data`.
#' @param n_quantiles An integer vector of quantile counts (e.g. `c(2, 3, 4)`).
#'   Each value `q` produces a column `featurename_gpq` with `q` groups.
#'   Required when `mode` is `"quantile"` or `"both"`.
#' @param data_normal A data frame of normal samples used to define interval
#'   boundaries for `mode = "normal_based"` or `"both"`. Must contain the
#'   same feature columns as `data`. Required when `mode` is `"normal_based"`
#'   or `"both"`.
#' @param mode Character. One of `"quantile"` (default), `"normal_based"`,
#'   or `"both"`. Controls which grouping strategy is applied:
#' \describe{
#'   \item{`"quantile"`}{Splits each feature into `n_quantiles` equal-frequency
#'     groups based on the distribution in `data`.}
#'   \item{`"normal_based"`}{Splits into 4 groups using the min, median, and
#'     max of the normal sample distribution as boundaries.}
#'   \item{`"both"`}{Applies both strategies, producing multiple new columns
#'     per feature.}
#' }
#' @param normal_labels A character vector of exactly 4 labels for the
#'   normal-based groups, ordered from lowest to highest. Default is
#'   `c("Below Normal", "Lower Normal", "Upper Normal", "Above Normal")`.
#'   Only used when `mode` is `"normal_based"` or `"both"`.
#'
#' @return The input `data` with new grouping columns appended:
#' \describe{
#'   \item{`featurename_gpQ`}{One column per value of `n_quantiles`, where
#'     `Q` is the number of quantile groups, present when `mode` is `"both"` or
#'     `"quantile"`. Returns a factor with levels ordered from lowest to highest
#'      values and interval labels.}
#'   \item{`featurename_gp_n`}{Normal-based grouping column, present when
#'     `mode` is `"normal_based"` or `"both"`. Returns a factor with
#'     with `normal_labels` as levels, ordered from lowest to highest.}
#' }
#'
#' @details
#' For quantile-based splitting, boundaries are computed from `data` itself
#' using equal-frequency breaks. The lower and upper bounds are set to
#' `-Inf` and `Inf` respectively to capture all values including extremes.
#'
#' For normal-based splitting, boundaries are derived from the min, median,
#' and max of `data_normal`, producing 4 biologically meaningful intervals
#' relative to the normal sample range. A warning is issued if the normal
#' range boundaries are not monotonically increasing for a given feature,
#' which may indicate that disease expression does not exceed normal levels.
#'
#' Existing grouping columns matching the pattern `featurename_gpX` will be
#' overwritten with a warning.
#'
#' @examples
#' set.seed(42)
#' my_data<-data.frame(
#'   id   = paste0("S",1:20),
#'   EGFR = rnorm(20,mean=2),
#'   AKT1 = rnorm(20,mean=1)
#' )
#' #Quantile split into 2 and 3 groups
#' get_quantile_groups(
#'   data       = my_data,
#'   features   = c("EGFR","AKT1"),
#'   n_quantiles = c(2,3),
#'   mode       = "quantile"
#' )
#'
#' @seealso [collapse_quantile_groups()], [determine_optimal_quantile()]
#'
#' @export
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

      #Safety check
      if (anyDuplicated(qt)) {
        stop("Feature '",feature,"' cannot be split into ",n_quantiles,
             " groups: because of too many tied values.",call.=FALSE)
      }

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

#' Collapse quantile groups into fewer named groups
#'
#' Recodes a quantile-grouped factor column into a smaller set of named
#' groups by mapping quantile level indices to new labels. Typically used
#' after [get_quantile_groups()], for instance, to collapse 4 quantiles into
#' "low" vs "high" grouping.
#'
#' @param data A data frame or tibble containing the factor column to recode.
#' @param feature A string. Name of the factor column to collapse. Must be a
#'   factor produced by [get_quantile_groups()].
#' @param groups A named list mapping new group labels to quantile level
#'   indices. Indices are ordered from lowest to highest expression
#'   (e.g. `list(low = c(1, 2), high = c(3, 4))`). All indices must be
#'   valid and non-duplicated. Unassigned indices become `NA` with a warning.
#' @param new_col A string. Name of the new variable to add to `data`.
#'   If `NULL` (default), the new variable name is auto-generated as
#'   `feature_gp_nK`,where `K` is the number of groups that will be made
#'   according to `groups` (e.g. `"EGFR_gp4_gp_n2"` for a 2-group collapse).
#' @param verbose Logical. If `TRUE` (default), prints a summary of how
#'   levels were recoded.
#'
#' @return The input `data` with a new factor column `new_col` appended,
#'   with levels in the order defined by `groups`.
#'
#' @details
#' Quantile level indices follow the ordering produced by [base::cut()]:
#' index `1` is always the lowest expression group and index `n` is the
#' highest. So `list(low = c(1, 2), high = c(3, 4))` maps the two lowest
#' quantiles to `"low"` and the two highest to `"high"`.
#'
#' @examples
#' set.seed(42)
#' my_data<-data.frame(
#'   id   = paste0("S",1:20),
#'   EGFR = rnorm(20,mean=2)
#' )
#' my_data <- get_quantile_groups(
#'   data        = my_data,
#'   features    = "EGFR",
#'   n_quantiles = 5,
#'   mode        = "quantile"
#' )
#' collapse_quantile_groups(
#'   data    = my_data,
#'   feature = "EGFR_gp5",
#'   groups  = list(low=c(1,2,3),high=c(4,5)),
#'   new_col = "EGFR_high_vs_low"
#' )
#'
#' @seealso [get_quantile_groups()], [determine_optimal_quantile()]
#'
#' @export
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

  #Recode Variables
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

#' Determine the optimal quantile split for survival analysis
#'
#' Systematically tests all possible quantile splits of a single feature
#' across a range of quantile counts, evaluating each split using a
#' log-rank test, with p-values are adjusted for multiple comparisons.
#' Returns the data with the best split applied and a full results table
#' ordered by q-value, then p-value, then n_quantiles.
#'
#' @param data A data frame or tibble containing the feature, time, and
#'   event variables.
#' @param feature A string. Name of the numeric feature column to split.
#' @param time_var A string. Name of the survival time column.
#' @param event_var A string. Name of the event indicator column
#'   (1 = event occurred, 0 = censored).
#' @param q_range An integer vector of quantile counts to test. Default
#'   is `2:6`, meaning splits from 2 to 6 quantiles are evaluated. A maximum of
#'   6 quantiles is used as default reduce the probability of
#'   generating very small comparison groups (e.g. < 15%).
#' @param verbose Logical. If `TRUE` (default), prints the best quantile split.
#' @param padj_method A string. Name of p-value correction method. Default method
#'   is `BH`. Possible methods are: `holm`,`hochberg`,`hommel`,`bonferroni`,
#'   `BH`,`BY`,`fdr`,`none`.
#'
#' @return The input `data` with two new columns appended:
#'     the intermediate quantile grouping (`feature_gpN`) and the
#'     optimal binary split (`feature_group`) with values `"Low"`
#'     and `"High"`, and an attribute named `optimal_split_results`,
#'     which consists in a data frame of all tested splits ordered by
#'     `q_value`, with columns: `n_quantiles`, `split_point`,`low_group`,
#'     `high_group`, `n_low_group`, `n_high_group`, `p_value`,`q_value`.
#'
#' @details
#' For each quantile count `q` in `q_range`, the feature is split into
#' `q` equal-frequency groups. All `q - 1` possible binary collapses are
#' then tested (e.g. for `q = 4`: 1vs3, 2vs2, 3vs1) using a log-rank test
#' via [survival::survdiff()]. The split with the lowest `q_value` is
#' selected and applied to the data. Ties are resolved by lower `p-value` and
#' then lower `n_quantile`. The `q_value` column represents the adjusted
#' `p-value` calculated using [stats::p.adjust()] using the desired method.
#'
#' Sometimes tied values or limited feature variability causes quantile
#' groups to be empty. Therefore, quantile splits that fail to produce
#' `"Low"` and `"High"` groups with at least one case are skipped.
#'
#' The output column `feature_group` always has two values: `"Low"` and
#' `"High"`, where `"Low"` contains the bottom quantiles and `"High"`
#' the top quantiles as determined by the best split.
#'
#' **Note:** This function performs multiple testing across all splits and
#' quantile counts. The resulting q-values and p-values should be interpreted
#' with caution and ideally validated in an independent cohort.
#'
#' @examples
#' \dontrun{
#' #Requires survival data with time and event columns
#' result<-determine_optimal_quantile(
#'  data      = my_data,
#'  feature   = "EGFR",
#'  time_var  = "surv_time",
#'  event_var = "event"
#' )
#' attr(result,"optimal_split_results")  #full table of tested splits
#' head(result)                          #data with optimal split applied
#' }
#'
#' @seealso [get_quantile_groups()], [collapse_quantile_groups()]
#'
#' @export
determine_optimal_quantile<-function(data,feature,time_var,event_var,
                                     q_range=2:6,padj_method="BH",verbose=TRUE) {

  #Safety checks
  padjust_methods<-c("holm","hochberg","hommel","bonferroni","BH","BY","fdr","none")
  if (!padj_method %in% padjust_methods) {
    stop("Invalid p-value adjustment method. Choose one of: ",
         paste(padjust_methods,collapse=", "))
  }
  if (!feature %in% names(data))
    stop("Feature '",feature,"' not found in data")

  results<-list()
  for ( q in q_range ) {

    #Create quantiles
    data_qt<-tryCatch(
      get_quantile_groups(
      data        = data,
      features    = feature,
      n_quantiles = q ,
      mode        = "quantile"
    ),error = function(e) NULL)
    if (is.null(data_qt)) next

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

      #Define test_id to align column names
      test_id<-paste( q , split , sep="_")

      if (length(grp_counts) < 2 || any(grp_counts==0)) {
        results[[test_id]]<-data.frame(
          n_quantiles   = q,
          split_point   = split,
          low_group     = paste0(split, "/", q),
          high_group    = paste0(q - split, "/", q),
          n_low_group   = NA,
          n_high_group  = NA,
          p_value       = 1.0,
          note          = "invalid split"
        )
        next
      }

      fit<-survival::survdiff(survival::Surv(
        data_qt_binary[[time_var]],data_qt_binary[[event_var]])~group,
        data=data_qt_binary)

      p_value<-stats::pchisq(fit$chisq, df = 1, lower.tail = FALSE)
      #df=1 means degrees of freadom of 1 (comparison between 2 groups: high vs low)

      results[[paste( q , split , sep="_")]]<-data.frame(
        n_quantiles   =  q ,
        split_point   = split,
        low_group     = paste0(split, "/", q ),
        high_group    = paste0( q - split, "/", q ),
        n_low_group   = sum(data_qt_binary$group=="Low" , na.rm = TRUE),
        n_high_group  = sum(data_qt_binary$group=="High", na.rm = TRUE),
        p_value       = p_value,
        note          = "valid split"
      )
    }
  }
  results_data<-dplyr::bind_rows(results) |>
    dplyr::filter(.data$note == "valid split") |>
    dplyr::mutate(q_value = stats::p.adjust(.data$p_value,method = padj_method)) |>
    dplyr::arrange(.data$q_value,.data$p_value,.data$n_quantiles)

  #Safety check
  if (nrow(results_data) == 0) {
    stop("No valid quantile split could be generated.",
         "Feature may contain too few unique values.")
  }

  #Apply best split to the actual data
  best<-results_data[1,]
  if(verbose) cat("\nBest split: ",best$low_group,
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


