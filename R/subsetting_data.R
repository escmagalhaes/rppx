### Subsetting data according to continuous and categorical variables ###

#' Generate logical group masks from a categorical variable
#'
#' Creates a named list of logical vectors (masks) that identify which rows
#' belong to each group defined by a categorical variable. Works on any column
#' whose values represent category labels, including character columns, factor
#' columns, and discretized numeric columns produced by [get_quantile_groups()].
#' Raw numeric columns are not supported. Masks can be used directly to subset
#' data frames via [get_subsets()].
#'
#' @param data A data frame or tibble containing the variable to group by.
#' @param var A string. Name of the categorical column to group by.Can be a
#' character, factor, or quantile factor column from [get_quantile_groups()].
#' @param var_labs A string. Display label for the variable, used as prefix
#'   in group names (e.g. `"MUT_status"` to `"MUT_status_pos"`). If `NULL`
#'   (default), `var` is used as the prefix.
#' @param ref_patterns A character vector of patterns matched against values
#'   in `var` to identify the reference/negative group. Matching is
#'   case-insensitive and whole-word. If `NULL` (default), no reference
#'   group is identified and only `all`, `missing`, and individual level
#'   masks are returned. Use `NULL` for variables with no natural reference
#'   (e.g. therapy arms, quantile groups).
#' @param ref_suffix A string. Suffix for the reference group name.
#'   Only used when `ref_patterns` is not `NULL`. Default is `NULL`.
#' @param alt_suffix A string. Suffix for the combined alternative group
#'   name. Only used when `ref_patterns` is not `NULL`. Default is `NULL`.
#' @param missing_suffix A string. Suffix for the missing values group
#'   name. Default is `"missing"`.
#' @param verbose Logical. If `TRUE` (default), prints a message when no
#'   alternative levels are detected.
#'
#' @return A named list of logical vectors, each of length `nrow(data)`:
#' \describe{
#'   \item{`all`}{All rows (`TRUE` for every row).}
#'   \item{`prefix_missing_suffix`}{Rows with `NA` in `var`.}
#'   \item{`prefix_ref_suffix`}{Reference group. Only present when`ref_patterns`
#'    is not `NULL`.}
#'   \item{`prefix_alt_suffix`}{All non-reference, non-missing rows combined.
#'    Only present when `ref_patterns` is not `NULL`.}
#'   \item{`prefix_value`}{One mask entry per unique non-reference value,
#'     with spaces replaced by `_` and special characters removed.}
#' }
#'
#' @details
#' When `ref_patterns = NULL`, the function skips reference group detection
#' entirely and returns `all`, `missing`, and one mask per unique value.
#' This is the correct behavior for variables with no natural reference
#' group, such as therapy arms, quantile interval factors, or any
#' nominal variable where all levels are equally valid comparators.
#'
#' When `ref_patterns` is provided, whole-word regex matching is used to detect
#' the values of `ref_patterns`. Values not matching any reference pattern and
#' not missing are treated as alternative levels. Each unique alternative value
#' gets its own mask in addition to the combined `alt_suffix` mask.
#'
#' Group names use `var_labs` as prefix, making them traceable back to
#' the source variable. This naming convention is relied upon by
#' [get_subsets()] and the wrapper [create_all_subsets()].
#'
#' @examples
#' my_data<-data.frame(
#'   mut_flt3=c("WT","ITD","wt","TKD",NA,"WT"),
#'   tx  = c("AraC","HMA","venetoclax","AraC","HMA",NA),
#'   stringsAsFactors=FALSE
#' )
#'
#' #With reference group (Binary variable e.g. mutation status)
#' get_var_groups(
#'   data         = my_data,
#'   var          = "mut_flt3",
#'   var_labs     = "FLT3",
#'   ref_patterns = c("WT","wt"),
#'   ref_suffix   = "WT",
#'   alt_suffix   = "mutated"
#' )
#'
#' #Without reference group (Categorical variable e.g. therapy)
#' get_var_groups(
#'   data     = my_data,
#'   var      = "tx",
#'   var_labs = "Therapy"
#' )
#'
#' @seealso [get_quantile_groups()], [get_subsets()], [create_all_subsets()]
#'
#' @export
get_var_groups<-function(data,var,var_labs=NULL
                         ,ref_patterns=NULL,ref_suffix=NULL,alt_suffix=NULL
                         ,missing_suffix="missing",verbose=TRUE) {

  #Safety check
  if (is.null(ref_patterns) && (!is.null(ref_suffix) || !is.null(alt_suffix)))
    warning("ref_suffix and alt_suffix are not used when ref_patterns is NULL")

  #Use var name as label if none provided
  if(is.null(var_labs)) var_labs<-var

  is_missing <-is.na(data[[var]])
  status_vec <-trimws(as.character(data[[var]]))
  prefix     <-paste0(var_labs,"_")
  var_groups <- list()

  #All samples
  var_groups$all<-rep(TRUE,length(status_vec))

  #Missing values
  var_groups[[paste0(prefix,missing_suffix)]]<-is_missing

  #Reference group detection
  if (!is.null(ref_patterns)) {
    pattern <-paste0("\\b(",paste(ref_patterns,collapse="|"),")\\b")
    is_neg  <-grepl(pattern,status_vec,ignore.case=TRUE)
    var_groups[[paste0(prefix,ref_suffix)]] <-is_neg
    idx     <-(!is_neg) & !is_missing
  } else {
    is_neg <-rep(FALSE,length(status_vec))
    idx    <-!is_missing
  }

  #Individual levels
  var_values<-unique(status_vec[idx])
  var_values<-var_values[!is.na(var_values)]

  if(length(var_values)==0){
    if(verbose) message("No alternative levels detected for variable: ",var)
    return(var_groups)
  }

  #Combined alternative group (only when ref_patterns provided)
  if (!is.null(ref_patterns)) {
    var_groups[[paste0(prefix,alt_suffix)]]<-idx
  }

  #Individual levels masks
  for(val in var_values){
    safe_val<-gsub("\\s+","_",val)
    safe_val<-gsub("[^A-Za-z0-9_]","",safe_val)
    var_groups[[paste0(prefix,safe_val)]]<-(status_vec==val) & !is_missing
  }
  return(var_groups)
}

#' Create subsets from two group logical vector lists (masks)
#'
#' Takes exactly two named lists of logical vector lists (masks) and produces
#' all non-empty pairwise combinations as subsetted data frames. Designed
#' to work with the output of [get_var_groups()].
#'
#' @param data A data frame or tibble to subset.
#' @param group_lists A named list of exactly 2 named lists of logical
#'   vectors. Each logical vector must have length `nrow(data)`.
#'   Typically produced by [get_var_groups()].
#' @param print_summary Logical. If `TRUE`, prints a summary table of
#' subset names and sizes ordered by number of rows.
#' @param exclude_NA Logical. If `TRUE` (default), excludes group
#'   combinations whose name contains `NA` or the `missing_suffix`.
#' @param missing_suffix String. Suffix for identifying missing values,
#'   used when `exclude_NA = TRUE`. Should match the `missing_suffix` used in
#'   [get_var_groups()]. Default is `"missing"`.
#'
#' @return A named list with two elements:
#' \describe{
#'   \item{`subsets`}{A named list of data frames, one per non-empty
#'     pairwise combination. Names follow the pattern
#'     `"group1name_group2name"`.}
#'   \item{`summary`}{A data frame with columns `subset_name` and
#'     `n_rows`, documenting the size of each subset.}
#' }
#'
#' @details
#' All pairwise combinations of masks from the two group lists are
#' evaluated. Empty combinations (0 rows) are silently dropped.
#'
#' The combination of `all` from both group lists is renamed to `"all"`
#' rather than `"all_all"`, representing the full dataset as the baseline
#' comparison group.
#'
#' When using [get_subsets()] with two group lists derived from the same
#' feature (e.g. expression level and mutation status of the same protein),
#' use short distinct `var_labs` for each call to avoid redundant subset
#' names (see @examples ).
#'
#' `exclude_NA = TRUE` excludes combinations whose name contains
#' `"missing"` or `"NA"`,  matching the naming convention of
#' [get_var_groups()]. Actual `NA` values in the data are not tracked.
#'
#' Currently limited to exactly 2 group lists. For single-dimension
#' subsetting or more than 2 dimensions, use [create_all_subsets()].
#'
#' @examples
#' set.seed(123)
#' my_data<-data.frame(
#'  FLT3=rnorm(20,mean=0.1),
#'  mut_flt3=rep(c("WT","ITD","wt","TKD",NA),4)
#'  )
#'
#' normal_data<-data.frame(FLT3=rnorm(20,mean=0))
#'
#' ptn_sets<-get_quantile_groups(
#'  data          = my_data,
#'  features      = "FLT3",
#'  data_normal   = normal_data,
#'  mode          = "normal_based",
#'  normal_labels = paste(c("below","lower","upper","above"),"normal",sep="_")
#'  )
#'
#' ptn_var_groups<-get_var_groups(
#'   data         = ptn_sets,
#'   var          = "FLT3_gp_n",
#'   var_labs     = "expr"
#' )
#'
#' mut_var_groups<-get_var_groups(
#'   data         = my_data,
#'   var          = "mut_flt3",
#'   var_labs     = "FLT3",
#'   ref_patterns = c("WT","wt"),
#'   ref_suffix   = "WT",
#'   alt_suffix   = "mutated"
#' )
#'
#' result <- get_subsets(
#'   data          = my_data,
#'   group_lists   = list(ptn = ptn_var_groups, mut = mut_var_groups),
#'   print_summary = FALSE
#' )
#' result$summary
#' result$subsets[["expr_lower_normal_FLT3_mutated"]]
#'
#' @seealso [get_quantile_groups()], [get_var_groups()], [create_all_subsets()]
#'
#' @export
get_subsets<-function(data,group_lists,print_summary=FALSE,
                      exclude_NA=TRUE,missing_suffix="missing") {

  if(length(group_lists) !=2)
    stop("group_lists must contain exactly 2 named lists")

  #Check all vectors in both group lists match nrow(data)
  n<-nrow(data)
  for( g in group_lists) {
    for(nm in names( g )) {
      if(length( g[[ nm ]]) != n )
        stop("Vector '",nm,"' has length ",length( g[[ nm ]]),
             " but data has ", n ," rows")
    }
  }

  #Helper to define names from possible NAs
  is_na_group_name<-function(x,missing_suffix) {
    grepl(paste0("(^|_)(NA|",missing_suffix,")(_|$)"), x)
  }

  g1 <-group_lists[[1]]
  g2 <-group_lists[[2]]
  subset_list<-list()
  summary_table<-data.frame(subset_name=character(),n_rows=integer()
                            ,stringsAsFactors=FALSE)
  for( n1 in names( g1 )) {
    if(exclude_NA && is_na_group_name(n1,missing_suffix)) next
    v1 <- g1[[ n1 ]]
    for( n2 in names( g2 )) {
      if(exclude_NA && is_na_group_name(n2,missing_suffix)) next
      v2 <- g2[[ n2 ]]
      subset_name<-paste0( n1 ,"_", n2 )

      #Rename all_all subset to just "all" for clarity
      if (n1=="all" && n2=="all") subset_name<-"all"

      mask<- v1 & v2
      mask[is.na(mask)]<-FALSE
      n_rows<-sum(mask)
      if(n_rows==0) next
      subset_list[[subset_name]]<-data[mask,,drop = FALSE]
      summary_table<-rbind(summary_table,data.frame(
        subset_name=subset_name,n_rows=n_rows))
    }
  }
  if(print_summary && nrow(summary_table) > 0) {
    cat("Subset sizes (total subsets =",nrow(summary_table),"):\n")
    summary_sorted<-summary_table[order(summary_table$n_rows),]
    for(i in seq_len(nrow(summary_sorted))) {
      cat(sprintf("  %-30s : %d rows\n",
                  summary_sorted$subset_name[ i ],
                  summary_sorted$n_rows[ i ]))
    }
  }
  return(list(subsets = subset_list,
              summary = summary_table))
}

#' Generate all subsets across two groups of variables
#'
#' Wrapper function that loops over all combinations of two variable groups,
#' calling [get_quantile_groups()], [get_var_groups()], [get_subsets()] and
#' [collapse_quantile_groups()] (if called) internally. Returns a nested list
#' of all non-empty subsets organized by variable combination.
#' Both variable groups can contain categorical or continuous variables.
#' Continuous variables are discretized first via [get_quantile_groups()].
#'
#' @param data A data frame or tibble to subset.
#' @param var1 A character vector of variable names for the first
#'   subsetting dimension.
#' @param var2 A character vector of variable names for the second
#'   subsetting dimension. If `NULL`, single-dimension subsetting is
#'   applied using `var1` only.
#' @param continuous_vars A character vector of variable names (from
#'   `var1` and/or `var2`) that require discretization via
#'   [get_quantile_groups()] before masking. If `NULL`, all variables
#'   are treated as already categorical.
#' @param data_normal A data frame of normal samples used to define
#'   reference-based interval boundaries. Required when `continuous_vars`
#'   is provided and `quantile_mode = "normal_based"`.
#' @param quantile_mode Character. One of `"normal_based"` (default) or
#'   `"quantile"`. Passed to [get_quantile_groups()]. Only relevant when
#'   `continuous_vars` is provided.
#' @param n_quantiles A single integer. Passed to [get_quantile_groups()].
#'   Required when `continuous_vars` is given and `quantile_mode = "quantile"`.
#' @param var1_labs A named list mapping names in `var1` to display
#'   labels used as prefix in mask names. If `NULL`, column names are
#'   used directly.
#' @param var2_labs A named list mapping names in `var2` to display
#'   labels. If `NULL`, column names are used directly.
#' @param ref_patterns A named list mapping variable names to reference
#'   patterns. Each element is a character vector passed to
#'   [get_var_groups()]. Variables absent from the list default to
#'   `NULL` (no reference group). Continuous variables are always
#'   masked with `ref_patterns = NULL`.
#' @param ref_suffixes A named list mapping variable names to reference
#'   group suffix strings. Defaults to `"neg"` for variables not in
#'   the list.
#' @param alt_suffixes A named list mapping variable names to alternative
#'   group suffix strings. Defaults to `"pos"` for variables not in
#'   the list.
#' @param missing_suffix A string. Suffix for missing value masks.
#'   Applied to all variables. Default is `"missing"`.
#' @param add_binary_sets Logical. If `TRUE` (default), adds combined
#'   `<label>_low` and `<label>_high` masks to continuous variable group
#'   lists by splitting quantile levels at `binary_split`. These are skipped
#'   when there are less than 3 groups.
#' @param binary_split A single integer. Level index at which split into `low`
#'   and `high` when `add_binary_sets = TRUE`. Default is `2`.
#' @param exclude_NA Logical. Passed to [get_subsets()]. Default `TRUE`.
#' @param print_summary Logical. Passed to [get_subsets()]. Default `FALSE`.
#'
#' @return A named list where each element corresponds to a variable
#'   combination (named `"var1_var2"`) and contains:
#' \describe{
#'   \item{`subsets`}{Named list of subsetted data frames.}
#'   \item{`summary`}{Data frame with columns `subset_name` and `n_rows`.}
#' }
#'
#' @details
#' `var1` and `var2` define the two subsetting dimensions. All pairwise
#' combinations across `var1 × var2` are evaluated. When `var2 = NULL`,
#' single-dimension subsetting is applied to `var1` only.
#'
#' Variables listed in `continuous_vars` are discretized via
#' [get_quantile_groups()] before being passed to [get_var_groups()].
#' If `add_binary_sets = TRUE` combined `low` and `high` masks are created
#' with [collapse_quantile_groups()], according to the `binary_split` value.
#' All other variables go directly to [get_var_groups()].
#'
#' `ref_patterns`, `ref_suffixes`, and `alt_suffixes` are named lists.
#' Each variable can have its own settings. Variables absent from these
#' lists use defaults (`NULL`, `"neg"`, `"pos"` respectively). Continuous
#' variables always use `ref_patterns = NULL` regardless of what is
#' specified.
#'
#' Note that continuous variables must be listed in `var1`/`var2` , which
#' determine the pairs iterated over, and also in `continuous_vars`, which
#' controls discretization. See @examples .
#'
#' @examples
#' \dontrun{
#' # Continuous x categorical
#' all_subsets<-create_all_subsets(
#'   data            = disease_data,
#'   data_normal     = normal_data,
#'   var1            = c("EGFR", "AKT1"),
#'   var2            = c("mut_flt3", "mut_npm1"),
#'   continuous_vars = c("EGFR", "AKT1"),
#'   var1_labs       = list(EGFR = "EGFR_expr", AKT1 = "AKT1_expr"),
#'   var2_labs       = list(mut_flt3 = "FLT3", mut_npm1 = "NPM1"),
#'   ref_patterns    = list(mut_flt3 = c("WT","wt"),
#'                          mut_npm1 = c("WT","wt")),
#'   ref_suffixes    = list(mut_flt3 = "WT", mut_npm1 = "WT"),
#'   alt_suffixes    = list(mut_flt3 = "mutated", mut_npm1 = "mutated")
#' )
#' all_subsets[["EGFR_mut_flt3"]]$summary
#'
#' # Two categoricals
#' cat_subsets <- create_all_subsets(
#'   data         = disease_data,
#'   var1         = c("mut_npm1", "mut_flt3"),
#'   var2         = c("tx", "cyto_risk"),
#'   var1_labs    = list(mut_npm1 = "NPM1", mut_flt3 = "FLT3"),
#'   var2_labs    = list(tx = "Therapy", cyto_risk = "CytoRisk"),
#'   ref_patterns = list(mut_npm1 = c("WT","wt"),
#'                       mut_flt3 = c("WT","wt"))
#' )
#' }
#'
#' @seealso [get_var_groups()], [get_quantile_groups()] and [get_subsets()]
#'
#' @export
create_all_subsets<-function(data,
                             var1,
                             var2             = NULL,
                             continuous_vars  = NULL,
                             data_normal      = NULL,
                             quantile_mode    = c("normal_based", "quantile"),
                             n_quantiles      = NULL,
                             var1_labs        = NULL,
                             var2_labs        = NULL,
                             ref_patterns     = NULL,
                             ref_suffixes     = NULL,
                             alt_suffixes     = NULL,
                             missing_suffix   = "missing",
                             add_binary_sets  = TRUE,
                             binary_split     = 2,
                             exclude_NA       = TRUE,
                             print_summary    = FALSE) {

  #Safety checks
  if (missing(var1) || is.null(var1))
    stop("'var1' must be provided")

  #Validate quantile_mode arguments only when continuous_vars provided
  if (!is.null(continuous_vars)) {

    quantile_mode<-match.arg(quantile_mode)
    if (quantile_mode=="quantile" && is.null(n_quantiles))
      stop("n_quantiles is required when quantile_mode='quantile'")
    if (quantile_mode=="normal_based" && is.null(data_normal))
      stop("data_normal is required when quantile_mode='normal_based'")

    if (add_binary_sets) {
      n_lev<-if (quantile_mode=="normal_based") 4L else n_quantiles
      if (length(n_lev) != 1)
        stop("'n_quantiles' must be a single integer")

      if (n_lev < 3) {
        message("Only ",n_lev," groups per continuous variable: ",
                "'low'/'high' masks would be redundant and were skipped")
        add_binary_sets<-FALSE
      } else if (binary_split < 1 || binary_split >= n_lev) {
        stop("'binary_split' must be between 1 and ",n_lev - 1,
             " for ",n_lev," groups")
      }
    }

    unknown_cont_var<-setdiff(continuous_vars,c(var1,var2))
    if (length(unknown_cont_var) > 0)
      stop("These continuous_vars are not var1 or var2: ",
           paste(unknown_cont_var,collapse=", "))
  }

  #Helper to add high/low binary masks
  add_binary_masks<-function(data,quantile_col,binary_split) {

    #Safety checks
    if (!quantile_col %in% names(data))
      stop("Column '",quantile_col,"' not found in data")
    if (!is.factor(data[[quantile_col]]))
      stop("Column '",quantile_col,"' must be a factor. ",
           "Use get_quantile_groups() first")

    n_levels<-nlevels(data[[quantile_col]])
    if (binary_split < 1 || binary_split >= n_levels)
      stop("binary_split must be between 1 and ",n_levels - 1,
           " for a ",n_levels,"-level factor")

    binary_col<-paste0(quantile_col,"_binary")

    collapsed<-collapse_quantile_groups(
      data    = data,
      feature = quantile_col,
      groups  = list(low  = 1:binary_split,
                     high = (binary_split + 1):n_levels),
      new_col = binary_col,
      verbose = FALSE
    )

    list(
      low  = !is.na(collapsed[[binary_col]]) & collapsed[[binary_col]]=="low",
      high = !is.na(collapsed[[binary_col]]) & collapsed[[binary_col]]=="high"
    )
  }

  #Helper function to build masks for a single variable
  build_masks<-function(var,labs) {
    is_continuous<-!is.null(continuous_vars) && var %in% continuous_vars

    if (is_continuous) {
      #Discretize first
      disc_data<-get_quantile_groups(
        data        = data,
        features    = var,
        data_normal = data_normal,
        mode        = quantile_mode,
        n_quantiles = n_quantiles
      )
      quantile_col<- if (quantile_mode=="normal_based") {
        paste0(var,"_gp_n")
      } else {
        paste0(var,"_gp",n_quantiles)
      }
      var_lab<-if (!is.null(labs) && !is.null(labs[[var]])) labs[[var]] else var

      masks<-get_var_groups(
        data           = disc_data,
        var            = quantile_col,
        var_labs       = var_lab,
        ref_patterns   = NULL,
        ref_suffix     = NULL,
        alt_suffix     = NULL,
        missing_suffix = missing_suffix
      )

      if (add_binary_sets) {
        binary_masks <-add_binary_masks(disc_data,quantile_col,binary_split)
        masks[[paste0(var_lab,"_low")]]  <-binary_masks$low
        masks[[paste0(var_lab,"_high")]] <-binary_masks$high
      }

    } else {

      #Categorical: direct masking
      var_lab     <-if (!is.null(labs)) labs[[var]] else var

      var_ref     <-if (!is.null(ref_patterns)) ref_patterns[[var]] else NULL

      var_ref_suf <-if (!is.null(var_ref)) {
        if (!is.null(ref_suffixes)) ref_suffixes[[var]] else "neg"
      } else NULL

      var_alt_suf <- if (!is.null(var_ref)) {
        if (!is.null(alt_suffixes)) alt_suffixes[[var]] else "pos"
      } else NULL

      masks<-get_var_groups(
        data           = data,
        var            = var,
        var_labs       = var_lab,
        ref_patterns   = var_ref,
        ref_suffix     = var_ref_suf,
        alt_suffix     = var_alt_suf,
        missing_suffix = missing_suffix
      )
    }
    masks
  }

  #Helper to define names from possible NAs
  is_na_group_name<-function(x,missing_suffix) {
    grepl(paste0("(^|_)(NA|",missing_suffix,")(_|$)"), x)
  }

  #Build iteration lists
  var2_use    <-if (is.null(var2)) "" else var2
  all_results <-list()

  for ( v1 in var1) {
    for ( v2 in var2_use) {

      group_lists<-list()

      #Build var1 masks
      group_lists$var1<-build_masks(v1,var1_labs)

      #Build var2 masks
      if (v2 != "") {
        group_lists$var2<-build_masks(v2,var2_labs)
      }

      #Create subsets
      if (length(group_lists)==1) {
        g1<-group_lists[[1]]
        if (exclude_NA) g1<-g1[!is_na_group_name(names(g1),missing_suffix)]
        subsets       <-lapply(g1, function(mask) data[mask,,drop=FALSE])
        subsets       <-Filter(function(s) nrow(s) > 0,subsets)
        summary_table <-data.frame(
          subset_name      = names(subsets),
          n_rows           = vapply(subsets,nrow,integer(1)),
          stringsAsFactors = FALSE
        )
      } else {
        result<-get_subsets(
          data           = data,
          group_lists    = group_lists,
          exclude_NA     = exclude_NA,
          missing_suffix = missing_suffix,
          print_summary  = print_summary
        )
        subsets       <-result$subsets
        summary_table <-result$summary
      }

      #Store results
      key<-paste(c(v1,v2)[c(v1,v2) !=""],collapse="_")
      all_results[[key]]<-list(subsets=subsets,summary=summary_table)
    }
  }
  return(all_results)
}


