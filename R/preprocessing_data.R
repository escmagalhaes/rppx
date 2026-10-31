### Preprocessing raw data ###

#' Get unique values from variables in a dataset
#'
#' Extracts unique values from one or more columns in a data frame,
#' returning a structured table for manual inspection and cleaning.
#' Useful for auditing raw categorical variables before recoding.
#'
#' @param data A data frame or tibble containing the variables to inspect.
#' @param features A character vector of column names to extract values from.
#'
#' @return A data frame with three columns:
#' \describe{
#'   \item{variable}{Column name from the input data.}
#'   \item{raw_value}{Unique values found in that column, lowercased,
#'     trimmed of whitespace, sorted, and with `NA` removed.}
#'   \item{clean_value}{Empty column (`NA`) for manual completion:
#'     fill this in to map raw values to standardized labels.}
#' }
#'
#' @details
#' The `clean_value` column is intentionally empty. The intended workflow
#' is to export the result, fill in `clean_value` manually, and use the
#' completed table as a recoding reference for downstream cleaning.
#'
#' @examples
#' dataframe<-data.frame(
#'   status = c("Alive","alive","Dead ","dead",NA),
#'   grade  = c("low","Low","HIGH","high","medium")
#' )
#' get_rawtable_values(data=dataframe,features=c("status","grade"))
#'
#' @seealso [safe_clean_raw_data()]
#'
#' @export
get_rawtable_values<-function(data,features) {

  #Safety checks
  if (!is.data.frame(data))
    stop("data must be a dataframe or tibble.")

  if (missing(features) || !is.character(features))
    stop("features must be a character vector of column names.")

  unknown_features<-setdiff(features,names(data))
  if (length(unknown_features) > 0)
    stop("Some variable(s) were not found in data: ",paste(unknown_features,collapse=", "))

  #Extract unique values
  data_rows<-lapply(features,function(var) {
    vals<-unique(tolower(trimws(as.character(data[[ var ]]))))
    vals<-sort(vals[!is.na(vals) & vals !="na"])
    data.frame(
      variable    = var,
      raw_value   = vals,
      clean_value = NA_character_,
      stringsAsFactors = FALSE )
  })
  do.call(rbind,data_rows)
}

#' Clean and recode raw variables using a corrections table
#'
#' Applies standardized value recoding to categorical variables using a
#' corrections table, and coerces common NA strings to `NA` in numeric
#' variables. Returns the cleaned data frame with a full change report
#' attached as an attribute.
#'
#' @param data A data frame or tibble containing the variables to clean.
#' @param corrections_table A data frame with three columns: `variable`
#'   (column name in `data`), `raw_value` (original value to replace),
#'   and `clean_value` (replacement value). Typically the
#'   output of [get_rawtable_values()] after manual completion.
#' @param numeric_vars A character vector of column names to treat as
#'   numeric. These are not recoded but have common NA strings coerced
#'   to `NA`. Default is `NULL`.
#' @param na_strings A character vector of strings to treat as `NA`.
#'   Matching is case-insensitive after whitespace trimming. Defaults
#'   to a set of common NA representations.
#'
#' @return The input `data` with cleaned variables, plus an attribute
#'   `safe_clean_report` which is a data frame that documents all changes:
#' \describe{
#'   \item{variable}{Column name that was processed.}
#'   \item{row}{Row index of the changed value.}
#'   \item{raw_value}{Original value before cleaning.}
#'   \item{cleaned_value}{Value after cleaning, or `NA` if coerced.}
#'   \item{status}{One of `"cleaned"` (recoded), `"unmatched"` (value
#'     not in corrections table), or `"coerced_na"` (NA string in a
#'     numeric variable).}
#' }
#'
#' @details
#' This function is designed to work in a two-step workflow with
#' [get_rawtable_values()]:
#' 1. Run [get_rawtable_values()] to extract unique raw values
#' 2. Fill in `clean_value` manually in the resulting table
#' 3. Pass the completed table here as `corrections_table`
#'
#' Unmatched values (present in `data` but not in `corrections_table`)
#' are left in lower case and with trimmed whitespace, and flagged in
#' the report with status `"unmatched"`. Values matching `na_strings` in
#' numeric variables are coerced to `NA` and flagged as `"coerced_na"`.
#' Retrieve the full report with `attr(result,"safe_clean_report")`.
#'
#' Variables appearing in both `corrections_table` and `numeric_vars`
#' trigger a warning stating that they will be processed as categorical
#' variables.
#'
#' @examples
#' my_data<-data.frame(
#'   status = c("Alive", "alive", "Dead ", "dead", NA),
#'   age    = c("45", "n/a", "67", "unknown", "52")
#' )
#' corrections<-data.frame(
#'   variable    = c("status", "status"),
#'   raw_value   = c("alive", "dead"),
#'   clean_value = c("Alive", "Dead")
#' )
#' result<-safe_clean_raw_data(
#'   data              = my_data,
#'   corrections_table = corrections,
#'   numeric_vars      = "age"
#' )
#' attr(result,"safe_clean_report")
#'
#' @seealso [get_rawtable_values()]
#'
#' @export
safe_clean_raw_data<-function(
    data,corrections_table,numeric_vars=NULL,
    na_strings=c("na","n/a","n.a","nan","null","none","missing","not available",
                 "not applicable","unknown","","--","---")) {

  #Safety checks and warnings
  if (!is.data.frame(data))
    stop("data must be a dataframe or tibble.")

  if (!is.data.frame(corrections_table))
    stop("corrections_table must be a dataframe.")

  na_clean<-is.na(corrections_table$clean_value)
  if (any(na_clean))
    warning("\nThese character(s) were not assigned a clean_value will remain unchanged:\n",
            paste(corrections_table$variable[na_clean],'==',
                  corrections_table$raw_value[na_clean],'\n'),call.=FALSE)

  unknown_vars<-setdiff(unique(corrections_table$variable),names(data))
  if (length(unknown_vars) > 0)
    warning("These variable(s) in corrections_table were not found in data and were skipped: ",
            paste(unknown_vars,collapse=", "),call.=FALSE)

  if (!is.null(numeric_vars)) {
    unknown_numeric<-setdiff(numeric_vars,names(data))
    if (length(unknown_numeric) > 0)
      warning("These variable(s) in numeric_vars were not found in data and were skipped: ",
              paste(unknown_numeric,collapse=", "), call.=FALSE)
  }

  overlap_vars<-intersect(unique(corrections_table$variable),numeric_vars)
  if (length(overlap_vars) > 0)
    warning("These variable(s) were found in both corrections_table and numeric_vars: ",
            paste(overlap_vars,collapse=", "),call.=FALSE)

  #Processing categorical/binary variables
  report_rows   <-list()
  result        <-data
  vars_to_clean <-intersect(unique(corrections_table$variable),names(data))

  for (var in vars_to_clean) {

    #Get var values
    raw_var_values<-as.character(data[[ var ]])

    #Lower case and NA coercion with white space trimming
    variable_values<-ifelse(is.na(raw_var_values),NA_character_,
                            tolower(trimws(raw_var_values)))
    variable_values[variable_values %in% na_strings]<-NA_character_

    #Corrections table lookup
    var_correction<-corrections_table[corrections_table$variable == var ,]
    for (i in seq_len(nrow(var_correction))) {
      match_var_value<- !is.na(variable_values) & variable_values == var_correction$raw_value[ i ]
      variable_values[match_var_value]<-var_correction$clean_value[ i ]
    }

    #Report changes
    changed   <- !is.na(raw_var_values) & (is.na(variable_values) | tolower(trimws(raw_var_values)) != variable_values)
    unmatched <- !is.na(variable_values) & !variable_values %in% var_correction$clean_value

    idx<-which(changed | unmatched)
    for ( i in idx) {
      report_rows[[length(report_rows)+1]]<-data.frame(
        variable         = var,
        row              = i ,
        raw_value        = raw_var_values[ i ],
        cleaned_value    = variable_values[ i ],
        status           = ifelse(unmatched[ i ],"unmatched","cleaned"),
        stringsAsFactors = FALSE
      )
    }
    result[[ var ]]<-variable_values
  }

  #Processing numeric variables
  numeric_vars_clean <-intersect(numeric_vars,names(data))

  for (var in numeric_vars_clean) {
    raw_var_values<-as.character(data[[ var ]])
    variable_values<-ifelse(is.na(raw_var_values),NA_character_,
                            tolower(trimws(raw_var_values)))
    variable_values[variable_values %in% na_strings]<-NA_character_

    changed<- !is.na(raw_var_values) & is.na(variable_values)
    idx<-which(changed)
    for ( i in idx) {
      report_rows[[length(report_rows)+1]]<-data.frame(
        variable      = var ,
        row           = i ,
        raw_value     = raw_var_values[ i ],
        cleaned_value = NA_character_,
        status        = "coerced_na",
        stringsAsFactors = FALSE
      )
    }
    result[[ var ]]<-variable_values
  }

  #Compile Report
  report<- if (length(report_rows) > 0) {
    do.call(rbind,report_rows)
  } else {
    data.frame(variable=character(),row=integer(),raw_value=character(),
               cleaned_value=character(),status=character(),
               stringsAsFactors=FALSE)
  }
  if (any(report$status=="unmatched"))
    warning(sum(report$status=="unmatched")," unmatched value(s).\n
            Check attr(result,'safe_clean_report')", call. = FALSE)

  attr(result,"safe_clean_report")<-report
  return(result)
}

#' Extract a feature matrix from a dataset
#'
#' Gets an ID column and a contiguous block of feature columns from a data frame,
#' given that the expected number of features matches the input data.
#' Originally designed for protein data, but works in any dataset where
#' features occupy a contiguous right-hand block of columns.
#'
#' @param data A data frame or tibble with both sample ID and feature columns.
#' @param n_features A single integer. Expected number of feature columns.
#' @param id_column A string. Name of the column containing sample
#'   identifiers (e.g. patient IDs). Retained in the output.
#' @param start_column A string. Name of the first feature column.
#'   All columns from this position to the last column are extracted.
#'
#' @return A data frame containing the `id_column` and all feature
#'   columns from `start_column` to the last column.
#'
#' @details
#' The function assumes features occupy a contiguous block of columns
#' starting at `start_column` and ending at the last column of data.
#' Columns that identify samples must appear before `start_column`.
#' Useful to create a feature-only dataset to compare against reference
#' features and filter out features using [remove_extra_terms()]
#'
#' The `n_features` check is intentional: it acts as a safeguard
#' against silent dataset changes (e.g. columns added or removed
#' upstream). If the count does not match, the function stops with
#' an informative error.
#'
#' @examples
#' my_data <- data.frame(
#'   patient_id = c("P1","P2","P3"),
#'   age        = c(45,67,52),
#'   ABL1       = c(1.2,0.8,1.5),
#'   AKT1       = c(0.9,1.1,0.7)
#' )
#' get_feature_dataset(
#'   data         = my_data,
#'   n_features   = 2,
#'   id_column    = "patient_id",
#'   start_column = "ABL1"
#' )
#'
#' @seealso [remove_extra_terms()]
#'
#' @export
get_feature_dataset<-function(data,n_features,id_column,start_column) {

  #Safety checks
  if (!is.data.frame(data))
    stop("data must be a data.frame or tibble.")

  if (missing(n_features) || !is.numeric(n_features) || length(n_features) != 1)
    stop("n_features is required and must be a single number.")

  if (!id_column %in% names(data))
    stop("Patient ID column (",id_column, ") not found in dataframe")

  if (!start_column %in% names(data))
    stop("Feature start column(", start_column,") not found in dataframe",
         "Has the first feature column been renamed or moved?")

  #Detect feature columns in dataframe
  first_feature_position<-which(names(data)==start_column)
  n_feature_data<-ncol(data) - first_feature_position + 1

  if (n_feature_data != n_features)
    stop(sprintf("Feature count mismatch: expected %d, found %d (from '%s' to last column).\n",
                 n_features,n_feature_data,start_column))

  #Extract id and all feature columns
  feature_cols<-names(data)[first_feature_position:ncol(data)]
  data_feature<-data |> dplyr::select(dplyr::all_of(c(id_column,feature_cols)))

  return(data_feature)
}

#' Add factor levels to a metadata table
#'
#' Identify features that are classified as factors in a data frame and
#' appends their levels to the corresponding rows in a metadata table.
#' Designed to enrich a manually curated metadata table after data has
#' been loaded and factors have been determined.
#'
#' @param data A data frame or tibble containing the variables described
#'   in `meta_table`. Factor columns will be inspected for their levels.
#' @param meta_table A data frame with at least a `variable` column
#'   containing column names from `data`. Expected to also contain
#'   `label`, `type`, and `class` columns. Typically created manually
#'   to describe all variables in the dataset (see Details below).
#'
#' @return The input `meta_table` with a new `factor_levels` column
#'   added, containing the factor levels for each variable as a
#'   pipe-separated string (e.g. `"low|medium|high"`). Non-factor
#'   variables receive `NA`.
#'
#' @details
#' The metadata table (`meta_table`) is a manually curated data frame
#' that acts as the central registry of all variables in the dataset.
#' Beyond `variable` and `label`, it typically contains columns such as:
#' \itemize{
#'   \item `type` — variable type (`"continuous"`, `"binary"`, `"categorical"`)
#'   \item `class` — R class (`"factor"`, `"numeric"`, `"integer"`)
#'   \item `role` — statistical role (`"predictor"`, `"event"`, `"time"`, `"exposure"`)
#'   \item `block_table`, `block_models` — grouping for tables and models
#'   \item `include_table`, `include_uv_model`, `include_mv_model` — inclusion flags
#'   \item `forced`, `priority` — model selection controls
#' }
#'
#' `get_factor_levels()` is typically called once during data setup,
#' after factors have been assigned, to save the level structure
#' into the metadata table for downstream use.
#'
#' @examples
#' my_data<-data.frame(
#'   grade  = factor(c("low","high","medium"),levels=c("low","medium","high")),
#'   status = factor(c("alive","dead","alive")),
#'   age    = c(45,67,52)
#' )
#' meta<-data.frame(
#'   variable = c("grade","status","age"),
#'   label    = c("Tumour grade","Vital status","Age (years)"),
#'   type     = c("categorical","binary","continuous"),
#'   class    = c("factor","factor","numeric")
#' )
#' get_factor_levels(data=my_data,meta_table=meta)
#'
#' @seealso [adjust_datatable()]
#'
#' @export
get_factor_levels<-function(data,meta_table) {

  variable<-label<-type<-factor_levels<-.data<-NULL

  missing_vars<-setdiff(meta_table$variable,names(data))
  if (length(missing_vars) > 0) {
    warning("These variables are in meta_table but not in data: "
            ,paste(missing_vars,collapse=", "))
  }
  meta_table<-meta_table |>
    dplyr::mutate(factor_levels=sapply(.data$variable,function(var) {
      if (var %in% names(data) && is.factor(data[[var]])) {
        paste(levels(data[[var]]),collapse="|")
      } else {
        NA_character_
      }
    })
    ) |> dplyr::select(variable,label,type,class,factor_levels,dplyr::everything())
  return(meta_table)
}

#' Remove features not present in a reference feature set
#'
#' Filters columns from a data frame by removing entries whose corresponding
#' feature names are not found in a reference feature set. Reports which
#' features were removed. Useful for aligning datasets that may have
#' received datasets with updated features upstream. Designed to use after
#' extracting a feature-only dataset with [get_feature_dataset()].
#'
#' @param data A data frame or tibble to filter.
#' @param features A character vector of feature names corresponding to
#'   the columns of `data`. Must be present in `data`.
#' @param reference_features A character vector of feature names to keep.
#'   Features who are not in `reference_features` are removed.
#' @param verbose Logical. If `TRUE` (default), prints how many features
#'   were removed and which ones.
#'
#' @return The input `data` with columns corresponding to extra features
#'   removed. Returns `data` unchanged if no extra features are found.
#'
#' @details
#' Extra features are defined as those in `features` but not in
#' `reference_features` (i.e. `setdiff(features, reference_features)`).
#'
#' The function is particularly useful when upstream data providers
#' send updated datasets with more features than the working dataset,
#' and you need to freeze the feature set to match a reference.
#'
#' @examples
#' my_data<-data.frame("EGFR"=c(1.2, 0.8),"AKT1"=c(1.5,0.9),"TP53"=c(0.3,1.4))
#'
#' remove_extra_terms(
#'   data               = my_data,
#'   features           = colnames(my_data),
#'   reference_features = c("EGFR","AKT1")
#' )
#'
#' @seealso [get_feature_dataset()]
#'
#' @export
remove_extra_terms<-function(data,features,reference_features,verbose=TRUE) {

  #Safety checks
  if (!is.data.frame(data))
    stop("data must be a data frame or tibble")

  missing_features<-setdiff(features,colnames(data))
  if (length(missing_features) > 0) {
    stop("The following features are missing from data: ",
         paste(missing_features,collapse=", "))
  }

  extra_terms<-setdiff(features,reference_features)

  if (length(extra_terms)==0) {
    if(verbose) cat("\nNo extra terms found. Returning dataset unchanged\n")
    return(data)
  }
  if(verbose) cat("\nRemoved",length(extra_terms),"extra term(s):",
                  paste(extra_terms,collapse=", "),"\n")

  data<-data |> dplyr::select(!any_of(extra_terms))

  return(data)
}

#' Apply data types and labels from a metadata table to a dataframe
#'
#' Enforces column classes (factor, numeric, integer, Date) and variable
#' labels defined in a metadata table onto a data frame. Designed to
#' standardize a freshly imported dataset using a manually curated
#' metadata table as the source of truth. Designed to use after data cleanup
#' with [get_rawtable_values()], [safe_clean_raw_data()] and adding factor
#' levels to metatable with [get_factor_levels()].
#'
#' @param datatable A data frame or tibble to adjust. Typically imported
#'   from a CSV or XLSX file before any type enforcement.
#' @param meta_table A data frame with at least four columns: `variable`
#'   (column names in `data`), `class` (target R class),
#'   `factor_levels` (pipe-separated factor levels for factor columns,
#'   as produced by [get_factor_levels()]), and `label` (display label).
#'@param date_format A string. Determines the format of dates that will be
#'  parsed with as.Date() function.
#'
#' @return The input `data` with column classes enforced and
#'   variable labels applied via [labelled::var_label()].
#'
#' @details
#' Columns in `data` not present in `meta_table` trigger a warning
#' but are left unchanged. Columns in `meta_table` not present in
#' `data` are silently skipped.
#'
#' Supported classes: `"factor"`, `"numeric"`, `"integer"`, `"Date"`.
#' `"Date"` default behavior expects ISO 8601 format (`"%Y-%m-%d"`) and
#' other formats require explicitly passing `date_formats` argument.
#' If multiple date columns require different formats, adjust them
#' upstream so they have the same format before using this function.
#' Factor levels are read from `meta_table$factor_levels` as a
#' pipe-separated string (e.g. `"low|medium|high"`), as produced by
#' [get_factor_levels()].
#'
#' Variable labels are applied using [labelled::var_label()], making
#' them accessible to downstream reporting functions such as
#' [gtsummary::tbl_summary()].
#'
#' @examples
#' my_data<-data.frame(
#'   status = c("alive","dead","alive"),
#'   age    = c("45","67","52"),
#'   stringsAsFactors = FALSE
#' )
#' meta<-data.frame(
#'   variable      = c("status","age"),
#'   label         = c("Vital status","Age (years)"),
#'   class         = c("factor","numeric"),
#'   factor_levels = c("alive|dead",NA),
#'   stringsAsFactors = FALSE
#' )
#' adjust_datatable(data = my_data, meta_table = meta)
#'
#' @seealso [get_factor_levels()], [get_rawtable_values()],
#' [safe_clean_raw_data()] and [labelled::var_label()]
#'
#' @export
adjust_datatable<-function(datatable,meta_table,date_format="%Y-%m-%d") {

  #Check for variables in datatable not in meta_table
  missing_vars<-setdiff(names(datatable),meta_table$variable)

  if (length(missing_vars) > 0) {
    warning("These variables are in datatable but not in meta_table: ",
            paste(missing_vars,collapse=", "))
  }

  #Apply classes and factor levels
  for (i in seq_len(nrow(meta_table))) {
    var<-meta_table$variable[i]
    if (!var %in% names(datatable)) next

    if (meta_table$class[ i ]=="factor") {
      levels<-strsplit(meta_table$factor_levels[ i ],"\\|")[[1]]
      datatable[[var]]<-factor(datatable[[var]],levels=levels)
    } else if (meta_table$class[ i ]=="numeric") {
      datatable[[var]]<-as.numeric(datatable[[var]])
    } else if (meta_table$class[ i ]=="integer") {
      datatable[[var]]<-as.integer(datatable[[var]])
    } else if (meta_table$class[ i ]=="Date") {
      datatable[[var]]<-as.Date(datatable[[var]],format=date_format)
    }
  }

  #Apply labels
  labelled::var_label(datatable)<-meta_table$label[match(names(datatable),meta_table$variable)]
  return(datatable)
}

#' Convert a binary membership data frame to a named list
#'
#' Transforms any binary membership data frame into a named list, in which
#' each element is a character vector of members belonging to a specific group.
#' Works for any rows × features membership structure: pathway annotations,
#' gene sets, protein complexes, or any other custom functional grouping.
#' Designed to work with data after data cleanup via [clean_protein_names()],
#' [get_rawtable_values()] and [safe_clean_raw_data()].
#'
#' @param data A data frame or tibble with groups as rows and features
#'   as columns. Must contain the column given by `membership_col`.
#'   All values should be standardized membership strings (e.g. `"yes"`
#'   or `"no"`) after cleaning with [safe_clean_raw_data()].
#' @param membership_col A string giving the name of the column from `data`
#'   containing the membership group identifiers. This column is used
#'   to name the elements of the output list.
#' @param yes_value A string. The value indicating membership.
#'   Default is `"yes"`.
#'
#' @return A named list where each item corresponds to a group (named by
#'   `membership_col`) and contains a sorted character vector of
#'   features that belong to that membership.
#'
#' @details
#' The expected input is a binary membership reference data frame, where rows
#' represent groups and columns represent features. One of `data` columns
#' (`membership_col`) is used as list element names in the output.
#' Values not equal to `yes_value` are treated as non-membership.
#'
#' @examples
#' membership_data<-data.frame(
#'   EGFR = c("yes","yes","no"),
#'   AKT1 = c("yes","no","no"),
#'   TP53 = c("no", "no","yes"),
#'   pathways = c("PI3K","MAPK","Apoptosis"),
#'   stringsAsFactors = FALSE
#' )
#' membership_to_list(data = membership_data, membership_col = "pathways")
#'
#' @seealso [get_rawtable_values()], [safe_clean_raw_data()]
#' and [clean_protein_names()]
#'
#' @export
membership_to_list<-function(data,membership_col,yes_value="yes") {

  #Safety checks
  if (!is.data.frame(data))
    stop("data must be a data frame or tibble")
  if (!is.character(membership_col) || length(membership_col) != 1)
    stop("membership_col must be a single column name")
  if (!membership_col %in% names(data))
    stop("membership_col must be a column of data")

  membership<-data[[membership_col]]

  if (any(is.na(membership)))
    stop("membership_col contains missing values")

  membership_data<-data |> dplyr::select(-dplyr::all_of(membership_col))
  result<-lapply(seq_len(nrow(membership_data)), function( i ) {
    row<-membership_data[ i , ]
    sort(names(row)[which(unlist(row==yes_value))])
  })
  names(result)<-membership
  return(result)
}

