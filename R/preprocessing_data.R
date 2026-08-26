### Preprocessing raw data ###

#' Function to get protein dataset
#' @export
#' @examples
#' ## Placeholder Example ##
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

#' Function to get factor levels from main datatable and add it to metatable
#' @export
#' @examples
#' ## Placeholder Example ##
get_factor_levels<-function(df,meta_table) {

  missing_vars<-setdiff(meta_table$variable,names(df))
  if (length(missing_vars) > 0) {
    warning("These variables are in meta_table but not in df: "
            ,paste(missing_vars,collapse=", "))
  }
  meta_table<-meta_table %>%
    dplyr::mutate(factor_levels=sapply(variable,function(var) {
      if (var %in% names(df) && is.factor(df[[var]])) {
        paste(levels(df[[var]]),collapse="|")
      } else {
        NA_character_
      }
    })
    ) %>% dplyr::select(variable,label,type,class,factor_levels,dplyr::everything())
  return(meta_table)
}

#' Function to remove extra proteins that are not in working dataset
#' @export
#' @examples
#' ## Placeholder Example ##
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

#' Function to get raw values from raw datatable
#' @export
#' @examples
#' ## Placeholder Example ##
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

#' Function to safe clean raw data
#' @export
#' @examples
#' ## Placeholder Example ##
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

#' Function used to transform data input into a list
#' @export
#' @examples
#' ## Placeholder Example ##
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

#' Function to adjust imported Main datatable from .CSV file
#' @export
#' @examples
#' ## Placeholder Example ##
adjust_datatable<-function(datatable,meta_table) {

  #Check for variables in df not in meta_table
  missing_vars<-setdiff(names(df),meta_table$variable)

  if (length(missing_vars) > 0) {
    warning("These variables are in df but not in meta_table: ",
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
      datatable[[var]]<-as.Date(datatable[[var]])
    }
  }

  #Apply labels
  labelled::var_label(datatable)<-meta_table$label[match(names(datatable),meta_table$variable)]
  return(datatable)
}
