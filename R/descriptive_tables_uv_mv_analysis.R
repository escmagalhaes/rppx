### Functions for descriptive tables, UV and MV models ###

#' Function to detect rare and sparse variables for removal before making descriptive tables
#' @export
#' @examples
#' ## Placeholder Example ##
rare_sparse_vars<-function(df,vars,group_var,mode=c("table","outcome")
                           ,min_n=5,max_levels=6,min_freq=0.03,min_epv=10) {

  mode<-match.arg(mode)

  #Validate group_var
  if (!group_var %in% names(df)) {
    stop(paste("group_var '",group_var,"' not found in dataframe."))
  }

  #Precompute for outcome mode
  if (mode=="outcome") {
    event_vals<-df[[group_var]]
    #Safety check
    if (!all(event_vals[!is.na(event_vals)] %in% c(0,1))) {
      stop("In outcome mode, group_var must be a binary 0/1 event indicator.")
    }
    n_events<-sum(event_vals==1,na.rm=TRUE)
  }

  sparse_vars<-c()
  rare_vars<-c()
  for ( v in vars) {
    if (!v %in% names(df)) {
      warning(paste("Variable", v ,"not found in dataframe. Skipping..."))
      next
    }
    col<-df[[ v ]]
    col_vals<-col[!is.na(col)]
    n_levels<-length(unique(col_vals))

    #Rare frequency filter (categorical/binary variables only)#
    if (n_levels == 2) {
      tab <- table(col, df[[group_var]], useNA = "no")
      if (ncol(tab) == 0 || nrow(tab) == 0) { next }
      col_sums <- colSums(tab)
      if (any(col_sums == 0)) { sparse_vars <- c(sparse_vars, v); next }
      group_freqs <- sweep(tab, 2, col_sums, "/")
      max_freq_per_level <- apply(group_freqs, 1, max)
      if (all(max_freq_per_level < min_freq)) {
        rare_vars <- c(rare_vars, v)
        next
      }
    }

    #Table mode: sparsity filter#
    if (mode == "table") {
      tab <- table(col, df[[group_var]])
      level_totals <- rowSums(tab)
      col_totals <- colSums(tab)
      n_tab_levels <- nrow(tab)

      if (n_tab_levels < 2) { sparse_vars <- c(sparse_vars, v); next }
      if (n_tab_levels > max_levels) { sparse_vars <- c(sparse_vars, v); next }

      # Skip if entire variable is absent in one group
      if (any(col_totals == 0)) { sparse_vars <- c(sparse_vars, v); next }

      # Skip min_n only if no perfect separation (0 in one group is informative)
      col_mins <- apply(tab, 1, min)
      if (!any(col_mins == 0) && any(level_totals < min_n)) {
        sparse_vars <- c(sparse_vars, v); next
      }
    }

    #Outcome mode: EPV filter (categorical variables only)
    if (mode=="outcome") {
      if (is.factor(col) || is.character(col)) {
        tab_events<-table(col,df[[group_var]])
        if (ncol(tab_events) < 2 || !"1" %in% colnames(tab_events)) {
          sparse_vars<-c(sparse_vars, v )
          next
        }
        event_totals<-tab_events[,"1"]
        if (any(event_totals < min_epv)) {
          sparse_vars<-c(sparse_vars, v )
          next
        }
      }
    }
  }
  return(list(sparse_vars=setdiff(unique(sparse_vars),group_var)
              ,rare_vars=setdiff(unique(rare_vars),group_var)))
}

#' Function to detect low sample size and indicate test
#' @export
#' @examples
#' ## Placeholder Example ##
define_stat_test<-function(df,grouping_variable,vars) {

  test_list<-list()
  test_args<-list()
  n_groups<-nlevels(factor(df[[grouping_variable]]))

  for (v in vars) {
    #Categorical variables
    if (is.factor(df[[ v ]]) || is.character(df[[ v ]])) {
      tab<-table(df[[ v ]], df[[grouping_variable]])
      expected<-suppressWarnings(chisq.test(tab)$expected)

      if (all(expected >= 5)) {
        test_list[[ v ]]<-"chisq.test"
      } else {
        test_list[[ v ]]<-"fisher.test"
        if (nrow(tab) > 2 || ncol(tab) > 2) {
          test_args[[ v ]]<-list(simulate.p.value=TRUE,B=10000)
        }
      }

      #Continuous variables
    } else if (is.numeric(df[[ v ]])) {

      #Normality check per group using Shapiro-Wilk
      groups<-split(df[[ v ]],df[[grouping_variable]])
      groups<-lapply(groups,function(x) x[!is.na(x)])

      #Shapiro-Wilk requires at least 3 observations per group
      normal<-all(sapply(groups,function(x) {
        if (length(x) < 3)  return(FALSE)
        if (length(x) > 5000) {
          warning(paste("Variable", v , ": group has more than 5000 observations, normality assumed."))
          return(TRUE)
        }
        if (var( x,na.rm=TRUE)==0) return(FALSE)
        shapiro.test( x )$p.value > 0.05
      }))

      test_list[[ v ]]<-if (normal) {
        ifelse(n_groups==2,"t.test","aov")
      } else {
        ifelse(n_groups==2,"wilcox.test","kruskal.test")
      }

      #Add exact=FALSE for wilcox.test to avoid ties warning
      if (!normal && n_groups==2) {
        test_args[[ v ]]<-list(exact=FALSE)
      }
    }
  }

  #Adjust 'test' formating for gtsummary handling
  tests_fmt<-setNames(lapply(names(test_list), function( v ) {
      reformulate(paste0('"',test_list[[ v ]],'"'),response = v )
    }),names(test_list) )

  return(list(tests=tests_fmt,test_args=test_args))
}

#' Function to make descriptive table with gtsummary
#' @export
#' @examples
#' ## Placeholder Example ##
get_descr_tab<-function(df,grouping_variable,vars,skip_pval=-9999,vars_to_test=NULL
                        ,tests=NULL,args=NULL,skipped_vars=NULL,exclude_sig=NULL
                        ,sig_type=NULL,sig_value=NULL,sig_labels=NULL
                        ,seed=123){

  #Helper function to adapt test list to gtsummary
  filter_cat_tests<-function(tests,df,vars) {
    tests[names(tests) %in% names(
      which(sapply(df[vars],function(x) is.factor(x) || is.character(x))))]
  }

  set.seed(seed) #For fisher's exact test with simulated p-values
  table<-df |> dplyr::select(all_of(c(grouping_variable,vars))) |>
    tbl_summary(by=grouping_variable
                ,statistic=list(all_continuous()~"{mean}({sd})"
                                ,all_categorical()~"{n}/{N}({p}%)")
                ,percent='column'
                ,digits=all_continuous()~1
                ,type=all_categorical()~"categorical"
                ,missing='ifany'
                ,missing_text='N/A')

  if(!is.null(vars_to_test) && !is.null(tests)){
    cat_vars_to_test <- vars_to_test[sapply(vars_to_test, function(v)
      is.factor(df[[v]]) || is.character(df[[v]]))]

    table<-table |> add_p(include=all_of(vars_to_test)
                           ,test=unname(tests[names(tests) %in% cat_vars_to_test])
                           ,test.args=args) |>
      modify_table_body(~.x |>
                          dplyr::mutate(p.value=case_when(
                            !is.null(skipped_vars)
                            & variable %in% skipped_vars
                            & row_type=="label"~skip_pval
                            ,TRUE~p.value))) |>
      modify_fmt_fun(p.value=function(x) {
        dplyr::case_when(x==skip_pval~"-"
                         ,is.na(x)~""
                         ,TRUE~ style_pvalue(x))})
  } else {
    table<-table |> add_p()
  }
  table<-table |> bold_p() |> add_overall() |> bold_labels()

  #Extract significant variables
  sig_vars<-table$table_body |> dplyr::filter(row_type=="label"
                                        & p.value < 0.05
                                        & p.value !=skip_pval) |> pull(variable)

  #Check if sig_vars is not NULL
  if (length(sig_vars)==0) {
    message("No significant variables found")
    return(list(extended_table=table,significant_table=NULL))
  }

  #Remove excluded variables
  sig_vars<-sig_vars[!sig_vars %in% exclude_sig]

  #Safety check in case all sig_vars may have been excluded
  if (length(sig_vars)==0) {
    message("All significant variables were excluded or not significant variables were detected. No significant table produced.")
    return(list(extended_table=table,significant_table=NULL))
  }

  #Generate table with significant variables only
  cat_sig_vars <- sig_vars[sapply(sig_vars, function(v)
    is.factor(df[[v]]) || is.character(df[[v]]))]

  set.seed(seed) #For fisher's exact test with simulated p-values
  sig_table<-suppressWarnings(
    df |> dplyr::select(all_of(c(grouping_variable,sig_vars))) |>
    tbl_summary(by=grouping_variable
                ,statistic=list(all_continuous()~"{mean}",all_categorical()~"{p}%")
                ,percent='column'
                ,digits=all_continuous()~1
                ,type=sig_type
                ,value=sig_value
                ,label=sig_labels
                ,missing='no'
                ,missing_text='N/A') |>
    add_p(include=all_of(sig_vars)
          ,test=unname(tests[names(tests) %in% cat_sig_vars])
          ,test.args=args[names(args) %in% sig_vars]
          ,pvalue_fun=~style_pvalue(.x,digits=2)) |> bold_p() |> bold_labels()
  )

  return(list(extended_table=table,significant_table=sig_table))
}

#' Function to merge descriptive tables of significant variables
#' @export
#' @examples
#' ## Placeholder Example ##
get_comb_sig_tables<-function(table_list,grouping_variable="group",exclude_sig=NULL
                              ,label_vars=NULL,dichotomous_vars=NULL,seed=123) {

  if (!is.list(table_list) || length(table_list) < 2) {
    stop("table_list must be a list of at least 2 tables.")
  }

  #Get union of significant variables across both pipelines
  sig_vars_combined<-unique(unlist(lapply(table_list,function(tab) {
    if (is.null(tab$descr_table_sig)) return(character(0))
    tab$descr_table_sig$table_body |> dplyr::filter(row_type=="label") |> dplyr::pull(variable)
  })))

  if (length(sig_vars_combined)==0) {
    message("No significant variables found across any table. Returning NULL.")
    return(NULL)
  }

  #Remove excluded variables
  if (!is.null(exclude_sig)) { sig_vars_combined<-setdiff(sig_vars_combined,exclude_sig) }

  #Reorder according to meta_table order (union across all tables)
  meta_order<-Reduce(union,lapply(table_list,function(tab) tab$meta_table$variable))
  sig_vars_combined<-meta_order[meta_order %in% sig_vars_combined]

  #Helper to rebuild sig table
  rebuild_sig<-function(tab) {
    df<-tab$processed_df
    meta<-tab$meta_table

    #Build type and value lists
    cat_vars<-meta |> dplyr::filter(type %in% c("binary","categorical"),include_table==TRUE,
             !role %in% c("exposure","event","time")) |> pull(variable)

    binary_tab_vars<-cat_vars[sapply(cat_vars, function( v ) {
      if (!v %in% names(df)) return(FALSE)
      levs <- unique(df[[ v ]][!is.na(df[[ v ]])])
      length(levs)==2 && all(tolower(as.character(levs)) %in% c("yes","no"))
    })]

    type_tab_list <-setNames(lapply(binary_tab_vars,function(v) "dichotomous"),binary_tab_vars)
    value_tab_list<-setNames(lapply(binary_tab_vars,function(v) "yes"),binary_tab_vars)

    if (!is.null(dichotomous_vars)) {
      for ( v in names(dichotomous_vars)) {
        if ( v %in% names(df)) {
          type_tab_list[[ v ]]<-"dichotomous"
          value_tab_list[[ v ]]<-dichotomous_vars[[ v ]]
        }
      }
    }

    #Filter to variables present in this df
    vars_present<-sig_vars_combined[sig_vars_combined %in% names(df)]

    #Categorical vars among present sig vars
    cat_present<-vars_present[sapply(vars_present, function( v )
      is.factor(df[[ v ]]) || is.character(df[[ v ]]))]

    #Rebuild stat tests
    stat_tests<-define_stat_test(df=df,grouping_variable=grouping_variable,vars=cat_present)

    #Rebuild sig table
    set.seed(seed) #For fisher's exact test with simulated p-values
    df |> dplyr::select(all_of(c(grouping_variable,vars_present))) |>
      tbl_summary(by=grouping_variable
                  ,statistic=list(all_continuous()~"{mean}",all_categorical()~"{p}%")
                  ,percent="column"
                  ,digits=all_continuous()~1
                  ,type=type_tab_list[names(type_tab_list) %in% vars_present]
                  ,value=value_tab_list[names(value_tab_list) %in% vars_present]
                  ,label=label_vars
                  ,missing="no") |>
      add_p(include=all_of(vars_present)
            ,test=unname(stat_tests$tests[names(stat_tests$tests) %in% cat_present])
            ,test.args=stat_tests$test_args[names(stat_tests$test_args) %in% vars_present]
            ,pvalue_fun=~ style_pvalue(.x,digits=2)) |>
      bold_p() |> bold_labels()
  }

  #Return tables as a named list
  sig_tables<-lapply(table_list,rebuild_sig)
  names(sig_tables)<- if(!is.null(names(table_list))) {
    paste0("sig_table_",names(table_list))
  } else {
    paste0("sig_table",seq_along(sig_tables))
  }
  return(c(sig_tables,list(sig_vars_used=sig_vars_combined)))
}

#' Function to compute UV models and extract best variables for MVA
#' @export
#' @examples
#' ## Placeholder Example ##
run_uv_coxph<-function(df,vars,time_var,event_var,uv_p_cutoff=0.20
                       ,labels=NULL,meta_table=NULL,group_var_label=NULL){

  #Temporarily assign to global environment to avoid issue with tbl_uvregression
  .uv_data<<-df
  .uv_surv<<-Surv(df[[time_var]],df[[event_var]])
  on.exit({ rm(list=c(".uv_data",".uv_surv"),envir=.GlobalEnv) },add=TRUE)

  #Detect binary vars for show_single_row
  binary_vars<-unique(c("group",vars))[sapply(unique(c("group",vars)),function( v ) {
    levs<-unique(df[[ v ]][!is.na(df[[ v ]])])
    length(levs)==2
  })]

  #Full UV table
  uv_table<-tbl_uvregression(data=.uv_data
                             ,y=.uv_surv
                             ,include=all_of(unique(c("group",vars)))
                             ,method=coxph,exponentiate=TRUE,hide_n=TRUE
                             ,add_estimate_to_reference_rows=TRUE
                             ,pvalue_fun=~style_pvalue(.x,digits=2)) |>
    bold_p() |> bold_labels()

  #Extract significant variables
  sig_uv_vars<-uv_table$table_body |>
    dplyr::filter(variable !="group",
           (row_type=="label" & !is.na(p.value) & p.value < 0.05)
      | (row_type=="level" & !is.na(p.value) & p.value < 0.05)) |>
    dplyr::pull(variable) |> unique()

  #Check if sig_uv_vars is not NULL
  if (length(sig_uv_vars)==1) {
    message("No significant variables found. Sig table will show group only.")
    return(list(uv_table=uv_table,sig_uv_table=NULL,selected_mv_vars=NULL))
  }

  #Grouping variable must appear always in sig table regardless of significance
  sig_uv_vars<-unique(c("group",sig_uv_vars))

  #Binary vars in sig table
  binary_sig_vars<-setdiff(intersect(binary_vars,sig_uv_vars),"group")

  #Create label list for significant table
  sig_labels<-setNames(as.list(meta_table$label[meta_table$variable %in% sig_uv_vars]),
    meta_table$variable[meta_table$variable %in% sig_uv_vars])
  if (!is.null(labels)) { for ( v in names(labels)) { sig_labels[[ v ]]<-labels[[ v ]] } }

  #Add group label if present
  if (!is.null(group_var_label)) { sig_labels[["group"]]<-group_var_label }

  # Update global surv object for sig table
  .uv_surv<<-Surv(df[[time_var]],df[[event_var]])

  #Significant only UV model (refit needed for show_single_row display)
  sig_uv_table<-tbl_uvregression(data=.uv_data
                                 ,y=.uv_surv
                                 ,include=all_of(sig_uv_vars)
                                 ,show_single_row=all_of(binary_sig_vars)
                                 ,method=coxph,exponentiate=TRUE,hide_n=TRUE
                                 ,add_estimate_to_reference_rows=TRUE
                                 ,label=if(length(sig_labels)>0) sig_labels else NULL
                                 ,pvalue_fun= ~style_pvalue(.x,digits=2)) |>
    bold_p() |> bold_labels()

  #Extract variables for MV model
  sel_mv_vars<-uv_table$table_body |>
    dplyr::filter(variable != "group",
           (row_type=="label" & !is.na(p.value) & p.value < uv_p_cutoff)
           | (row_type=="level" & !is.na(p.value) & p.value < uv_p_cutoff)) |>
    dplyr::pull(variable) |> unique()

  return(list(uv_table=uv_table,sig_uv_table=sig_uv_table,selected_mv_vars=sel_mv_vars))
}

#' Function to perform forward selection in MV models
#' @export
#' @examples
#' ## Placeholder Example ##
mv_variable_selection<- function(data,candidate_vars,mv_vars_pool,
                                 time_var,event_var,grouping_variable,
                                 max_predictors, meta_table,min_epv=10,
                                 cluster_var_p_cutoff=0.05,
                                 lrt_p_cutoff=0.05,block_vars=TRUE,
                                 selection_strategy=c("strict","two_phase","cie"),
                                 cie_threshold=0.10,cie_final_cleanup=FALSE,
                                 all_levels_sig=FALSE,
                                 forced_selected=NULL,verbose=TRUE) {

  selection_strategy<-match.arg(selection_strategy)

  #Helper functions#
  # Helper: get minimum group p-value from a fitted model
  get_group_p <- function(model) {
    broom::tidy(model) |>
      dplyr::filter(startsWith(term, grouping_variable)) |>
      dplyr::pull(p.value) |>
      min(na.rm = TRUE)
  }

  # Helper: get group HR from a fitted model
  get_group_hr <- function(model) {
    broom::tidy(model) |>
      dplyr::filter(startsWith(term, grouping_variable)) |>
      dplyr::pull(estimate) |>
      mean(na.rm = TRUE)
  }

  # Helper: apply blocking/exclusion after variable is selected
  apply_blocks <- function(v, mv_vars_local) {
    v_meta         <- meta_table |> dplyr::filter(variable == v)
    v_excl_group   <- v_meta |> dplyr::pull(exclusive_group)
    v_blocks_group <- if ("blocks_group" %in% names(meta_table)) {
      v_meta |> pull(blocks_group)
    } else NA

    related_excl <- if (!is.na(v_excl_group)) {
      setdiff(
        meta_table |> dplyr::filter(exclusive_group == v_excl_group) |>
          dplyr::pull(variable), v
      )
    } else c()

    related_blocked <- if (!is.na(v_blocks_group)) {
      setdiff(
        meta_table |> dplyr::filter(exclusive_group == v_blocks_group) |>
          dplyr::pull(variable), v
      )
    } else c()

    all_related <- unique(c(related_excl, related_blocked))

    if (length(all_related) > 0) {
      if (isTRUE(block_vars)) {
        mv_vars_local <- setdiff(mv_vars_local, all_related)
        if (verbose) cat("  Excluding (blocked by '", v, "'): ",
                         paste(all_related, collapse = ", "), "\n", sep = "")
      } else {
        if (verbose) cat("  Deprioritizing (soft block by '", v, "'): ",
                         paste(all_related, collapse = ", "), "\n", sep = "")
      }
    }
    mv_vars_local
  }

  #Initialize
  selected      <- if (!is.null(forced_selected)) forced_selected else c()
  mv_vars_local <- mv_vars_pool
  continue      <- TRUE
  current_model <- coxph(
    as.formula(paste0("Surv(", time_var, ",", event_var, ") ~ ", grouping_variable)),
    data = data
  )

  # Build starting model: includes forced variables if any
  if (length(selected) > 0) {
    init_formula <- as.formula(
      paste0("Surv(", time_var, ",", event_var, ") ~ ",
             paste(c(grouping_variable, selected), collapse = " + "))
    )
    current_model <- tryCatch(
      coxph(init_formula, data = data),
      error = function(e) {
        warning("Could not fit initial model with forced variables. ",
                "Starting from base model.")
        coxph(as.formula(
          paste0("Surv(", time_var, ",", event_var, ") ~ ", grouping_variable)
        ), data = data)
      }
    )
    if (verbose) {
      cat("\nInitial model with forced variables: group p =",
          round(get_group_p(current_model), 3), "\n")
    }
  } else {
    #No forced variables: start from base model
    current_model <- coxph(
      as.formula(paste0("Surv(", time_var, ",", event_var,
                        ") ~ ", grouping_variable)),
      data = data
    )
  }

  #Phase 1: Forward selection
  if (verbose) cat("\nPhase 1: Forward selection (strategy: '",
                   selection_strategy, "')\n", sep = "")

  while (continue && length(selected) < max_predictors) {
    added <- FALSE

    # Sort candidates: deprioritize exclusive/blocked groups
    selected_excl_groups <- meta_table |>
      dplyr::filter(variable %in% selected, !is.na(exclusive_group)) |>
      dplyr::pull(exclusive_group) |> unique()

    selected_blocks_groups <- if ("blocks_group" %in% names(meta_table)) {
      meta_table |>
        dplyr::filter(variable %in% selected, !is.na(blocks_group)) |>
        dplyr::pull(blocks_group) |> unique()
    } else c()

    all_deprioritized_groups <- if (isTRUE(block_vars)) {
      selected_excl_groups
    } else {
      unique(c(selected_excl_groups, selected_blocks_groups))
    }

    sorted_candidates <- candidate_vars |>
      dplyr::filter(variable %in% mv_vars_local) |>
      mutate(group_overlap = ifelse(
        !is.na(exclusive_group) & exclusive_group %in% all_deprioritized_groups, 1, 0
      )) |>
      dplyr::arrange(group_overlap, mv_order) |>
      dplyr::pull(variable)

    if (verbose) cat("\nTesting candidates in order:",
                     paste(sorted_candidates, collapse = ", "), "\n")

    for (v in sorted_candidates) {

      # EPV check
      if (is.factor(data[[v]])) {
        tab_events <- table(data[[v]], data[[event_var]])
        if (ncol(tab_events) < 2)             next
        if (any(tab_events[, "1"] < min_epv)) next
      }

      # Complete case subset for this candidate
      complete_v      <- stats::complete.cases(
        data[, intersect(c(time_var, event_var, grouping_variable, selected, v),
                         names(data))]
      )
      data_v          <- data[complete_v, ]

      if (verbose) cat("  Testing:", v,
                       "| n complete:", sum(complete_v),
                       "| current model n:", current_model$n, "\n")

      # Refit current model on same subset
      current_formula_v <- as.formula(
        paste0("Surv(", time_var, ",", event_var, ") ~ ",
               paste(c(grouping_variable, selected), collapse = " + "))
      )
      current_model_v <- tryCatch(
        coxph(current_formula_v, data = data_v),
        error = function(e) NULL
      )
      if (is.null(current_model_v)) next

      # Fit candidate model on same subset
      new_vars    <- c(grouping_variable, selected, v)
      new_formula <- as.formula(
        paste0("Surv(", time_var, ",", event_var, ") ~ ",
               paste(new_vars, collapse = " + "))
      )
      new_model <- tryCatch(coxph(new_formula, data = data_v), error = function(e) NULL)
      if (is.null(new_model)) next

      #Retention criterion depends on strategy
      retain <- FALSE

      if (selection_strategy == "strict") {

        #Gate 1: LRT = variable improves model fit
        lrt   <- tryCatch(
          anova(current_model_v, new_model, test = "LRT"),
          error = function(e) {
            if (verbose) cat("  LRT error for", v, ":", conditionMessage(e), "\n")
            NULL
          }
        )
        p_lrt <- if (!is.null(lrt)) lrt$"Pr(>|Chi|)"[2] else NA

        if (verbose) cat("  Testing:", v, "| LRT p:",
                         ifelse(!is.na(p_lrt), round(p_lrt, 4), "NA"), "\n")

        if (is.na(p_lrt) || p_lrt >= lrt_p_cutoff) next

        #Gate 2: Wald = variable individually significant
        var_wald_p <- broom::tidy(new_model) |>
          dplyr::filter(startsWith(term, v)) |>
          dplyr::pull(p.value) |>
          min(na.rm = TRUE)

        if (verbose) cat("  Wald p for '", v, "': ", round(var_wald_p, 4), "\n", sep = "")

        if (is.na(var_wald_p) || var_wald_p >= lrt_p_cutoff) next

        #Gate 3: group remains globally significant
        formula_no_cluster <- as.formula(
          paste0("Surv(", time_var, ",", event_var, ") ~ ",
                 paste(c(selected, v), collapse = " + "))
        )
        fit_no_cluster   <- tryCatch(
          coxph(formula_no_cluster, data = data_v),
          error = function(e) NULL
        )
        cluster_global_p <- if (!is.null(fit_no_cluster)) {
          lrt_cl <- tryCatch(
            anova(fit_no_cluster, new_model, test = "LRT"),
            error = function(e) NULL
          )
          if (!is.null(lrt_cl)) lrt_cl$"Pr(>|Chi|)"[2] else NA
        } else NA

        cluster_sig <- !is.na(cluster_global_p) && cluster_global_p < cluster_var_p_cutoff
        if (!cluster_sig) {
          if (verbose) cat("\n\nGroup global p:",
                           ifelse(!is.na(cluster_global_p), round(cluster_global_p, 4), "NA"),
                           "| Group not significant: skipping\n")
          next
        }

        # Gate 4 (optional): all cluster levels individually significant
        if (isTRUE(all_levels_sig)) {
          level_ps <- broom::tidy(new_model) |>
            dplyr::filter(startsWith(term, grouping_variable)) |>
            dplyr::select(term, p.value)

          all_sig <- all(level_ps$p.value < cluster_var_p_cutoff, na.rm = TRUE)

          if (verbose) {
            cat("  Cluster level p-values:\n")
            for (i in seq_len(nrow(level_ps))) {
              cat("    ", level_ps$term[i], ":", round(level_ps$p.value[i], 4),
                  ifelse(level_ps$p.value[i] < cluster_var_p_cutoff, "OK", "NOT OK"), "\n")
            }
          }

          if (!all_sig) {
            if (verbose) cat("\n\nNot all cluster levels significant: skipping '",
                             v, "'\n", sep = "")
            next
          }
        }

        retain <- TRUE

        if (verbose) cat("  --> Group global p:", round(cluster_global_p, 4),
                         "| All gates passed | Retain: TRUE\n")
      } else if (selection_strategy == "two_phase") {
        lrt   <- tryCatch(
          anova(current_model_v, new_model, test = "LRT"),
          error = function(e) NULL
        )
        p_lrt  <- if (!is.null(lrt)) lrt$"Pr(>|Chi|)"[2] else NA
        retain <- !is.na(p_lrt) && p_lrt < lrt_p_cutoff

        if (verbose) cat("  Testing:", v, "| LRT p:",
                         ifelse(!is.na(p_lrt), round(p_lrt, 4), "NA"),
                         "| Retain:", retain, "\n")

      } else if (selection_strategy == "cie") {
        hr_current <- get_group_hr(current_model_v)
        hr_new     <- get_group_hr(new_model)
        pct_change <- abs(exp(hr_new) - exp(hr_current)) / exp(hr_current)
        retain     <- pct_change >= cie_threshold

        if (verbose) cat("  Testing:", v, "| HR change:",
                         round(pct_change * 100, 1), "%",
                         "| Retain:", retain, "\n")
      }

      if (retain) {
        selected      <- c(selected, v)
        mv_vars_local <- setdiff(mv_vars_local, v)
        mv_vars_local <- apply_blocks(v, mv_vars_local)
        # Update current_model on FULL data
        current_model <- tryCatch(
          coxph(as.formula(
            paste0("Surv(", time_var, ",", event_var, ") ~ ",
                   paste(c(grouping_variable, selected), collapse = " + "))
          ), data = data),
          error = function(e) current_model
        )
        added <- TRUE
        break
      }
    }
    if (!added) continue <- FALSE
  }

  #Phase 2: Backward pruning (two_phase only)
  if (selection_strategy == "two_phase" && length(selected) > 0) {

    group_p <- get_group_p(current_model)

    # Check all stopping criteria
    group_sig  <- group_p < cluster_var_p_cutoff

    all_vars_sig <- if (isTRUE(all_levels_sig)) {
      all(broom::tidy(current_model) |>
            dplyr::filter(!startsWith(term, grouping_variable)) |>
            dplyr::pull(p.value) < lrt_p_cutoff, na.rm = TRUE)
    } else TRUE

    all_cluster_levels_sig <- if (isTRUE(all_levels_sig)) {
      all(broom::tidy(current_model) |>
            dplyr::filter(startsWith(term, grouping_variable)) |>
            dplyr::pull(p.value) < cluster_var_p_cutoff, na.rm = TRUE)
    } else TRUE

    if (group_sig && all_vars_sig && all_cluster_levels_sig) {
      if (verbose) cat("\nAll criteria met after Phase 1. No pruning needed.",
                       "\n  Group p =", round(group_p, 3), "\n")
    } else {

      if (verbose) {
        cat("\nPhase 2: Backward pruning\n")
        cat("  Group p =", round(group_p, 3),
            "| Group significant:", group_sig, "\n")
        if (isTRUE(all_levels_sig)) {
          cat("  All vars significant:", all_vars_sig, "\n")
          cat("  All cluster levels significant:", all_cluster_levels_sig, "\n")
        }
      }

      pruned         <- selected
      pruned_model   <- current_model

      # Continue pruning until all criteria met or no variables left
      while (length(pruned) > 0) {

        # Check stopping criteria on current pruned model
        group_p_pruned <- get_group_p(pruned_model)
        group_sig_now  <- group_p_pruned < cluster_var_p_cutoff

        all_vars_sig_now <- if (isTRUE(all_levels_sig)) {
          # Check non-forced variables only
          non_forced_terms <- broom::tidy(pruned_model) |>
            dplyr::filter(!startsWith(term, grouping_variable)) |>
            dplyr::filter(if (!is.null(forced_selected)) {
              !sapply(term, function(t)
                any(sapply(forced_selected, function(v) startsWith(t, v))))
            } else TRUE)

          nrow(non_forced_terms) == 0 ||
            all(non_forced_terms$p.value < lrt_p_cutoff, na.rm = TRUE)
        } else TRUE

        all_cluster_levels_sig_now <- if (isTRUE(all_levels_sig)) {
          all(broom::tidy(pruned_model) |>
                dplyr::filter(startsWith(term, grouping_variable)) |>
                dplyr::pull(p.value) < cluster_var_p_cutoff, na.rm = TRUE)
        } else TRUE

        # Stop if all criteria met
        if (group_sig_now && all_vars_sig_now && all_cluster_levels_sig_now) {
          if (verbose) cat("\nAll criteria met after pruning.",
                           "\n  Group p =", round(group_p_pruned, 3),
                           "| Predictors remaining:", length(pruned), "\n")
          break
        }

        # Find worst variable to remove
        # Priority: non-significant variables first, then least significant
        model_terms <- broom::tidy(pruned_model) |>
          dplyr::filter(!startsWith(term, grouping_variable))

        # Never remove forced variables
        if (!is.null(forced_selected)) {
          model_terms <- model_terms |>
            dplyr::filter(!sapply(term, function(t)
              any(sapply(forced_selected, function(v) startsWith(t, v)))))
        }

        if (nrow(model_terms) == 0) break

        worst_term <- model_terms |>
          dplyr::arrange(desc(p.value)) |>
          dplyr::slice(1)

        # Map term to variable name
        non_forced <- if (!is.null(forced_selected)) {
          setdiff(pruned, forced_selected)
        } else pruned

        worst_var <- non_forced[sapply(non_forced, function(v)
          startsWith(worst_term$term, v))]

        if (length(worst_var) == 0) {
          if (verbose) cat("  Could not map term '", worst_term$term,
                           "' to variable. Stopping pruning.\n", sep = "")
          break
        }

        if (verbose) cat("  Removing '", worst_var,
                         "' (p = ", round(worst_term$p.value, 4), ")\n", sep = "")

        pruned <- setdiff(pruned, worst_var)

        # Refit pruned model
        if (length(pruned) == 0) {
          pruned_model <- coxph(
            as.formula(paste0("Surv(", time_var, ",", event_var,
                              ") ~ ", grouping_variable)),
            data = data
          )
        } else {
          pruned_formula <- as.formula(
            paste0("Surv(", time_var, ",", event_var, ") ~ ",
                   paste(c(grouping_variable, pruned), collapse = " + "))
          )
          pruned_model <- tryCatch(
            coxph(pruned_formula, data = data),
            error = function(e) NULL
          )
          if (is.null(pruned_model)) break
        }
      }

      # Final check after pruning
      final_group_p <- get_group_p(pruned_model)

      if (final_group_p < cluster_var_p_cutoff) {
        if (verbose) cat("\nPruning complete. Final model has",
                         length(pruned), "predictor(s)",
                         "| Group p =", round(final_group_p, 3), "\n")
        selected      <- pruned
        current_model <- pruned_model
      } else {
        if (verbose) cat("\nWarning: group remains non-significant (p =",
                         round(final_group_p, 3),
                         ") even after complete backward pruning.",
                         "Returning base model.\n")
        selected      <- if (!is.null(forced_selected)) forced_selected else c()
        current_model <- coxph(
          as.formula(paste0("Surv(", time_var, ",", event_var,
                            ") ~ ", grouping_variable,
                            if (!is.null(forced_selected))
                              paste0(" + ", paste(forced_selected, collapse = " + "))
                            else "")),
          data = data
        )
      }
    }
  }

  #CIE cleanup (cie strategy + cie_final_cleanup = TRUE only)
  if (selection_strategy == "cie" && isTRUE(cie_final_cleanup) && length(selected) > 0) {

    if (verbose) cat("\nCIE cleanup: removing non-significant variables",
                     "that don't meaningfully affect group HR...\n")

    cleanup <- TRUE
    base_hr <- get_group_hr(current_model)

    while (cleanup && length(selected) > 0) {
      cleanup <- FALSE

      # Find non-significant variable with highest p-value
      nonsig_term <- broom::tidy(current_model) |>
        dplyr::filter(!startsWith(term, grouping_variable), p.value >= 0.05) |>
        dplyr::arrange(desc(p.value)) |>
        dplyr::slice(1) |>
        dplyr::pull(term)

      if (length(nonsig_term) == 0 || is.na(nonsig_term)) break

      # Map term back to variable name
      worst_var <- selected[sapply(selected, function(v) startsWith(nonsig_term, v))]
      if (length(worst_var) == 0) break

      # Try removing it
      trial <- setdiff(selected, worst_var)

      if (length(trial) == 0) {
        trial_model <- coxph(
          as.formula(paste0("Surv(", time_var, ",", event_var,
                            ") ~ ", grouping_variable)),
          data = data
        )
      } else {
        trial_formula <- as.formula(
          paste0("Surv(", time_var, ",", event_var, ") ~ ",
                 paste(c(grouping_variable, trial), collapse = " + "))
        )
        trial_model <- tryCatch(
          coxph(trial_formula, data = data),
          error = function(e) NULL
        )
      }
      if (is.null(trial_model)) break

      # Check if group HR changes by less than cie_threshold after removal
      trial_hr   <- get_group_hr(trial_model)
      pct_change <- abs(exp(trial_hr) - exp(base_hr)) / exp(base_hr)

      worst_p <- broom::tidy(current_model) |>
        dplyr::filter(term == nonsig_term) |>
        dplyr::pull(p.value)

      if (pct_change < cie_threshold) {
        #Safe to remove: not a true confounder
        if (verbose) cat("  Removing '", worst_var,
                         "' (p = ", round(worst_p, 3),
                         ", HR change = ", round(pct_change * 100, 1),
                         "%: below threshold)\n", sep = "")
        selected      <- trial
        current_model <- trial_model
        base_hr       <- trial_hr
        cleanup       <- TRUE
      } else {
        #Keep: true confounder despite non-significance
        if (verbose) cat("  Keeping '", worst_var,
                         "' (p = ", round(worst_p, 3),
                         ", HR change = ", round(pct_change * 100, 1),
                         "%: true confounder)\n", sep = "")
      }
    }

    if (verbose) cat("\nCIE cleanup complete. Final model has",
                     length(selected), "predictor(s).\n")
  }

  #Final group significance note
  if (length(selected) > 0) {
    final_group_p  <- get_group_p(current_model)
    final_group_hr <- round(exp(get_group_hr(current_model)), 3)
    if (verbose) cat("\nFinal model: ", length(selected), "predictor(s) | ",
                     "group HR =", final_group_hr,
                     "| group p =", round(final_group_p, 3), "\n")
  }

  return(list(selected = selected, final_model = current_model))
}

#' Function to generate stepwise forward MV models
#' @export
#' @examples
#' ## Placeholder Example ##
run_mv_coxph<-function(df,meta_table,vars,time_var,event_var,grouping_variable
                       ,priority_levels=c("core","optional","low")
                       ,selection_strategy=c("strict","two_phase","cie")
                       ,cie_threshold=0.10,cie_final_cleanup=FALSE
                       ,priority_override=NULL,min_epv=10
                       ,max_vars=NULL,cluster_var_p_cutoff=0.05,labels=NULL
                       ,max_missing_mv=0.20,warn_missing_n=50,lrt_p_cutoff=0.05
                       ,block_vars=TRUE,force_mv_vars=NULL,all_levels_sig=FALSE
                       ,run_bootstrap=FALSE,n_bootstrap=200
                       ,bootstrap_cutoff=0.50,parallel=FALSE,seed=123
                       ,n_cores=parallel::detectCores()-1,verbose=TRUE) {

  selection_strategy<-match.arg(selection_strategy)

  ## Initial setup ##
  #Generate a variable candidate list from meta_table
  candidate_vars<-meta_table |>
    dplyr::filter(variable %in% vars,
           include_mv_model==TRUE,
           !is.na(priority),
           priority %in% priority_levels) |> arrange(mv_order)

  #Apply priority override
  if (!is.null(priority_override)) {
    for (v in names(priority_override)) {
      if (v %in% candidate_vars$variable) {
        candidate_vars<-candidate_vars |>
          dplyr::mutate(
            priority=ifelse(variable== v ,priority_override[[ v ]],priority))
        if(verbose) cat("  Priority override applied: '", v, "' -> '"
                        ,priority_override[[ v ]],"'\n",sep="")
      } else {
        warning(paste("priority_override: variable '", v
                      ,"' not found in candidate_vars. Skipping."))
      }
    }
    # Re-filter by priority_levels after override
    candidate_vars<-candidate_vars |>  dplyr::filter(priority %in% priority_levels)
  }

  #Filter Missingness
  #Compute missingness for each candidate variable
  missing_rates<-sapply(candidate_vars$variable, function( v ) { mean(is.na(df[[ v ]])) })
  excluded_missing<-candidate_vars$variable[missing_rates > max_missing_mv]

  if (length(excluded_missing) > 0) {
    excluded_labels<-meta_table$label[meta_table$variable %in% excluded_missing]
    if(verbose) cat("\nExcluding",length(excluded_missing),"variable(s) from MV model due to >",
                    scales::percent(max_missing_mv),"missingness:\n")
    if(verbose) cat(paste(" ",excluded_labels,paste0("(",round(missing_rates[excluded_missing] * 100, 1)
                                                     ,"% missing)"),collapse="\n"),"\n")
  }

  #Keep only variables below missingness threshold
  candidate_vars<-candidate_vars |> dplyr::filter(!variable %in% excluded_missing)
  if (nrow(candidate_vars)==0) {
    stop("No candidate variables remain after missingness filter. ",
         "Consider increasing max_missing_mv.")
  }

  #Compute n_events on complete cases only (Avoid overcounting dropped observations due to missingness)
  model_vars<-c(time_var,event_var,grouping_variable,candidate_vars$variable)
  complete_idx<-complete.cases(df[,intersect(model_vars,names(df))])
  n_complete<-sum(complete_idx)
  n_dropped<-nrow(df) - n_complete
  n_events<-sum(df[[event_var]][complete_idx]==1,na.rm=TRUE)

  # Warning if too many observations are dropped after missingness filter
  if (n_dropped > warn_missing_n) {

    warning("\nWarning:",n_dropped,"observation(s) dropped due to remaining missingness",
            "in the MV model analytic sample (n =",n_complete,",","events =",n_events,").\n")

    #Identify which variables still have missing values
    remaining_missing<-candidate_vars$variable[sapply(candidate_vars$variable,function( v ) {
      any(is.na(df[[v]]))
    })]
    if (length(remaining_missing) > 0) {
      remaining_labels<-meta_table$label[meta_table$variable %in% remaining_missing]
      if(verbose) cat("Variables with remaining missingness:\n")
      if(verbose) cat(paste(" ",remaining_labels,paste0("(",round(sapply(remaining_missing,function( v )
        mean(is.na(df[[ v ]]))) * 100, 1),"% missing)"),collapse="\n"),"\n")
    }
  }

  #Compute max_predictors
  max_predictors<-if (is.null(max_vars)) floor(n_events / min_epv) else max_vars
  if(verbose) cat("\nMV model: n =",n_complete,"| events =",n_events,"| max predictors =",max_predictors,"\n")

  #Define initial possible candidate variables
  mv_vars_init<-setdiff(candidate_vars$variable,grouping_variable)

  #Handle forced variables
  if (!is.null(force_mv_vars)) {
    #Validate forced vars exist in df and candidate_vars
    missing_forced<-setdiff(force_mv_vars, names(df))
    if (length(missing_forced) > 0) {
      warning("force_mv_vars: variables not found in df and will be skipped: ",
              paste(missing_forced,collapse=", "))
      force_mv_vars<-intersect(force_mv_vars,names(df))
    }

    not_in_candidates<-setdiff(force_mv_vars, candidate_vars$variable)
    if (length(not_in_candidates) > 0 && verbose) {
      cat("\nNote: forced variables not in candidate list","(will be added regardless):",
          paste(not_in_candidates,collapse=", "),"\n")
    }

    if (verbose) cat("\nForcing variables into model before selection:",
                     paste(force_mv_vars,collapse=", "),"\n")

    #Remove forced vars from candidate pool
    mv_vars_init<-setdiff(mv_vars_init,force_mv_vars)
    candidate_vars<-candidate_vars |> dplyr::filter(!variable %in% force_mv_vars)
  }

  #Run forward selection model
  main_run<-mv_variable_selection(data=df
                                  ,candidate_vars=candidate_vars
                                  ,mv_vars_pool=mv_vars_init
                                  ,time_var=time_var
                                  ,event_var=event_var
                                  ,grouping_variable=grouping_variable
                                  ,meta_table=meta_table
                                  ,min_epv=min_epv
                                  ,lrt_p_cutoff=lrt_p_cutoff
                                  ,max_predictors=max_predictors
                                  ,cluster_var_p_cutoff=cluster_var_p_cutoff
                                  ,block_vars=block_vars
                                  ,selection_strategy=selection_strategy
                                  ,cie_threshold=cie_threshold
                                  ,cie_final_cleanup=cie_final_cleanup
                                  ,forced_selected=force_mv_vars
                                  ,all_levels_sig=all_levels_sig
                                  ,verbose=verbose)

  selected<-main_run$selected        #Selected variable for MV model
  best_model<-main_run$final_model   #Final (best) MV model extrated from analysis

  #Generate step summary table
  summary_table <- tibble()
  if (length(selected) > 0) {

    # Start from forced vars model if any, otherwise base model
    if (!is.null(force_mv_vars) && length(force_mv_vars) > 0) {
      init_formula  <- as.formula(
        paste0("Surv(", time_var, ",", event_var, ") ~ ",
               paste(c(grouping_variable, force_mv_vars), collapse = " + "))
      )
      current_model <- tryCatch(
        coxph(init_formula, data = df),
        error = function(e) coxph(as.formula(
          paste0("Surv(", time_var, ",", event_var, ") ~ ", grouping_variable)
        ), data = df)
      )
    } else {
      current_model <- coxph(as.formula(
        paste0("Surv(", time_var, ",", event_var, ") ~ ", grouping_variable)
      ), data = df)
    }

    # Build label lookup once before the loop
    var_labels <- setNames(meta_table$label, meta_table$variable)

    step_counter <- 0
    for (i in seq_along(selected)) {
      v <- selected[i]

      #Skip forced variables
      if (!is.null(force_mv_vars) && v %in% force_mv_vars) {
        if (verbose) cat("  Skipping summary row for forced variable '", v, "'\n", sep = "")
        next
      }
      step_counter <- step_counter + 1
      v_label <- if (!is.na(var_labels[v])) var_labels[v] else v #fallback to variable name

      new_vars<-c(grouping_variable, selected[seq_len( i )])
      new_formula<-as.formula(
        paste0("Surv(", time_var,",",event_var,") ~ ",paste(new_vars,collapse=" + ")))

      new_model<-coxph(new_formula,data=df)
      lrt<-tryCatch(anova(current_model,new_model,test="LRT"), error=function(e) NULL)
      p_lrt<-if (!is.null(lrt)) lrt$"Pr(>|Chi|)"[2] else NA

      formula_no_cluster<-as.formula(
        paste0("Surv(",time_var,",",event_var,") ~ ",paste(selected[seq_len(i)],collapse=" + ")))

      fit_no_cluster<-tryCatch(coxph(formula_no_cluster,data=df),error=function(e) NULL)
      cluster_global_p<-if (!is.null(fit_no_cluster)) {
        lrt_cl<-tryCatch(anova(fit_no_cluster,new_model,test="LRT"),error=function(e) NULL)
        if (!is.null(lrt_cl)) lrt_cl$"Pr(>|Chi|)"[2] else NA
      } else NA

      summary_table <-bind_rows(
        summary_table,
        tibble(
          step= step_counter ,
          added_variable= v_label ,
          variables=paste(c(grouping_variable,var_labels[selected[seq_len( i )]]),collapse=", "),
          n_vars= i ,
          AIC=AIC(new_model),
          LRT_p=round(p_lrt,4),
          cluster_global_p=round(cluster_global_p,4),
          epv=round(n_events / ( i + 1), 1)))

      current_model<-new_model
    }
  }

  ## (optional) Bootstrap model stability ##
  boot_stability<-NULL
  if (run_bootstrap) {

    if(verbose) cat("\nRunning bootstrap stability analysis (n =",n_bootstrap,")...\n")
    set.seed(seed)

    #Helper function for bootstrapping
    boot_one<-function( i ) {
      boot_df<-df[sample(nrow(df),nrow(df),replace=TRUE), ]
      tryCatch(mv_variable_selection(data=boot_df
                                     ,candidate_vars=candidate_vars
                                     ,mv_vars_pool=mv_vars_init
                                     ,time_var=time_var
                                     ,event_var=event_var
                                     ,grouping_variable=grouping_variable
                                     ,meta_table=meta_table
                                     ,min_epv=min_epv
                                     ,max_predictors=max_predictors
                                     ,cluster_var_p_cutoff=cluster_var_p_cutoff
                                     ,lrt_p_cutoff=lrt_p_cutoff
                                     ,block_vars=block_vars
                                     ,selection_strategy=selection_strategy
                                     ,cie_threshold=cie_threshold
                                     ,cie_final_cleanup=cie_final_cleanup
                                     ,forced_selected=force_mv_vars
                                     ,all_levels_sig=all_levels_sig
                                     ,verbose=FALSE
      )$selected,error=function(e) character(0))}

    if (parallel) {
      cl<-makeCluster(n_cores)
      on.exit(stopCluster(cl),add=TRUE)
      clusterExport(cl,varlist=c("mv_variable_selection","df","meta_table",
                                    "candidate_vars","grouping_variable","time_var",
                                    "event_var","min_epv","max_predictors",
                                    "cluster_var_p_cutoff","mv_vars_init",
                                    "lrt_p_cutoff","block_vars","selection_strategy",
                                    "cie_threshold","cie_final_cleanup","force_mv_vars"),
                    envir=environment())
      registerDoParallel(cl)
      registerDoRNG(seed=seed)

      boot_results<-foreach( i = seq_len(n_bootstrap)
                             ,.packages=c("dplyr", "survival")) %dopar% { boot_one( i ) }
    } else {
      boot_results<-lapply(seq_len(n_bootstrap), function( i ) {
        if ((i == 1 || i %% 50 == 0) && verbose)
          cat("  Bootstrap iteration", i , "/",n_bootstrap,"\n")
        boot_one( i )})
    }

    #Summarise stability
    all_boot_vars<-unique(unlist(boot_results))
    boot_stability<-tibble(
      variable=all_boot_vars,
      inclusion_rate=sapply(all_boot_vars, function( v ) {
        mean(sapply(boot_results, function(x) v %in% x))}),
      in_main_model=all_boot_vars %in% selected) |>
      dplyr::arrange(desc(inclusion_rate)) |>
      dplyr::mutate(stability=case_when(
        inclusion_rate >= 0.75~"High",
        inclusion_rate >= bootstrap_cutoff~"Moderate",
        TRUE~"Low"))

    boot_stability<-boot_stability |>
      dplyr::left_join(meta_table |> dplyr::select(variable,label),by="variable") |>
      dplyr::select(label,inclusion_rate,in_main_model,stability) |>
      dplyr::rename(variable=label)

    if(verbose) {
      cat("\nBootstrap complete. Variable stability summary:\n")
      print(boot_stability)
    }
  }

  ## Outputs ##
  mv_tbl<-if (!is.null(best_model)) {
    #Detect binary and 2-level categorical variables from final model using selected MV vars for show_single_row
    single_row_mv_vars<-meta_table |>
      dplyr::filter(variable %in% selected, type %in% c("binary","categorical"),
             !role %in% c("exposure","event","time")) |>
      dplyr::pull(variable) |> setdiff(grouping_variable) |> #grouping_variable should show all levels
    #Variables with more than two levels should show all levels
      intersect(selected[sapply(selected, function( v ) {
          if (!v %in% names(df)) return(FALSE)
          levs<-unique(df[[ v ]][!is.na(df[[ v ]])])
          length(levs)==2})])

  #Generate table of best model
    tbl_regression(best_model
                   ,exponentiate=TRUE
                   ,show_single_row=all_of(single_row_mv_vars)
                   ,label= if (length(labels) > 0) labels else NULL
                   ,pvalue_fun= ~ style_pvalue(.x,digits=2)) |>
      bold_p() |> bold_labels()
  } else NULL

  #Selection report if bootstrapping
  mv_selection_report <- if (run_bootstrap && !is.null(boot_stability) && nrow(summary_table) > 0) {
    summary_table |> dplyr::left_join(
      boot_stability |>  dplyr::select(variable,inclusion_rate,stability)
      ,by=c("added_variable"="variable")) |>
      dplyr::mutate(inclusion_rate=paste0(round(inclusion_rate * 100), "%")) |>
      dplyr::select(step,added_variable,AIC,LRT_p,cluster_global_p,epv,
                    inclusion_rate,stability) |>
      dplyr::rename(
        'Step'=step,'Variable'=added_variable,'AIC'=AIC,'LRT p-value'=LRT_p
             ,'Group p-value'=cluster_global_p,'Events per variable'=epv
             ,'Bootstrap inclusion'=inclusion_rate,'Stability'=stability)
  } else NULL

  return(list(
    summary=summary_table,
    best_model=best_model,
    mv_table=mv_tbl,
    selected_vars=selected,
    boot_stability=boot_stability,
    mv_selection_report=mv_selection_report))
}

#' Function to define best MVA strategy
#' @export
#' @examples
#' ## Placeholder Example ##
run_mv_sensitivity<-function(df,meta_table,time_var,event_var,vars=NULL,
                             grouping_variable,group_var_label="Group",
                             strategies=c("strict","two_phase","cie"),
                             cie_threshold=0.10,cie_final_cleanup=FALSE,
                             priority_levels=c("core","optional","low"),
                             priority_override=NULL,block_vars=TRUE,
                             min_epv=10,max_vars=NULL,max_missing_mv=0.20,
                             warn_missing_n=50,cluster_var_p_cutoff=0.05,
                             lrt_p_cutoff=0.05,label_vars=NULL,
                             relevel_vars=NULL,dichotomous_vars=NULL,
                             force_mv_vars=NULL,mv_independent=TRUE,
                             run_bootstrap=FALSE,n_bootstrap=200,
                             bootstrap_cutoff=0.50,parallel=FALSE,seed=123,
                             n_cores=parallel::detectCores()-1,
                             verbose=TRUE,all_levels_sig=FALSE) {

  # Validate strategies
  valid_strategies <- c("strict", "two_phase", "cie")
  invalid <- setdiff(strategies, valid_strategies)
  if (length(invalid) > 0) {
    stop("Invalid strategy/strategies: ", paste(invalid, collapse = ", "),
         ". Valid options are: ", paste(valid_strategies, collapse = ", "))
  }

  # Use meta_table$variable as default if vars not provided
  if (is.null(vars)) vars <- meta_table$variable

  if (verbose) cat("Running MV sensitivity analysis across",
                   length(strategies), "strategies:",
                   paste(strategies, collapse = ", "), "\n")

  # Run pipeline for each strategy
  results <- lapply(setNames(strategies, strategies), function(s) {
    if (verbose) cat("\n", strrep("=", 50), "\n",
                     "Strategy: ", toupper(s), "\n",
                     strrep("=", 50), "\n", sep = "")

    tryCatch(
      run_tab_uv_mv_pipeline(
        df                   = df,
        meta_table           = meta_table,
        vars                 = vars,
        time_var             = time_var,
        event_var            = event_var,
        grouping_variable    = grouping_variable,
        group_var_label      = group_var_label,
        descr_tab            = FALSE,
        uv_mod               = TRUE,
        mv_mod               = TRUE,
        mv_independent       = mv_independent,
        selection_strategy   = s,
        cie_threshold        = cie_threshold,
        cie_final_cleanup    = if (s == "cie") cie_final_cleanup else FALSE,
        priority_levels      = priority_levels,
        priority_override    = priority_override,
        block_vars           = block_vars,
        min_epv              = min_epv,
        max_vars             = max_vars,
        max_missing_mv       = max_missing_mv,
        warn_missing_n       = warn_missing_n,
        cluster_var_p_cutoff = cluster_var_p_cutoff,
        lrt_p_cutoff         = lrt_p_cutoff,
        label_vars           = label_vars,
        relevel_vars         = relevel_vars,
        dichotomous_vars     = dichotomous_vars,
        force_mv_vars        = force_mv_vars,
        all_levels_sig       = all_levels_sig,
        run_bootstrap        = run_bootstrap,
        n_bootstrap          = n_bootstrap,
        bootstrap_cutoff     = bootstrap_cutoff,
        parallel             = parallel,
        seed                 = seed,
        n_cores              = n_cores,
        verbose              = verbose
      ),
      error = function(e) {
        cat("  ERROR in strategy '", s, "': ", conditionMessage(e), "\n", sep = "")
        warning("Strategy '", s, "' failed with error: ", conditionMessage(e))
        NULL
      }
    )
  })

  # Remove failed strategies
  failed <- names(results)[sapply(results, is.null)]
  if (length(failed) > 0) {
    warning("The following strategies failed and are excluded from comparison: ",
            paste(failed, collapse = ", "))
    results <- results[!names(results) %in% failed]
  }

  if (length(results) == 0) {
    stop("All strategies failed. Cannot build comparison.")
  }

  # Helper to safely extract summary column: handles both renamed and original names
  get_summary_col <- function(summary, new_name, old_name) {
    if (new_name %in% names(summary)) return(summary[[new_name]])
    if (old_name %in% names(summary)) return(summary[[old_name]])
    return(NULL)
  }

  # Helper to check if model has selected variables
  has_predictors <- function(r) {
    if (is.null(r$mv_summary) || nrow(r$mv_summary) == 0) return(FALSE)
    n_col <- get_summary_col(r$mv_summary, "N variables", "n_vars")
    !is.null(n_col) && length(n_col) > 0 && max(n_col) > 0
  }

  # Build comparison summary
  comparison <- tibble(
    Strategy             = names(results),

    'N predictors'       = sapply(results, function(r) {
      if (is.null(r$mv_summary) || nrow(r$mv_summary) == 0) return(0L)
      n_col <- get_summary_col(r$mv_summary, "N variables", "n_vars")
      if (is.null(n_col)) return(0L)
      as.integer(max(n_col))
    }),

    'Selected variables' = sapply(results, function(r) {
      if (is.null(r$mv_summary) || nrow(r$mv_summary) == 0) return("None")
      var_col  <- get_summary_col(r$mv_summary, "Variable", "added_variable")
      n_col    <- get_summary_col(r$mv_summary, "N variables", "n_vars")
      vars_col <- get_summary_col(r$mv_summary, "Model variables", "variables")
      if (is.null(var_col) || is.null(n_col)) return("None")
      # Get the variables string from the last step
      if (!is.null(vars_col)) {
        gsub("\\bgroup\\b", group_var_label, vars_col[which.max(n_col)])
      } else {
        paste(var_col, collapse = ", ") }
    }),

    !!paste(group_var_label,'HR')  := sapply(results, function(r) {
      if (is.null(r$mv_model_final)) return(NA_real_)
      if (!has_predictors(r))        return(NA_real_)
      broom::tidy(r$mv_model_final) |>
        dplyr::filter(startsWith(term, "group")) |>
        dplyr::summarise(hr = round(exp(mean(estimate, na.rm = TRUE)), 3)) |>
        dplyr::pull(hr)
    }),

    !!paste(group_var_label,'p-value') := sapply(results, function(r) {
      if (is.null(r$mv_model_final)) return(NA_real_)
      if (!has_predictors(r))        return(NA_real_)
      broom::tidy(r$mv_model_final) |>
        dplyr::filter(startsWith(term, "group")) |>
        dplyr::pull(p.value) |>
        min(na.rm = TRUE) |>
        round(4)
    }),

    !!paste(group_var_label,'significant')  := sapply(results, function(r) {
      if (is.null(r$mv_model_final)) return("#N/A")
      if (!has_predictors(r))        return("NULL model")
      p <- broom::tidy(r$mv_model_final) |>
        dplyr::filter(startsWith(term, "group")) |>
        dplyr::pull(p.value) |>
        min(na.rm = TRUE)
      ifelse(p < cluster_var_p_cutoff, "Yes", "No")
    }),

    'Stability summary'  = sapply(results, function(r) {
      if (!run_bootstrap)
        return("Bootstrap not run")
      if (is.null(r$mv_model_stability) || nrow(r$mv_model_stability) == 0)
        return("#N/A")
      r$mv_model_stability |>
        dplyr::filter(in_main_model == TRUE) |>
        dplyr::mutate(summary=paste0(variable," = ",stability," (",round(inclusion_rate*100),"%)")) |>
        dplyr::pull(summary) |>
        paste(collapse = "; ")
    }),

    'N observations '    = sapply(results, function(r) {
      if (is.null(r$mv_model_final)) return(NA_integer_)
      as.integer(r$mv_model_final$n)
    }),

    'N events'           = sapply(results, function(r) {
      if (is.null(r$mv_model_final)) return(NA_integer_)
      as.integer(r$mv_model_final$nevent)
    })
  )

  if (verbose) {
    cat("\n", strrep("=", 50), "\n",
        "Sensitivity Analysis Summary\n",
        strrep("=", 50), "\n", sep = "")
    print(comparison)

    # Highlight agreement/disagreement
    sig_results <- comparison[[paste(group_var_label, "significant")]]
    sig_results <- sig_results[!is.na(sig_results) &
                                 sig_results != "NULL model" &
                                 sig_results != "Bootstrap not run"]

    if (length(sig_results) == 0) {
      cat("\nNo valid models to compare.\n")
    } else if (length(unique(sig_results)) == 1) {
      cat("\nAll strategies agree: group is",
          ifelse(sig_results[1] == "Yes", "SIGNIFICANT", "NOT SIGNIFICANT"),
          "across all selection approaches.\n")
    } else {
      cat("\nStrategies disagree on group significance: results are sensitive",
          "to the selection approach. Review individual models carefully.\n")
    }
  }

  return(list(
    results    = results,
    comparison = comparison
  ))
}


