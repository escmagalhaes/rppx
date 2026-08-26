#' Function to wrap descriptive tables, UVA and MVA in a single pipeline
#' @export
#' @examples
#' ## Placeholder Example ##
run_tab_uv_mv_pipeline<-function(df,meta_table,vars,time_var,event_var,
                                 grouping_variable,group_var_label="Group",
                                 descr_tab=TRUE,uv_mod=TRUE,mv_mod=TRUE,
                                 min_n=5,max_levels=8,min_freq=0.03,
                                 skip_pval=-9999,exclude_sig=NULL,
                                 min_freq_var=0.03,min_epv=10,uv_p_cutoff=0.20,
                                 priority_levels=c("core","optional","low"),
                                 priority_override=NULL,block_vars=TRUE,
                                 selection_strategy=c("strict","two_phase","cie"),
                                 cie_threshold=0.10,cie_final_cleanup=FALSE,
                                 max_vars=NULL,cluster_var_p_cutoff=0.05,
                                 lrt_p_cutoff=0.05,label_vars=NULL,all_levels_sig=FALSE,
                                 relevel_vars= NULL,dichotomous_vars=NULL,
                                 max_missing_mv=0.20,warn_missing_n=50,
                                 mv_independent=TRUE,force_mv_vars=NULL,
                                 run_bootstrap=FALSE,n_bootstrap=200,
                                 bootstrap_cutoff=0.50,parallel=FALSE,seed=123,
                                 n_cores=parallel::detectCores()-1,verbose=TRUE){

  #Match MV strategy function argument
  selection_strategy<-match.arg(selection_strategy)

  #Safety check
  if (!grouping_variable %in% names(df)) {
    stop(paste("grouping_variable '", grouping_variable, "' not found in df."))
  }

  #Apply releveling for specified variables
  if (!is.null(relevel_vars)) {
    for ( v in names(relevel_vars)) {
      if ( v %in% names(df)) {
        df[[ v ]]<-relevel(factor(df[[ v ]]),ref=relevel_vars[[ v ]])
      } else {
        warning(paste("relevel_vars: variable '", v ,"' not found in df. Skipping."))
      }
    }
  }

  #Apply labels
  names(df)[names(df)==grouping_variable]<-"group"
  df<-set_var_labs(df,meta_table)
  if (!is.null(group_var_label)) { labelled::var_label(df$group)<-group_var_label }

  ### Descriptive table ###
  if(isTRUE(descr_tab)){

    #Convert categorical/binary variables to factor
    cat_vars_tab<-meta_table |>
      dplyr::filter(type %in% c("binary","categorical"),include_table==TRUE
             ,!role %in% c("exposure","event","time")) |> dplyr::pull(variable)
    for( var in intersect(cat_vars_tab,names(df))){ df[[var]]<-as.factor(df[[var]]) }

    #All table variables
    all_tab_vars<-meta_table |>
      dplyr::filter(include_table==TRUE,!role %in% c("exposure","event","time")) |>
      dplyr::pull(variable)

    #Detect variables that have too low frequency or are too sparse to test
    flag_tab_vars<-rare_sparse_vars(df=df
                                    ,group_var="group"
                                    ,vars=cat_vars_tab
                                    ,mode="table"
                                    ,min_n=min_n
                                    ,max_levels=max_levels
                                    ,min_freq=min_freq)

    #Select all table variables
    all_tab_vars<-meta_table |>
      dplyr::filter(include_table==TRUE,!role %in% c("exposure","event","time")) |>
      dplyr::pull(variable)

    #Variables to test: all table vars minus grouping_variable and flagged vars
    tab_vars_to_test<-setdiff(all_tab_vars,unlist(flag_tab_vars))

    #Report skipped variables by label
    sparse_tab_vars_labs<-meta_table[meta_table$variable %in% flag_tab_vars$sparse_vars,'label',drop=TRUE]
    rare_tab_vars_labs<-meta_table[meta_table$variable %in% flag_tab_vars$rare_vars,'label',drop=TRUE]

    if (length(unlist(flag_tab_vars)) > 0 && verbose){
      cat("\nSkipping statistical tests for these variables:\n")
      cat("\nSparse: ",paste(sparse_tab_vars_labs,collapse=", "),"\n")
      cat("\nRare: ",paste(rare_tab_vars_labs,collapse=", "),"\n")
    }

    #Define statistical test for each variable
    stat_tab_tests<-define_stat_test(df=df
                                     ,grouping_variable="group"
                                     ,vars=tab_vars_to_test)

    #Detect binary yes/no variables in all categorical variables
    binary_tab_vars<-cat_vars_tab[sapply(cat_vars_tab, function( v ) {
      levs<-unique(df[[ v ]][!is.na(df[[ v ]])])
      length(levs)==2 && all(tolower(levs) %in% c("yes","no"))
    })]

    #Generate type and value lists for gtsummary
    type_tab_list<-setNames(lapply(binary_tab_vars,function(v) "dichotomous"),binary_tab_vars)
    value_tab_list<-setNames(lapply(binary_tab_vars,function(v) "yes"),binary_tab_vars)

    #Add custom dichotomous variables (non-yes/no binary variables)
    if (!is.null(dichotomous_vars)) {
      for (v in names(dichotomous_vars)) {
        if (v %in% names(df)) {
          type_tab_list[[ v ]]<-"dichotomous"
          value_tab_list[[ v ]]<-dichotomous_vars[[ v ]]
        } else {
          warning(paste("dichotomous_vars: variable '", v , "' not found in dataframe. Skipping."))
        }
      }
    }

    #Descriptive tables
    set.seed(seed) #For fisher's exact test with simulated p-values
    descriptive_tables<-suppressWarnings(get_descr_tab(df=df
                                      ,grouping_variable="group"
                                      ,vars=all_tab_vars
                                      ,vars_to_test=tab_vars_to_test
                                      ,tests=stat_tab_tests$tests
                                      ,args=stat_tab_tests$test_args
                                      ,skipped_vars=unlist(flag_tab_vars)
                                      ,exclude_sig=exclude_sig
                                      ,sig_type=type_tab_list
                                      ,sig_value=value_tab_list
                                      ,sig_labels=label_vars
                                      ,skip_pval=skip_pval))

    descriptive_table_full<-descriptive_tables$extended_table
    descriptive_table_significant<-descriptive_tables$significant_table
  } else {
    descriptive_table_full<-NULL
    descriptive_table_significant<-NULL
  }

  ### UV Coxph table ###
  if(isTRUE(uv_mod)){

    #Convert categorical/binary variables to factor
    cat_vars_uv<-meta_table |>
      dplyr::filter(type %in% c("binary","categorical"),include_uv_model==TRUE
             ,!role %in% c("event","time")) |> dplyr::pull(variable)
    for( var in intersect(cat_vars_uv,names(df))){ df[[var]]<-as.factor(df[[var]]) }

    #All UV variables
    all_uv_vars<-meta_table |>
      dplyr::filter(include_uv_model==TRUE,!role %in% c("event","time")) |>
      dplyr::pull(variable)

    #Detect variables that have too low frequency or are too sparse to test
    flag_uv_vars<-rare_sparse_vars(df=df
                                   ,group_var=event_var
                                   ,vars=cat_vars_uv
                                   ,mode="outcome"
                                   ,min_freq=min_freq_var
                                   ,min_epv=min_epv)

    #Select all table variables
    all_uv_vars<-meta_table |>
      dplyr::filter(include_uv_model==TRUE,!role %in% c("exposure","event","time")) |>
      dplyr::pull(variable)

    #Variables to test: all UV vars minus flagged vars
    uv_vars_to_test<-setdiff(all_uv_vars,unlist(flag_uv_vars))

    #Report skipped variables by label
    sparse_uv_vars_labs<-meta_table[meta_table$variable %in% flag_uv_vars$sparse_vars,'label',drop=TRUE]
    rare_uv_vars_labs<-meta_table[meta_table$variable %in% flag_uv_vars$rare_vars,'label',drop=TRUE]

    if (length(unlist(flag_uv_vars)) > 0 && verbose){
      cat("\nSkipping statistical tests for these variables:\n")
      cat("\nSparse: ",paste(sparse_uv_vars_labs,collapse=", "),"\n")
      cat("\nRare: ",paste(rare_uv_vars_labs,collapse=", "),"\n")
    }

    #Run UV CoxPH
    uv_model<-run_uv_coxph(df=df
                           ,vars=uv_vars_to_test
                           ,time_var=time_var
                           ,event_var=event_var
                           ,uv_p_cutoff=uv_p_cutoff
                           ,labels=label_vars
                           ,meta_table=meta_table
                           ,group_var_label=group_var_label)

    uv_table_full<-uv_model$uv_table
    uv_table_significant<-uv_model$sig_uv_table
    selected_mv_vars<-uv_model$selected_mv_vars
  } else {
    uv_table_full<-NULL
    uv_table_significant<-NULL
    selected_mv_vars<-NULL
  }

  ### MV forward selection model ###
  if(isTRUE(mv_mod)){

    if(!isTRUE(uv_mod)){
      stop("MV model variable selection depends on UV models. Please, run UV models first.\n")
    }

    # Allow MV without UV only if mv_independent = TRUE
    if (!isTRUE(uv_mod) && !isTRUE(mv_independent)) {
      stop("MV model variable selection depends on UV models. ",
           "Set mv_independent=TRUE to run MV without UV.\n")
    }

    # Determine MV candidates
    if (isTRUE(mv_independent)) {

      #Use all meta_table candidates regardless of UV results
      if (verbose) cat("\nMV independent mode: using all eligible meta_table candidates.\n")

      mv_candidates<-meta_table |>
        dplyr::filter(include_mv_model==TRUE
               ,!role %in% c("exposure","event","time")
               ,!is.na(priority)
               ,priority %in% priority_levels) |> dplyr::pull(variable)

      } else if (length(selected_mv_vars)==0) {

        #Automatic fallback with warning
        warning("No UV-significant variables found at p < ", uv_p_cutoff,
                ". Falling back to all eligible meta_table MV candidates. ",
                "Consider using mv_independent=TRUE explicitly.")

        mv_candidates<-meta_table |>
          dplyr::filter(include_mv_model==TRUE
                 ,!role %in% c("exposure","event","time")
                 ,!is.na(priority),
                 priority %in% priority_levels) |> dplyr::pull(variable)

      } else {
      mv_candidates<-selected_mv_vars
    }

    #Convert categorical variables to factor
    cat_vars_mv<-meta_table |>
      dplyr::filter(type %in% c("binary","categorical")
                    ,include_mv_model==TRUE
                    ,!role %in% c("exposure","event","time")) |> dplyr::pull(variable)
    for ( v in intersect(cat_vars_mv,names(df))) df[[ v ]]<-as.factor(df[[ v ]])

    #Run MV CoxPH
    mv_model<-run_mv_coxph(df=df,
                           meta_table=meta_table,
                           vars=mv_candidates,
                           time_var=time_var,
                           event_var=event_var,
                           grouping_variable="group",
                           priority_levels=priority_levels,
                           min_epv=min_epv,
                           max_vars=max_vars,
                           cluster_var_p_cutoff=cluster_var_p_cutoff,
                           labels=label_vars,
                           max_missing_mv=max_missing_mv,
                           warn_missing_n=warn_missing_n,
                           lrt_p_cutoff=lrt_p_cutoff,
                           priority_override=priority_override,
                           block_vars=block_vars,
                           run_bootstrap=run_bootstrap,
                           n_bootstrap=n_bootstrap,
                           bootstrap_cutoff=bootstrap_cutoff,
                           parallel=parallel,
                           seed=seed, #for bootstrapping
                           n_cores=n_cores,
                           selection_strategy=selection_strategy,
                           cie_threshold=cie_threshold,
                           cie_final_cleanup=cie_final_cleanup,
                           force_mv_vars=force_mv_vars,
                           all_levels_sig=all_levels_sig,
                           verbose=verbose)

  mv_summary <-mv_model$summary
  final_model<-mv_model$best_model
  mv_table<-mv_model$mv_table
  mv_model_stability<-mv_model$boot_stability
  mv_selection_report<-mv_model$mv_selection_report
  } else {
    mv_summary <-NULL
    final_model<-NULL
    mv_table<-NULL
    mv_model_stability<-NULL
    mv_selection_report<-NULL
  }

  return(list(
    processed_df=df,
    meta_table=meta_table,
    descr_table=descriptive_table_full,
    descr_table_sig=descriptive_table_significant,
    uv_table_full=uv_table_full,
    uv_table_sig=uv_table_significant,
    mv_summary=mv_summary,
    mv_model_final=final_model,
    mv_table=mv_table,
    mv_model_stability=mv_model_stability,
    mv_selection_report=mv_selection_report
  ))
}
