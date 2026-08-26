### Machine Learning ###

#' Function to run one outer fold in nested CV
#' @export
#' @examples
#' ## Placeholder Example ##
run_outer_fold<-function(split,alpha_grid,n_inner,features,
                         time_var,event_var,lambda=c("lambda.min","lambda.1se")) {

  lambda<-match.arg(lambda)

  tr  <-rsample::analysis(split)
  tst <-rsample::assessment(split)

  X_tr <-as.matrix(dplyr::select(tr,  dplyr::all_of(features)))
  y_tr <-survival::Surv(tr[[time_var]],tr[[event_var]])

  X_ts <-as.matrix(dplyr::select(tst, dplyr::all_of(features)))
  y_ts <-survival::Surv(tst[[time_var]],tst[[event_var]])

  #Inner CV: sweep alpha_grid sequentially within this worker
  inner<-lapply(alpha_grid, function(a) {
    cv<-glmnet::cv.glmnet(
      X_tr,
      y_tr,
      family      = "cox",
      alpha       = a,
      nfolds      = n_inner,
      standardize = TRUE,
      parallel    = FALSE,   #No nested parallelism
      cox.ties    = "efron"
    )
    list(alpha = a, lambda = cv[[lambda]], cvm = min(cv$cvm))
  })
  best <-inner[[which.min(sapply(inner, `[[`, "cvm"))]]

  fit  <-glmnet::glmnet(
    X_tr,
    y_tr,
    family      = "cox",
    alpha       = best$alpha,
    lambda      = best$lambda,
    standardize = TRUE,
    cox.ties    = "efron"
  )

  risk<-as.numeric(predict(fit,newx=X_ts,type="link"))
  fold_ci<-survival::concordance(y_ts ~ risk,reverse=TRUE)$concordance

  list(risk        = risk,
       fold_cindex = fold_ci,
       best_alpha  = best$alpha,
       best_lambda = best$lambda,
       test_idx    = rownames(tst)
       )
}

#' Function to assemble out-of-fold risk vector from fold results
#' @export
#' @examples
#' ## Placeholder Example ##
assemble_risk<-function(fold_results,df) {
  cv_risk<-rep(NA_real_,nrow(df))
  names(cv_risk)<-rownames(df)
  for (r in fold_results) cv_risk[r$test_idx]<-r$risk
  cv_risk
}

#' Function to perform nested CV on train data with optional parallelization
#' @export
#' @examples
#' ## Placeholder Example ##
run_nested_cv<-function(data,features,time_var,event_var,id_col=NULL,
                        lambda=c("lambda.min","lambda.1se"),
                        alpha_grid=seq(0.1,1.0,by=0.1),
                        n_outer=5,n_inner=5,seed=42,
                        verbose=TRUE,n_cores=1,parallelize=TRUE){

  lambda<-match.arg(lambda)

  #Safety checks
  stopifnot(is.data.frame(data) || is.matrix(data))
  stopifnot(all(features %in% colnames(data)))
  stopifnot(time_var %in% colnames(data))
  stopifnot(event_var %in% colnames(data))
  stopifnot(n_cores >= 1)
  if (parallelize && n_cores == 1) {
    warning("parallelize = TRUE but n_cores = 1. Running sequentially.")
    parallelize<-FALSE
  }

  #Reproducibility
  set.seed(seed)

  if (parallelize) {
    #Close any lingering connections before starting new cluster
    if (foreach::getDoParWorkers() > 1) {
      doParallel::stopImplicitCluster()
    }
    cl<-parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)
    doRNG::registerDoRNG(seed)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    if (verbose) cat("\nRunning in parallel with",n_cores,"cores\n")
  } else {
    foreach::registerDoSEQ()   #ensures %dorng% works sequentially too
    if (verbose) cat("\nRunning sequentially\n")
  }

  # Create dataframe
  if (!is.null(id_col)) {
    #Use meaningful ID as rownames
    stopifnot(id_col %in% colnames(data))
    stopifnot(is.character(data[[id_col]]))
    train_df           <-as.data.frame(data[,c(id_col,features,time_var,event_var)])
    rownames(train_df) <-data[[id_col]]
    train_df<-train_df |> dplyr::select(-dplyr::all_of(id_col))
  } else {
    #Fall back to sequential integers
    train_df           <-as.data.frame(data[,c(features,time_var,event_var)])
    rownames(train_df) <-as.character(seq_len(nrow(train_df)))
  }

  #Get nested folds
  nested_folds_full<-rsample::nested_cv(
    data    = train_df,
    outside = rsample::vfold_cv(v = n_outer, strata = event_var),
    inside  = rsample::vfold_cv(v = n_inner, strata = event_var)
  )

  #Run outer folds
  #Inner alpha sweep runs sequentially within each worker
  fold_res_full<-foreach::foreach(
    s         = nested_folds_full$splits,
    .packages = c("rsample","glmnet","survival","dplyr"),
    .export   = "run_outer_fold"
  ) %dorng% {
    run_outer_fold(
      split      = s,
      alpha_grid = alpha_grid,
      n_inner    = n_inner,
      features   = features,
      event_var  = event_var,
      time_var   = time_var,
      lambda     = lambda
    )
  }

  ##Compile results##
  fold_cindex_full <-sapply(fold_res_full,`[[`,"fold_cindex")
  best_alphas_full <-sapply(fold_res_full,`[[`,"best_alpha")
  most_freq_alpha  <-as.numeric(names(which.max(table(best_alphas_full))))
  y_train          <-survival::Surv(train_df[[time_var]],train_df[[event_var]])
  cv_risk_full     <-assemble_risk(fold_results=fold_res_full,df=train_df)
  global_cindex    <-survival::concordance(y_train ~ cv_risk_full,reverse=TRUE)$concordance

  if(verbose){
    cat("  Features used:      ", length(features), "\n")
    cat("  Mean fold C-index:  ", round(mean(fold_cindex_full), 3),
        "+/-", round(sd(fold_cindex_full), 3), "\n")
    cat("  Global C-index:     ", round(global_cindex, 3), "\n")
    cat("  Most frequent alpha:", most_freq_alpha, "\n")
    cat("  Alpha per fold:     ", best_alphas_full, "\n")
  }

  #Adjust dataset to export
  train_df_exp<-train_df |>
    tibble::rownames_to_column(var = if (!is.null(id_col)) id_col else "row_id") |>
    tibble::as_tibble()

  return(list(
    fold_res_full   = fold_res_full,
    cv_risk         = cv_risk_full,
    global_cindex   = global_cindex,
    fold_cindex     = fold_cindex_full,
    mean_cindex     = mean(fold_cindex_full),
    sd_cindex       = sd(fold_cindex_full),
    best_alphas     = best_alphas_full,
    most_freq_alpha = most_freq_alpha,
    y_train         = y_train,
    data_model      = train_df_exp,
    lambda          = lambda,
    proteins_model  = features
  ))
}

#' Function to perform feature stability selection on train data with optional parallelization
#' @export
#' @examples
#' ## Placeholder Example ##
run_stability_selection<-function(data,features,time_var,event_var,alpha,lambda,
                                  stab_cutoff=0.75,b=100,seed=42,n_cores=1,
                                  parallelize=TRUE,verbose=TRUE) {



  #Safety checks
  stopifnot(is.data.frame(data) || is.matrix(data))
  stopifnot(all(features %in% colnames(data)))
  stopifnot(time_var %in% colnames(data))
  stopifnot(event_var %in% colnames(data))
  stopifnot(alpha >= 0 && alpha <= 1)
  stopifnot(stab_cutoff > 0 && stab_cutoff < 1)
  stopifnot(b > 0)
  if (parallelize && n_cores == 1) {
    warning("parallelize = TRUE but n_cores = 1. Running sequentially.")
    parallelize <- FALSE
  }

  #Build X and y
  X <-as.matrix(data[,features])
  y <-survival::Surv(data[[time_var]],data[[event_var]])

  #Parallelization setup
  if (parallelize) {
    #Close any lingering connections before starting new cluster
    if (foreach::getDoParWorkers() > 1) {
      doParallel::stopImplicitCluster()
    }
    cl<-parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)
    doRNG::registerDoRNG(seed)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    if (verbose) cat("\nRunning stability selection in parallel with",n_cores,"cores\n")
  } else {
    foreach::registerDoSEQ()
    if (verbose) cat("\nRunning stability selection sequentially\n")
  }

  if (verbose) cat("\nAlpha:",alpha,"| Cutoff:",stab_cutoff,"| B:",b,"| Lambda:",lambda,"\n")

  #Manual bootstrap loop
  #Each iteration: subsample 50% without replacement, fit glmnet, record selected features
  stab_matrix <- foreach::foreach(
    b_i       = seq_len(b),
    .packages = c("glmnet", "survival"),
    .combine  = "cbind"
  ) %dorng% {

    #50% subsampling without replacement (Meinshausen & Buhlmann 2010)
    idx   <-sample(nrow(X),size=floor(nrow(X)/2),replace=FALSE)
    X_sub <-X[idx,,drop=FALSE]
    y_sub <-y[idx,]

    #Fit glmnet with cross-validated lambda
    cv<-glmnet::cv.glmnet(
      X_sub,
      y_sub,
      family      = "cox",
      alpha       = alpha,
      standardize = TRUE,
      parallel    = FALSE,
      cox.ties    = "efron"
    )

    fit<-glmnet::glmnet(
      X_sub,
      y_sub,
      family      = "cox",
      alpha       = alpha,
      lambda      = cv[[lambda]],
      standardize = TRUE,
      cox.ties    = "efron"
    )

    #Record which features were selected (non-zero coefficient)
    coef_mat<-as.matrix(coef(fit))
    as.numeric(features %in% rownames(coef_mat)[coef_mat[,1] != 0])
  }

  #Compute selection frequencies
  rownames(stab_matrix)  <-features
  stab_freq_all          <-rowSums(stab_matrix) / b
  names(stab_freq_all)   <-features

  #Apply cutoff to select stable features
  selected_features <-names(stab_freq_all[stab_freq_all >= stab_cutoff])
  stab_freq         <-stab_freq_all[selected_features]

  if (verbose) {
    cat("  Features selected:",length(selected_features),"/",length(features),"\n")
    cat("  Top 10 by frequency:\n")
    top10<-sort(stab_freq_all,decreasing=TRUE)[1:min(10,length(features))]
    print(round(top10,2))
  }

  return(list(
    stab_matrix       = stab_matrix,       # full B x p selection matrix
    stab_freq_all     = stab_freq_all,     # frequency for ALL features
    selected_features = selected_features, # features above cutoff
    stab_freq         = stab_freq,         # frequency for selected features only
    alpha             = alpha,
    stab_cutoff       = stab_cutoff,
    lambda            = lambda,
    b                 = b
  ))
}

#' Function to rank features and run panel reduction
#' @export
#' @examples
#' ## Placeholder Example ##
run_panel_reduction<-function(data,features,stab_freq,time_var,event_var,
                              alpha,coefs,lambda,min_size=5,tolerance=0.02,
                              n_outer=5,n_inner=5,seed=42,id_col=NULL,
                              verbose=TRUE,n_cores=1,parallelize=TRUE) {

  #Safety checks
  stopifnot(is.data.frame(data) || is.matrix(data))
  stopifnot(all(features %in% colnames(data)))
  stopifnot(time_var %in% colnames(data))
  stopifnot(event_var %in% colnames(data))
  stopifnot(all(names(coefs) %in% features))
  stopifnot(all(features %in% names(stab_freq)))
  stopifnot(min_size >= 1)
  stopifnot(tolerance > 0 && tolerance < 1)
  if (parallelize && n_cores == 1) {
    warning("parallelize = TRUE but n_cores = 1. Running sequentially.")
    parallelize<-FALSE
  }

  #Reset rownames
  data<-as.data.frame(data)
  if (!is.null(id_col)) {
    stopifnot(id_col %in% colnames(data))
    stopifnot(is.character(data[[id_col]]))
    rownames(data)<-data[[id_col]]
    data<-data |> dplyr::select(-dplyr::all_of(id_col))
  } else {
    rownames(data)<-as.character(seq_len(nrow(data)))
  }

  #Rank features by importance score
  importance      <-abs(coefs) * stab_freq[features]
  ranked          <-names(sort(importance,decreasing=TRUE))
  candidate_sizes <-seq(min_size,length(ranked))

  if (verbose) {
    cat("\nRunning panel reduction\n")
    cat("  Stable features:    ",length(features),"\n")
    cat("  Candidate sizes:    ",min_size,"to",length(ranked),"\n")
    cat("  Tolerance:          ",tolerance,"\n")
    cat("  Top 10 by importance:\n")
    print(round(sort(importance,decreasing=TRUE)[1:min(10,length(importance))],4))
  }

  #Parallelization setup
  if (parallelize) {
    cl<-parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)
    doRNG::registerDoRNG(seed)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    if (verbose) cat("\nRunning in parallel with", n_cores, "cores\n")
  } else {
    foreach::registerDoSEQ()
    if (verbose) cat("\nRunning sequentially\n")
  }

  #Perform Nested CV at each panel size
  #Panel sizes run in parallel, nested CV within each worker runs sequentially
  results<-foreach::foreach(
    i         = seq_along(candidate_sizes),
    .packages = c("rsample","glmnet","survival","dplyr"),
    .export   = c("run_outer_fold","assemble_risk","id_col")
    ) %dorng% {
    N       <-candidate_sizes[ i ]
    panel_N <-ranked[seq_len(N)]

    #Build data for this panel size
    df_N           <-as.data.frame(data[,c(panel_N,time_var,event_var)])
    rownames(df_N) <-rownames(data)

    #Get nested folds
    folds_i<-rsample::nested_cv(
      data    = df_N,
      outside = rsample::vfold_cv(v = n_outer, strata = event_var),
      inside  = rsample::vfold_cv(v = n_inner, strata = event_var)
    )

    #Run outer folds sequentially within worker
    fold_res_i<-lapply(
      folds_i$splits,
      run_outer_fold,
      alpha_grid = alpha,
      n_inner    = n_inner,
      features   = panel_N,
      time_var   = time_var,
      event_var  = event_var,
      lambda     = lambda
      )

    fc<-sapply(fold_res_i,`[[`,"fold_cindex")
    lambda_i<-median(sapply(fold_res_i,`[[`,"best_lambda"))

    cv_risk_i<-assemble_risk(fold_results=fold_res_i,df=df_N)

    y_i        <-survival::Surv(df_N[[time_var]],df_N[[event_var]])
    global_ci  <-survival::concordance(y_i ~ cv_risk_i,reverse=TRUE)$concordance

    list(size        = N,
         mean_cindex = mean(fc),
         sd_cindex   = sd(fc),
         cv_risk     = cv_risk_i,
         y_model     = y_i,
         lambda      = lambda_i,
         global_ci   = global_ci)

  }

  #Compile reduction curve
  res_df<-do.call(rbind,lapply(results,function(x) {
    data.frame(size        = x$size,
               mean_cindex = x$mean_cindex,
               sd_cindex   = x$sd_cindex,
               lambda      = x$lambda,
               global_ci   = x$global_ci)
  }))

  #Find optimal size: smallest panel within tolerance of best C-index
  C_best          <-max(res_df$mean_cindex)
  C_floor         <-C_best - tolerance
  eligible_idx    <-which(res_df$mean_cindex >= C_floor)    #only models that did not pass tolerance
  optimal_idx     <-eligible_idx[which.min(res_df$size[eligible_idx])] #models with minimal number of proteins
  y_model         <-results[[optimal_idx]]$y_model
  panel_size      <-res_df$size[optimal_idx]
  panel_cv_risk   <-results[[optimal_idx]]$cv_risk
  panel_cindex    <-res_df$mean_cindex[optimal_idx]
  panel_global_ci <-res_df$global_ci[optimal_idx]
  panel_sd        <-res_df$sd_cindex[optimal_idx]
  panel_lambda    <-res_df$lambda[optimal_idx]
  final_panel     <-ranked[seq_len(panel_size)]

  if (verbose) {
    cat("\nPanel reduction complete\n")
    cat("  Best C-index:      ",round(C_best,3),"\n")
    cat("  Optimal size:      ",panel_size,"\n")
    cat("  Panel C-index:     ",round(panel_cindex,3),"\n")
    cat("  Global C-index:    ",round(panel_global_ci, 3), "\n")
    cat("  SD:                ",round(panel_sd,3),"\n")
    cat("  Final panel:\n")
    print(final_panel)
  }

  #Plot decision graph
  panel_plot<-ggplot2::ggplot(res_df, aes(x = size, y = mean_cindex)) +
    geom_ribbon(aes(ymin = mean_cindex - sd_cindex,ymax = mean_cindex + sd_cindex),
                alpha = 0.15, fill = "steelblue") +
    geom_line(colour = "steelblue", linewidth = 0.8) +
    geom_point(colour = "steelblue", size = 2) +
    geom_hline(yintercept = C_floor, linetype = "dashed",
               colour = "red", linewidth = 0.7) +
    geom_hline(yintercept = C_best, linetype = "dotted",
               colour = "darkgreen", linewidth = 0.7) +
    geom_vline(xintercept = panel_size, linetype = "dashed",
               colour = "orange", linewidth = 0.7) +
    geom_point(data = res_df[res_df$size==panel_size,],
               aes(x = size, y = mean_cindex),
               colour = "orange", size = 4, shape = 18) +
    annotate("text", x = panel_size + 0.5, y = min(res_df$mean_cindex),
             label = paste0("Optimal\nn = ", panel_size),
             hjust = 0, colour = "orange", size = 3.5) +
    annotate("text", x = max(res_df$size),
             y = C_floor + 0.002,
             label = paste0("Tolerance floor (", round(C_floor, 3), ")"),
             hjust = 1, colour = "red", size = 3.5) +
    annotate("text", x = max(res_df$size),
             y = C_best + 0.002,
             label = paste0("Best C-index (", round(C_best, 3), ")"),
             hjust = 1, colour = "darkgreen", size = 3.5) +
    scale_x_continuous(breaks = seq(min(res_df$size),
                                    max(res_df$size), by = 2)) +
    labs(x     = "Number of proteins in panel",
         y     = "Nested CV C-index",
         title = "Panel reduction curve",
         subtitle = paste0("Optimal panel: ", panel_size, " proteins",
                           " (C-index = ", round(panel_cindex, 3), ")")) +
    theme_classic() +
    theme(plot.title    = element_text(face = "bold"),
          plot.subtitle = element_text(colour = "grey40"),
          axis.text     = element_text(size = 10),
          axis.title    = element_text(size = 11))

  #Adjust dataset to export
  data_exp<-data |>
    tibble::rownames_to_column(var = if (!is.null(id_col)) id_col else "row_id") |>
    tibble::as_tibble()

  return(list(
    res_df        = res_df,
    data_model    = data_exp,
    best_alpha    = alpha,
    best_lambda   = panel_lambda,
    rank_features = ranked,
    importance    = importance,
    y_panel       = y_model,
    panel_size    = panel_size,
    cv_risk       = panel_cv_risk,
    final_panel   = final_panel,
    best_cindex   = C_best,
    global_cindex = panel_global_ci,
    panel_cindex  = panel_cindex,
    panel_sd      = panel_sd,
    panel_plot    = panel_plot
  ))
}

#' Function to test which seed is best to separate train and test data
#' @export
#' @examples
#' ## Placeholder Example ##
get_best_seed<-function(train_data,test_data,event_var,max_seed=2026){
  seed_results<-do.call(rbind, lapply(1:max_seed, function(s) {
    set.seed(s)
    idx      <- caret::createDataPartition(train_data[[event_var]], p = 0.8, list = FALSE)
    train_er <- mean(train_data[[event_var]][idx])
    test_er  <- mean(test_data[[event_var]][-idx])
    diff     <- round(abs(train_er - test_er), 3)
    data.frame(seed     = s,
               train_er = round(train_er, 3),
               test_er  = round(test_er, 3),
               diff     = diff)
  }))

  #Sort by smallest difference
  seed_results<-seed_results[order(seed_results$diff),]

  return(seed_results)
}

#' Function to compute Univariate CoxPH
#' @export
#' @examples
#' ## Placeholder Example ##
run_coxph<-function(data,features,time_var,event_var,feat_labs,verbose=TRUE){

  #Build formula
  formula_str<-paste0("survival::Surv(",time_var,",",event_var,") ~ ",
                      paste(features,collapse=" + "))

  #Fit model
  cox_mod <-survival::coxph(as.formula(formula_str),data=data)
  sum_mod <-summary(cox_mod)

  #Extract metrics into tidy dataframe
  coef_mat <-sum_mod$coefficients
  ci_mat   <-exp(confint(cox_mod))

  results_df<-data.frame(
    feature  = rownames(coef_mat),
    hr       = exp(coef_mat[,"coef"]),
    ci_lower = ci_mat[,1],
    ci_upper = ci_mat[,2],
    pval     = coef_mat[,"Pr(>|z|)"],
    stringsAsFactors = FALSE
  )
  cindex<-sum_mod$concordance[1]

  if (verbose) {
    cat("\nCoxPH model\n")
    cat("  Features:",paste(features,collapse=", "), "\n")
    cat("  C-index:",round(cindex, 3),"\n")
    results_df |>
      dplyr::mutate(
        hr       = round(hr, 3),
        ci_lower = round(ci_lower, 3),
        ci_upper = round(ci_upper, 3),
        pval     = format.pval(pval, digits = 3)
      ) |> print()
  }

  return(list(
    model      = cox_mod,
    results_df = results_df,
    cindex     = cindex,
    summary    = sum_mod
  ))
}


