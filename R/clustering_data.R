### Clustering functions ###

#' Function to generate clusters with progenyClust
#' @export
#' @examples
#' ## Placeholder Example ##
progeny_clusters<-function(df,proteins,n_clusters,cluster_prefix='C',
                           ncluster_pgcl=2:10,size_pgcl=10,iteration_pgcl=100,
                           repeats_pgcl=10,nrandom_pgcl=10,method_pgcl='gap',
                           plot_gap=TRUE,time_var='surv_time',
                           event_var='status',relapse_var='relapse'){

  #Initial safety checks
  missing_proteins<-setdiff(proteins,colnames(df))
  if (length(missing_proteins) > 0)
    stop("\nProteins not found in df: ",paste(missing_proteins,collapse=", "),"\n")

  missing_vars<-setdiff(c(time_var,event_var),colnames(df))
  if (length(missing_vars) > 0)
    stop("\nOutcome variables not found in df: ",paste(missing_vars,collapse=", "),"\n")

  missing_relapse<-setdiff(relapse_var,colnames(df))
  if (length(missing_relapse) > 0)
    warning("\nRelapse variable not found. Dataframe with complete cases for relapse will not be created\n")

  if (!all(sapply(df[,proteins,drop=FALSE],is.numeric)))
    stop("Some protein columns are not numeric")

  #Adjust tibble for clustering
  df<-as.data.frame(df)
  pg_mtx<-as.matrix(df[,proteins,drop=FALSE])

  #Cluster with progenyClust
  pgcl<-progenyClust(data=pg_mtx,FUNclust=hclust.progenyClust,
                     method=method_pgcl,ncluster=ncluster_pgcl,
                     size=size_pgcl,iteration=iteration_pgcl,
                     repeats=repeats_pgcl,nrandom=nrandom_pgcl)

  if (isTRUE(plot_gap)) {
    plot(pgcl)
    gap_plot<-recordPlot()
  } else {
    gap_plot<-NULL
  }

  #Select desired clusters
  cluster_cols<-paste0("C",n_clusters)

  #Clustering safety checks
  missing_cols<-setdiff(cluster_cols,colnames(pgcl$cluster))
  if (length(missing_cols) > 0)
    stop("Some cluster columns were not found: ",paste(missing_cols,collapse=", "))

  if (nrow(pgcl$cluster) != nrow(df))
    stop("Row mismatch between input df and progenyClust output. Check for NAs in protein columns")

  #Attach clusters to main df
  clusters<-as.data.frame(pgcl$cluster[,cluster_cols,drop=FALSE])
  colnames(clusters)<-paste0('cluster_',cluster_prefix,n_clusters)
  df<-dplyr::bind_cols(df,clusters)
  df<-tibble::as_tibble(df)

  if (length(cluster_cols)==1){
    names(df)[grepl("^cluster_",names(df))]<-"cluster"
    df<-prog_relevel(df=df,group_var='cluster',time_var=time_var,event_var=event_var)
    levels(df$cluster)<-paste0(cluster_prefix,seq_along(levels(df$cluster)))
  }
  if (length(cluster_cols) > 1){
    cluster_names<-names(df)[grepl("^cluster_",names(df))]
    for( i in seq_along(cluster_names)){
      df<-prog_relevel(df=df,group_var=cluster_names[ i ],time_var=time_var,event_var=event_var)
      levels(df[[cluster_names[ i ]]])<-paste0(cluster_prefix,seq_along(levels(df[[cluster_names[ i ]]])))
    }
  }

  #Create df for remission/relapse analysis
  if (relapse_var %in% names(df)) {
    df_rem<-df[!is.na(df[[relapse_var]]),]
  } else {
    df_rem<-NULL
  }

  return(list(df_full=df,df_rem=df_rem,gap_plot=gap_plot))
}

#' Function to relevel grouping variable according to prognnosis
#' @export
#' @examples
#' ## Placeholder Example ##
prog_relevel<-function(df,group_var,time_var='surv_time',event_var='status') {

  #Convert to character to handle both numeric and character group vars
  df[[group_var]]<-as.character(df[[group_var]])
  levels_vec<-sort(unique(df[[group_var]]))

  #Determine the median outcome (e.g. survival) of each factor level
  median_prog<-sapply(levels_vec, function(lvl) {
    group_level<-df[[group_var]]==lvl
    fit<-survfit(Surv(df[[time_var]][group_level],df[[event_var]][group_level])~1)
    tbl_fit<-summary(fit)$table
    med<-unname(tbl_fit["median"])
    if (is.na(med)) {
      warning(paste("Median survival NA for level",lvl,": using rmean as fallback"))
      unname(tbl_fit["rmean"])
    } else med
  })

  ordered_levels<-names(sort(median_prog,decreasing=TRUE))
  df[[group_var]]<-factor(df[[group_var]],levels=ordered_levels)
  return(df)
}


