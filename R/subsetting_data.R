### Subsetting data according to continuous and categorical variables ###

#' Function to make variable groups
#' @export
#' @examples
#' ## Placeholder Example ##
get_var_groups<-function(df,var,var_labs=NULL
                         ,ref_patterns=c("WT","wt","no")
                         ,ref_suffix="neg",alt_suffix="pos"
                         ,missing_suffix="missing"
                         ,verbose=TRUE) {

  #Use var name as label if none provided
  if(is.null(var_labs)){ var_labs<-var }
  is_missing<-is.na(df[[var]])
  status_vec<-trimws(as.character(df[[var]]))
  prefix<-paste0(var_labs, "_")

  #Detect reference / negative group
  pattern<-paste0("\\b(",paste(ref_patterns,collapse="|"),")\\b")
  is_neg<-grepl(pattern,status_vec,ignore.case=TRUE)
  var_groups<-list()

  #all samples
  var_groups$all<-rep(TRUE,length(status_vec))

  #Reference group
  var_groups[[paste0(prefix,ref_suffix)]]<-is_neg

  #Track missing values
  var_groups[[paste0(prefix,missing_suffix)]]<- is_missing

  #Alternative samples
  idx<-(!is_neg) & !is_missing
  var_values<-unique(status_vec[idx])
  var_values<-var_values[!is.na(var_values)]

  if(length(var_values)==0){
    if(verbose) message("No alternative levels detected for variable: ",var)
    return(var_groups)
  }

  #Combined alternative group
  var_groups[[paste0(prefix,alt_suffix)]]<-idx

  #Individual levels
  for(val in var_values){
    safe_val<-gsub("\\s+","_",val)
    safe_val<-gsub("[^A-Za-z0-9_]","",safe_val)
    var_groups[[paste0(prefix,safe_val)]]<-(status_vec==val) & !is_missing
  }
  return(var_groups)
}

#' Function to make protein groups
#' @export
#' @examples
#' ## Placeholder Example ##
get_ptn_groups<-function(df_normal,protein,df_disease=NULL
                         ,df_apply=NULL,normal_pattern="CD34") {
  #bins=c(min_AML,min_CD34+,median_CD34+(=0),max_CD34+,max_AML)
  #Expected levels: always CD34 first and then AML

  if(is.null(df_disease)) df_disease<-df_normal
  if(is.null(df_apply)) df_apply<-df_disease

  if(!protein %in% names(df_apply))
    stop("Protein not found in df_apply")

  #Create variable vectors according to dataframes
  var_name<-names(df_normal)[1]                      #df_normal=df with normal samples
  disease_samples<-as.numeric(df_disease[[protein]])  #df_disease=df with all disease samples
  apply_values<-as.numeric(df_apply[[protein]])      #df_applyf=subset of df with all disease samples (eg, df after sample qc)

  #Safety check
  normal_idx<-grepl(normal_pattern,df_normal[[var_name]],ignore.case=TRUE)
  normal_samples<-as.numeric(df_normal[[protein]][normal_idx])
  if(length(normal_samples)==0)
    stop("No normal samples found using pattern ",normal_pattern)

  #Compute bins
  bins<-c(
    -Inf,
    round(min(normal_samples),4),
    round(median(normal_samples),4),
    round(max(normal_samples),4),
    Inf
  )

  #Report logical vectors
  list(
    all   = rep(TRUE,length(apply_values)),
    high  = apply_values >= bins[3],
    low   = apply_values <  bins[3],
    above = apply_values >  bins[4],
    upper = apply_values >= bins[3] & apply_values <= bins[4],
    lower = apply_values >= bins[2] & apply_values <  bins[3],
    below = apply_values <  bins[2]
  )
}

#' Function to create patient subsets
#' @export
#' @examples
#' ## Placeholder Example ##
get_subsets<-function(df,group_lists,print_summary=TRUE,exclude_NA=TRUE) {
  # group_lists: named list of lists of logical vectors (max 2 lists)
  # e.g., list(ptn=ptn_groups,mut=mut_groups)
  # e.g., ptn_groups<-list(gp1=high_ptnX=df$ptnX > 1),
  #       mut_groups<-list(mutY_pos=df$mutY=="yes",mutY_neg=df$mutY=="no"))

  if(length(group_lists) !=2)
    stop("group_lists must contain exactly 2 named lists")

  #Check all vectors in both group lists match nrow(df)
  n<-nrow(df)
  for(g in group_lists) {
    for(nm in names(g)) {
      if(length(g[[nm]]) != n)
        stop("Vector '",nm,"' has length ",length(g[[nm]])," but df has ", n," rows")
    }
  }

  g1<-group_lists[[1]]
  g2<-group_lists[[2]]
  subset_list<-list()
  summary_table<-data.frame(subset_name=character(),n_rows=integer(),stringsAsFactors=FALSE)
  for(n1 in names(g1)) {
    if(exclude_NA && grepl("missing|NA",n1,ignore.case=TRUE)) next
    v1<-g1[[n1]]
    for(n2 in names(g2)) {
      if(exclude_NA && grepl("missing|NA",n2,ignore.case=TRUE)) next
      v2 <- g2[[n2]]
      subset_name<-paste0(n1, "_", n2)
      mask<- v1 & v2
      mask[is.na(mask)]<-FALSE
      n_rows<-sum(mask)
      if(n_rows==0) next
      subset_list[[subset_name]]<-df[mask,,drop = FALSE]
      summary_table<-rbind(summary_table,data.frame(subset_name=subset_name,n_rows=n_rows))
    }
  }
  if(print_summary && nrow(summary_table) > 0) {
    cat("Subset sizes (total subsets =",nrow(summary_table),"):\n")
    summary_sorted<-summary_table[order(summary_table$n_rows),]
    for(i in seq_len(nrow(summary_sorted))) {
      cat(sprintf("  %-30s : %d rows\n",summary_sorted$subset_name[ i ],summary_sorted$n_rows[ i ]))
    }
  }
  return(list(subsets=subset_list,summary=summary_table))
}

#' Function wrapper to generate subsets with many proteins and variable groups
#' @export
#' @examples
#' ## Placeholder Example ##
create_all_subsets<-function(df_normal,df_disease,df_apply
                             ,protein_list=NULL,variable_list=NULL,var_labs=NULL
                             ,normal_pattern="CD34",ref_patterns=c("WT","wt","no")
                             ,exclude_NA=TRUE,print_summary=TRUE) {
  all_results<-list()
  #Ensure at least one type of grouping
  if(is.null(protein_list) && is.null(variable_list))
    stop("Need at least one of protein_list or variable_list for subsetting")

  #If both are present, do full cross-product
  protein_list_use<-if(is.null(protein_list)) "" else protein_list
  variable_list_use<-if(is.null(variable_list)) "" else variable_list

  for(ptn in protein_list_use) {
    for(var in variable_list_use) {

      group_lists<-list()

      #Protein groups
      if(ptn !="") {
        if(is.null(df_normal)) stop("df_normal is required for protein-level subsetting")
        ptn_groups<-get_ptn_groups(
          df_normal=df_normal,
          protein=ptn,
          df_disease=df_disease,
          df_apply=df_apply,
          normal_pattern=normal_pattern)

        group_lists$ptn<-ptn_groups
      }

      #Variable groups
      if(var != "") {
        var_groups<-get_var_groups(
          df=df_apply,
          var=var,
          ref_patterns=ref_patterns,
          var_labs=if(!is.null(names(var_labs))) var_labs[[var]] else var_labs)

        group_lists$var<-var_groups
      }

      #Create subsets
      if(length(group_lists)==1) {

        #Only 1 grouping dimension
        name1<-names(group_lists)[1]
        g1<-group_lists[[name1]]

        #Apply exclude_NA consistently with get_subsets
        if(exclude_NA) g1<-g1[!grepl("missing|NA",names(g1),ignore.case=TRUE)]

        subsets<-lapply(g1,function(mask) df_apply[mask,,drop=FALSE])
        subsets<-Filter(function(s) nrow(s) > 0, subsets)     #drop empty
        summary_table<-data.frame(subset_name=names(subsets)
                                  ,n_rows=sapply(subsets,nrow),stringsAsFactors=FALSE)

      } else if(length(group_lists)==2) {

        #Two grouping dimensions
        result<-get_subsets(df_apply,group_lists,exclude_NA=exclude_NA,print_summary=print_summary)
        subsets<-result$subsets
        summary_table<-result$summary
      }

      # Store results
      key<-paste(c(ptn,var)[c(ptn,var) != ""], collapse="_")
      all_results[[key]]<-list(subsets=subsets,summary=summary_table)
    }
  }
  return(all_results)
}
