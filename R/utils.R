### General utilities ###

#' Function to adjust protein notation
#' @export
#' @examples
#' ## Placeholder Example ##
clean_protein_names<-function(protein_names,strict=TRUE){

  #Validate input
  if (is.null(protein_names)) { stop("protein_names cannot be NULL") }
  if (is.vector(protein_names) && !is.list(protein_names)) { protein_names<-list(group1=as.character(protein_names)) }
  if (!is.list(protein_names)) {
    stop("protein_names must be either a character vector or a list of character vectors")
  }

  #Convert list to df format
  group_names<-names(protein_names)
  if (is.null(group_names)) { group_names<-seq_along(protein_names) }
  protein_names_df<-do.call(rbind,lapply(seq_along(protein_names), function( i ){
    data.frame(group_id=group_names[ i ],original=protein_names[[ i ]],stringsAsFactors=FALSE)}))

  mapping_list<-list()
  for(i in seq_len(nrow(protein_names_df))){

    group_id<-protein_names_df$group_id[ i ] #ID of protein list (e.g. up vs down, RPPA vs MS, etc) mapped numerically
    original<-protein_names_df$original[ i ]   #Original dataset name
    base_name<-NA    #Original names without PTMs
    ptm<-NA          #Store PTMs
    export_name<-NA  #retrieves the name that should be used for analysis (STRINGDB,ENRICHR,etc.)
    reason<-NA       #map how hard was the change for each case

    #Remove Histone marks and set as histone_marks
    # Detect histone marks
    if (grepl("^H[1-4][A-Z]?K[0-9]+", original)) {
      reason<-"histone_mark"
      base_name<-sub("^(H[1-4](?:[AB])?).*","\\1",original,perl=TRUE) #Extract histone core (H3,H4, etc.)
      ptm<-sub("^[^K]+", "", original) #Extract PTM
      mapping_list[[i]] <- data.frame(
        group_id=group_id,
        original=original,
        base_name=base_name,
        ptm=ptm,
        export_name=NA,   # keep NA so it won't go to STRING/HGNC
        mapping_level=reason,
        stringsAsFactors=FALSE
      )
      next
    }

    #Remove and store PTM suffixes
    base_name<-sub("\\..*", "",original)
    ptm<-ifelse(grepl("\\.",original),sub("^[^.]+\\.","",original),NA)

    #Expand complex names
    if (grepl("_", base_name)){
      name_parts<-unlist(strsplit(base_name,"_"))
      root<-name_parts[1]
      suffixes<-name_parts[-1]
      expanded_names<-paste0(root,suffixes)
      reason<-"complex"
    } else {
      expanded_names<-base_name
      reason<-"simple"
    }

    #Create one row per analysed proteins in the final df
    mapping_list[[ i ]]<-data.frame(
      group_id=rep(group_id, length(expanded_names)),
      original=original,
      base_name=base_name,
      ptm=ptm,
      export_name=expanded_names,
      mapping_level=reason,
      stringsAsFactors=FALSE
    )
  }

  #Create mapping dataframe
  mapping_df<-dplyr::bind_rows(mapping_list)

  #Create PTM type variable
  mapping_df$ptm_type<-factor(
    dplyr::case_when(
      is.na(mapping_df$ptm)~"total",
      grepl("[STY][0-9]+", mapping_df$ptm)~'phospho',
      grepl("^[KR][0-9]+me[123]?$",mapping_df$ptm,ignore.case=TRUE)~"methyl",
      grepl("ub$", mapping_df$ptm,ignore.case = TRUE)~"ubiquitin",
      grepl("ac$", mapping_df$ptm,ignore.case=TRUE)~"acetyl",
      grepl("cle$", mapping_df$ptm,ignore.case=TRUE)~"cleaved",
      TRUE~"other"),
    levels=c("total","phospho","methyl","ubiquitin","acetyl","cleaved","other"))

  #Validate names with HGNChelper
  valid_names<-unique(na.omit(mapping_df$export_name))
  name_check<-HGNChelper::checkGeneSymbols(valid_names)
  colnames(name_check)<-c('export_name','HGNC_check','suggestion')
  mapping_df<-dplyr::left_join(mapping_df,name_check,by='export_name')
  mapping_df<-mapping_df[,c(colnames(mapping_df)[-1],colnames(mapping_df)[1])]
  mapping_df$final_names<-ifelse(is.na(mapping_df$HGNC_check),NA
                                 ,ifelse(mapping_df$HGNC_check==FALSE
                                         ,mapping_df$suggestion
                                         ,mapping_df$export_name))
  #Detect invalid names
  if (any(is.na(mapping_df$final_names) & mapping_df$mapping_level !='histone_mark')) {
    message1<-"\nSome protein names are invalid and have no HGNC suggestion. Double-check mapping table.\n"
    if (strict) stop(message1,call.=FALSE)
    else warning(message1,call.=FALSE)
  }

  if (any(mapping_df$mapping_level=='histone_mark')) {
    warning("\nHistone marks were delivered as input. Not valid as HGNC input.\n")
  }

  if (any(mapping_df$HGNC_check==FALSE & !is.na(mapping_df$suggestion),na.rm=TRUE)){
    message2<-"\nSome protein names were auto-corrected using HGNC suggestions. Double-check mapping table.\n"
    if (strict) stop(message2,call.=FALSE)
    else warning(message2,call.=FALSE)
  }

  return(list(name_mapping=mapping_df,
              original_names=unique(mapping_df$original),
              clean_names=unique(na.omit(mapping_df$final_names))
  ))
}

#' Function to get colors from pre-defined sets or used-defined sets
#' @export
#' @examples
#' ## Placeholder Example ##
get_mycolors<-function(n_colors=NULL,mode=c("list","vector")
                       ,color_list=NULL,color_vector=NULL) {

  mode<-match.arg(mode)

  #Safety checks
  if (!is.null(color_list) && !is.list(color_list)) {
    stop("color_list must be a list. If not passing a list, use color_vector instead.")
  }
  if (!is.null(color_vector) && is.list(color_vector)) {
    stop("color_vector cannot be a list. Use color_list instead.")
  }

  #Default colors
  default_color_vector<-c(
    'red2','blue3','green4','purple2','darkorange','darkcyan',
    'darkgoldenrod3','yellow2','deeppink2','lightsalmon2','darkturquoise'
  )

  default_color_list<-list(
    '2' =c('blue2','red2'),
    '3' =c('blue2','springgreen3','red2'),
    '4' =c('blue2','springgreen3','darkorange','red2'),
    '5' =c('purple2','blue2','springgreen3','darkorange','red2'),
    '6' =c('purple2','blue2','springgreen3','yellow3','darkorange','red2'),
    '7' =c('purple2','blue2','springgreen3','darkgreen','yellow3','darkorange','red2'),
    '8' =c('purple2','blue2','springgreen3','darkgreen','yellow3','darkgoldenrod3','darkorange','red2'),
    '9' =c('purple2','blue2','deepskyblue','springgreen3','darkgreen','yellow3','darkgoldenrod3','darkorange','red2'),
    '10'=c('purple2','blue2','deepskyblue','springgreen3','darkgreen','khaki3','yellow3','darkgoldenrod3','darkorange','red2')
  )

  #Apply defaults if user does not provide custom colors
  if (is.null(color_vector)) color_vector<-default_color_vector
  if (is.null(color_list)) color_list<-default_color_list

  if (mode=="vector") {
    if (is.null(n_colors)) return(color_vector)
    if (n_colors > length(color_vector)) {
      color_vector<-grDevices::colorRampPalette(color_vector)(n_colors)
    }
    return(color_vector[seq_len(n_colors)])
  }

  if (mode=="list"){
  if (is.null(n_colors)) return(color_list)
    color_key<-as.character(n_colors)
    color_palette<-color_list[[color_key]]
    if (is.null(color_palette)) {
      largest_palette<-color_list[[as.character(max(as.integer(names(color_list))))]]
      color_palette<-grDevices::colorRampPalette(largest_palette)(n_colors)
    }
    return(color_palette)
  }
}

#' Function to store labels
#' @export
#' @examples
#' ## Placeholder Example ##
set_var_labs<-function(data,meta_table) {
  labelled::var_label(data)<-meta_table$label[match(names(data),meta_table$variable)]
  return(data)
}











