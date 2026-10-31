### General utilities functions ###

#' Clean and standardize protein names and notations
#'
#' Most useful for names with a post-translational modification (PTM) as suffix.
#' Processes protein names by removing and storing PTM suffixes,
#' expanding complex protein notations, detecting histone marks,
#' and validating protein names against HGNC symbols.
#'
#' @param protein_names A character vector or named list of character vectors
#'   containing raw protein names as they appear in input data.
#' @param strict Logical. If `TRUE` (default), stops on invalid or
#'   auto-corrected HGNC names. If `FALSE`, issues a warning instead.
#'
#' @return A named list with three elements:
#' \describe{
#'   \item{name_mapping}{A data frame with full mapping details. Key columns:
#'     \itemize{
#'       \item `original` — raw input name
#'       \item `base_name` — name after PTM removal
#'       \item `ptm` — extracted PTM suffix if present
#'       \item `ptm_type` — factor: `total`, `phospho`, `methyl`,
#'         `ubiquitin`, `acetyl`, `cleaved`, or `other`
#'       \item `export_name` — name prepared for external tools
#'         (STRINGdb, EnrichR)
#'       \item `HGNC_check` — logical: `TRUE` = valid HGNC symbol,
#'         `FALSE` = invalid but suggestion might be available,
#'         `NA` = not found (usually histone marks)
#'       \item `suggestion` — HGNC-suggested correction when
#'         `HGNC_check == FALSE`
#'       \item `final_names` — **always use this column downstream**.
#'         Consolidates all cases: valid names, auto-corrected names.
#'         `NA` for histone marks and unresolvable names.
#'     }
#'   }
#'   \item{original_names}{Character vector of unique original input names.}
#'   \item{clean_names}{Character vector of unique valid final names after
#'     HGNC correction. Ready for use with STRINGdb and EnrichR.}
#' }
#'
#' @details
#' The function handles four naming cases:
#' - **Simple names**: validated directly against HGNC
#' - **PTM suffixes**: split at the first `.` into base name and PTM
#' - **Complex notations**: `BASENAME_1_2` expanded to `BASENAME1`, `BASENAME2`
#' - **Histone marks**: detected by pattern `H[1-4]K[0-9]+`, excluded
#'   from HGNC validation. `final_names` will be `NA` for these.
#'
#' The primary use cases are: (1) stripping PTMs correctly before
#' submitting to pathway analysis or STRINGdb, and (2) correcting
#' non-standard protein names in raw input data.
#'
#' @examples
#' clean_protein_names(c("EGFR", "AKT_1_2_3.pS473", "GSK3_A_B", "NOTCH1.cle1744"))
#' clean_protein_names(c("MAPK1", "H3K27Me3"), strict = FALSE)
#'
#' @seealso [HGNChelper::checkGeneSymbols()]
#'
#' @export
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
    #Detect histone marks
    if (grepl("^H[1-4][A-Z]?K[0-9]+", original)) {
      reason<-"histone_mark"
      base_name<-sub("^(H[1-4](?:[AB])?).*","\\1",original,perl=TRUE) #Extract histone core (H3,H4, etc.)
      ptm<-sub("^[^K]+", "", original) #Extract PTM
      mapping_list[[ i ]]<-data.frame(
        group_id         = group_id,
        original         = original,
        base_name        = base_name,
        ptm              = ptm,
        export_name      = NA,   # keep NA so it won't go to STRING/HGNC
        mapping_level    = reason,
        stringsAsFactors = FALSE
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
      group_id         = rep(group_id,length(expanded_names)),
      original         = original,
      base_name        = base_name,
      ptm              = ptm,
      export_name      = expanded_names,
      mapping_level    = reason,
      stringsAsFactors = FALSE
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
      grepl("^ub", mapping_df$ptm,ignore.case = TRUE)~"ubiquitin",
      grepl("^ac", mapping_df$ptm,ignore.case=TRUE)~"acetyl",
      grepl("^cle", mapping_df$ptm,ignore.case=TRUE)~"cleaved",
      TRUE~"other"),
    levels=c("total","phospho","methyl","ubiquitin","acetyl","cleaved","other"))

  #Validate names with HGNChelper
  valid_names <-unique(stats::na.omit(mapping_df$export_name))
  name_check  <-HGNChelper::checkGeneSymbols(valid_names)
  colnames(name_check)<-c('export_name','HGNC_check','suggestion')
  mapping_df  <-dplyr::left_join(mapping_df,name_check,by='export_name')
  mapping_df  <-mapping_df[,c(colnames(mapping_df)[-1],colnames(mapping_df)[1])]

  #Remember for HGNC_check: TRUE = valid, FALSE = corrected via suggestion,
  #NA = not found (usually histone marks)
  #Always use final_names downstream since it consolidates all three cases
  mapping_df$final_names<-ifelse(is.na(mapping_df$HGNC_check),NA
                                 ,ifelse(mapping_df$HGNC_check==FALSE
                                         ,mapping_df$suggestion
                                         ,mapping_df$export_name))
  #Detect invalid names
  if (any(is.na(mapping_df$final_names) & mapping_df$mapping_level !='histone_mark')) {
    message1<-"\nSome protein names are invalid and have no HGNC suggestion.
    Double-check mapping table.\n"
    if (strict) stop(message1,call.=FALSE)
    else warning(message1,call.=FALSE)
  }

  if (any(mapping_df$mapping_level=='histone_mark')) {
    warning("\nHistone marks were delivered as input. Not valid as HGNC input.\n")
  }

  if (any(mapping_df$HGNC_check==FALSE & !is.na(mapping_df$suggestion),na.rm=TRUE)){
    message2<-"\nSome protein names were auto-corrected using HGNC suggestions.
    Double-check mapping table.\n"
    if (strict) stop(message2,call.=FALSE)
    else warning(message2,call.=FALSE)
  }

  return(list(name_mapping   = mapping_df,
              original_names = unique(mapping_df$original),
              clean_names    = unique(stats::na.omit(mapping_df$final_names))
  ))
}

#' Get colors from pre-defined or user-defined color sets
#'
#' Returns a color vector or named list of color palettes for use in plots.
#' Includes defaults for up to 10 groups, with automatic interpolation
#' beyond that range via [grDevices::colorRampPalette()].
#'
#' @param n_colors Integer or `NULL`. Number of colors to return. If `NULL`,
#'   returns the full default vector or list depending on `mode`.
#' @param mode Character. One of `"list"` (default) or `"vector"`.
#'   `"vector"` returns a flat character vector of colors.
#'   `"list"` returns a named list where each element is a palette
#'   that corresponds to the number of groups (e.g. `list[["3"]]` for 3 groups).
#' @param color_list A named list of color vectors that overrides the default
#'   list palettes. Names should be the number of groups as character
#'   (e.g. `list("2" = c("blue", "red"))`). Ignored when `mode = "vector"`.
#' @param color_vector A character vector of colors to override the default
#'   vector palette. Ignored when `mode = "list"`.
#'
#' @return Either a character vector of colors (`mode = "vector"`) or a
#'   named list of color vectors (`mode = "list"`). If `n_colors` exceeds
#'   the available palette size, colors are interpolated automatically.
#'
#' @details
#' Default palettes are perceptually ordered from blue/purple (low) to
#' red/orange (high), suitable for ordered categorical variables. When
#' `n_colors` exceeds the largest pre-defined palette, colors are
#' interpolated using [grDevices::colorRampPalette()].
#'
#' Use `mode = "list"` when color order carries biological meaning —
#' typically when groups represent a gradient such as quantile-based
#' expression levels (e.g. low to high protein expression). Each palette
#' is hand-tuned to span blue/purple (low) through orange/red (high).
#' Use `mode = "vector"` for categorical groups where color order is
#' arbitrary (e.g. sample subtypes, gender) or when you prefer to assign
#' colors manually to ordered groups.
#'
#' @examples
#' #Get default vector of all colors
#' get_colors(mode="vector")
#'
#' #Get 5 colors as a vector
#' get_colors(n_colors=5,mode="vector")
#'
#' #Get optimized palette for 3 groups
#' get_colors(n_colors=3,mode="list")
#'
#' #Use custom colors
#' get_colors(n_colors=3,mode="vector",
#'              color_vector=c("navy","grey50","firebrick"))
#'
#' @export
get_colors<-function(n_colors=NULL,mode=c("list","vector")
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
    'red2','blue3','green4','purple2','darkorange','darkturquoise',
    'darkgoldenrod3','yellow2','deeppink2','lightsalmon2','darkcyan'
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











