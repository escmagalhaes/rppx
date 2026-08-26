### Mining data from online repositories ###

#' Function to download, decompress and load studies from cBioPortal
#' @export
#' @examples
#' ## Placeholder Example ##
get_cbio_data<-function(study_id,main_cache_path,
                        cache=TRUE,verbose=TRUE,force=FALSE) {

  #Create study directory
  if(verbose) cat("\nCreating study folders for",study_id,"\n")
  study_directory<-file.path(main_cache_path,study_id)
  if (!dir.exists(study_directory)) {
    dir.create(study_directory,recursive=TRUE)
  }

  #Return Cached files if requested and existant
  if (cache) {
    cache_version<-.cbio_cache_version
    cache_file<-file.path(study_directory,paste0(study_id, "_",cache_version,".qs2"))
    if (!force && file.exists(cache_file)) {
      if (verbose) cat("\nLoading cached study:",study_id,"\n")
      return(qs2::qs_read(cache_file))
    }
  }

  #Create temporary extraction sub-directory
  extract_dir<-file.path(study_directory,"extracted_files")
  if (!dir.exists(extract_dir)) {
    dir.create(extract_dir,recursive=TRUE)
  }

  #Download .tar.gz files (automatically skipped if file already exists)
  if(verbose) cat("\nDownloading",study_id,"\n")
  compressed_study_file<-cBioPortalData::downloadStudy(
    cancer_study_id = study_id,
    use_cache       = study_directory)

  #Decompress studies
  if(verbose) cat("\nDecompressing",study_id,"\n")
  study_files<-cBioPortalData::untarStudy(
    cancer_study_file = compressed_study_file,
    exdir             = extract_dir)

  #Load studies
  if(verbose) cat("\nLoading",study_id,"\n")
  loaded_study<-cBioPortalData::loadStudy(study_files,cleanup=FALSE)

  #Safe saving full study as qs2 file (avoid mid-workflow crash issues) if cache=TRUE
  if (cache) {
    if (verbose) cat("\nCaching study:",study_id,"\n")
    temp_cache_file<-paste0(cache_file,".tmp")
    qs2::qs_save(loaded_study,temp_cache_file)
    successful_workflow<-file.rename(
      from = temp_cache_file,
      to   = cache_file)

    if (!successful_workflow)
      stop("Failed to finalize cache file")

    if (verbose) cat("\nCleaning up raw files for",study_id,"\n")
    unlink(extract_dir,recursive=TRUE)
    unlink(compressed_study_file)
  } else {
    if (verbose) cat("\nNo caching requested. Extracted files retained in:",extract_dir,"\n")
  }
  return(loaded_study)
}

#' Function to extract relevant data from cBioPortal and create a single dataset
#' @export
#' @examples
#' ## Placeholder Example ##
get_study_data<-function(data,study_id,main_cache_path,
                         assay_names=NULL,assay_labels=NULL,
                         patient_id_col="PATIENT_ID",sample_id_col="SAMPLE_ID",
                         verbose=TRUE,cache=TRUE,force=FALSE) {

  if (cache) {
    if (is.null(main_cache_path))
      stop("main_cache_path must be provided when cache = TRUE")
    cache_path <-file.path(main_cache_path,study_id)
    if (!dir.exists(cache_path)) { dir.create(cache_path,recursive=TRUE) }
    cache_key  <-digest::digest(sort(assay_names),algo="xxhash32")
    cache_file <-file.path(cache_path,paste0(study_id,"_",cache_key,".qs2"))

    if (!force && file.exists(cache_file)) {
      if (verbose) cat("\nLoading cached extraction for", study_id, "\n")
      return(qs2::qs_read(cache_file))
    }
  }

  #Safety checks
  if(!is(data,"MultiAssayExperiment"))
    stop("data is not a 'MultiAssayExperiment'. Double check input")
  if(!is.null(assay_names)){
    missing_assay_names<-setdiff(assay_names,names(data))
    if(length(missing_assay_names) > 0)
      stop("These assay_names were not found in data:",missing_assay_names,
           "- Double check 'ExperimentList' item names")
  } else {
    warning("\nArgument assay_names not specified. Returning full data\n")
    return(data)
  }
  if(!is.null(assay_labels) && !is.null(assay_names)){
    if (length(assay_names) != length(assay_labels))
      stop("The number of assay_labels (",length(assay_labels),
           ") does not match the number of assay_names (",length(assay_names),
           "). Please make sure all assay_names are appropriately labelled")
    if (length(unique(assay_labels)) != length(assay_labels))
      stop("assay_labels must be unique because duplicate labels would cause column naming issues")
  }
  if(is.null(assay_labels) && !is.null(assay_names)){
    warning("\nArgument assay_labels not specified. Generic labels will be applied")
    assay_labels<-paste0("assay",seq_along(assay_names))
  }

  #Get sample and patient IDs mapping
  data_map<-MultiAssayExperiment::sampleMap(data) |> as_tibble()

  #Extract and transpose expression data (features x samples) and add patient IDs
  expr_list<-purrr::map2(assay_names,assay_labels,function(types,labs){
    if (!is(data[[types]],"RaggedExperiment")) {
      raw_mtx<-MultiAssayExperiment::assay(data[[types]])
    } else {
      raw_mtx<-MultiAssayExperiment::assay(data[[types]],i="Variant_Classification")
    }
    #Adjust data matrix
    raw_mtx_adj<-raw_mtx |> as.data.frame() |> t() |> as.data.frame() |>
      rename_with(~paste(labs,.x,sep="_")) |>
      tibble::rownames_to_column(var = "colname") |> as_tibble()

    #Filter data_map according to assay
    data_map_filt<-data_map |> dplyr::filter(assay == types) |> select(-assay)

    #Add Patient IDs to data matrix
    final_mtx<-right_join(data_map_filt,raw_mtx_adj,by=c("colname"="colname")) |>
      dplyr::rename(!!patient_id_col := primary, !!sample_id_col := colname)
  })
  names(expr_list)<-assay_names

  #Extract Clinical Data
  clinical_data<-MultiAssayExperiment::colData(data) |> as_tibble()

  result<-c(list(clinical_data = clinical_data),expr_list)

  if (cache) {
    if (verbose) cat("\nSaving cached extraction for",study_id,"\n")
    temp_cache_file<-paste0(cache_file,".tmp")
    qs2::qs_save(result,temp_cache_file)
    successful_workflow<-file.rename(
      from = temp_cache_file,
      to   = cache_file)

    if (!successful_workflow)
      stop("Failed to finalize cache file")
  }

  return(result)
}
