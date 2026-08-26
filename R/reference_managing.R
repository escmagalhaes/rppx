### References configuration using Zotero ###

#' Function to pass zotero headers
#' @export
#' @examples
#' ## Placeholder Example ##
zotero_headers<-function(api_key,zotero_version="3") list(
  `Zotero-API-Key`     = api_key,
  `Content-Type`       = "application/json",
  `Zotero-API-Version` = zotero_version)

#' Function to create Zotero collection
#' @export
#' @examples
#' ## Placeholder Example ##
create_zotero_collection<-function(name,user_id,api_key,main_cache_path,
                                   zotero_version="3",parent_key=NULL) {

  #POST temporary collection and delete after guards
  temp_resp<-request(paste0("https://api.zotero.org/users/", user_id ,"/collections")) |>
    req_method("POST") |>
    req_headers(!!!zotero_headers(api_key,zotero_version)) |>
    req_body_json(list(list(name = "temporary"))) |>
    req_perform()
  temp_res <-resp_body_json(temp_resp)
  temp_key <-temp_res$success$`0`

  if (is.null(temp_key)) {
    stop("\nTemporary collection creation failed. Check if Zotero POST response shape changed.\nRaw response: ",
         resp_body_string(temp_resp),"\n")
  }

  #GET collections
  collections<-request(paste0("https://api.zotero.org/users/", user_id ,"/collections")) |>
    req_method("GET") |>
    req_headers(!!!zotero_headers(api_key,zotero_version)) |>
    req_perform() |>
    resp_body_json()

  collection_ids <- tryCatch({
    keys  <-vapply(collections, \(x) x$data$key,  character(1))
    names <-vapply(collections, \(x) x$data$name, character(1))
    stats::setNames(keys, names)
  }, error = function(e) {
    stop("\nCould not parse collection list. Check if Zotero GET response shape changed.\n")
  })

  #DELETE temporary collection
  temp_version<-temp_res$successful$`0`$version %||% 1
  request(
    paste0("https://api.zotero.org/users/", user_id ,"/collections/", temp_key)) |>
    req_method("DELETE") |>
    req_headers(
      `Zotero-API-Key`     = api_key,
      `Zotero-API-Version` = zotero_version,
      `If-Unmodified-Since-Version` = as.character(temp_version)) |>
    req_perform()

  #Return early if already exists
  if(name %in% names(collection_ids)) {
    cat("\nCollection already exists\n")
    return(collection_ids[[name]])
  }

  #Create collection if non-existant
  body<-list(list(name=name,parentCollection=if (is.null(parent_key)) FALSE else parent_key))
  resp<-request(paste0("https://api.zotero.org/users/", user_id ,"/collections")) |>
    req_method("POST") |>
    req_headers(!!!zotero_headers(api_key,zotero_version)) |>
    req_body_json(body) |>
    req_perform()
  result  <-resp_body_json(resp)
  new_key <-result$success$`0`

  if (is.null(new_key)) {
    stop("\nCollection creation failed.\nRaw response: ",resp_body_string(resp),"\n")
  }

  #Create local cache directory
  collection_dir<-file.path(main_cache_path,name)
  if (!dir.exists(collection_dir)) {
    dir.create(collection_dir,recursive=TRUE,showWarnings = FALSE)
  }
  cat("\nCollection successfully created:",name,"=",new_key,"\n")
  cat("Local cache directory:",collection_dir,"\n")

  return(new_key)
}

#' Function to extract and validate identifiers from text
#' @export
#' @examples
#' ## Placeholder Example ##
extract_ids<-function(text,input_pattern,output_pattern,known_ids=character(0)) {
  stringr::str_extract_all(text,input_pattern)[[1]] |>
    stringr::str_extract_all(output_pattern) |> unlist() |>
    stringr::str_remove("\\.$") |> unique() |> setdiff(known_ids)
}

#' Function to validate PMIDs
#' @export
#' @examples
#' ## Placeholder Example ##
validate_pmid<-function(pmid) {

  rec<-tryCatch(
    rentrez::entrez_summary(db="pubmed",id=pmid),
    error=function(e) NULL)

  if (is.null(rec)) {
    return(tibble::tibble(
      pmid         = pmid,
      doi          = NA_character_,
      title        = NA_character_,
      first_author = NA_character_,
      year         = NA_character_,
      valid        = FALSE,
      note         = "NCBI request failed"))
  }

  if (is.null(rec$title) || rec$title=="") {
    return(tibble::tibble(
      pmid         = pmid,
      doi          = NA_character_,
      title        = NA_character_,
      first_author = NA_character_,
      year         = NA_character_,
      valid        = FALSE,
      note         = "Not found in PubMed"))
  }

  first_author <-tryCatch(rec$authors$name[1],error=function(e) NA_character_) %||% NA_character_
  year         <-stringr::str_extract(rec$pubdate %||% "", "\\d{4}")
  doi          <-rec$elocationid |> stringr::str_extract("(?<=doi: )\\S+") %||% NA_character_

  return(tibble::tibble(
    pmid         = pmid,
    doi          = doi,
    title        = rec$title,
    first_author = first_author,
    year         = year,
    valid        = TRUE,
    note         = "OK"))
}

#' Function to validate DOIs
#' @export
#' @examples
#' ## Placeholder Example ##
validate_doi<-function(doi,crossref_email=NULL) {

  #Resolve crossref email if passed (speeds up Crossref request)
  if (!is.null(crossref_email)) options(rcrossref_email = crossref_email)

  rec<-tryCatch(
    rcrossref::cr_works(dois=doi)$data,
    error=function(e) NULL)

  if (is.null(rec) || nrow(rec)==0) {
    return(tibble::tibble(
      pmid         = NA_character_,
      doi          = doi,
      title        = NA_character_,
      first_author = NA_character_,
      year         = NA_character_,
      valid        = FALSE,
      note         = "Not found in Crossref"))
  }

  title        <-rec$title %||% NA_character_
  first_author <-tryCatch(rec$author[[1]]$family[1],error=function(e) NA_character_) %||% NA_character_
  year         <-as.character(rec$issued %||% NA_character_)

  if (is.null(title) || is.na(title)) {
    return(tibble::tibble(
      pmid         = NA_character_,
      doi          = doi,
      title        = NA_character_,
      first_author = NA_character_,
      year         = NA_character_,
      valid        = FALSE,
      note         = "Not found in Crossref"))
  }

  #Fetch PMID via id_converter if available
  pmid<-tryCatch({
    res<-rcrossref::id_converter(doi)
    as.character(res$records$pmid %||% NA_character_)
  },error=function(e) NA_character_)

  return(tibble::tibble(
    pmid         = pmid,
    doi          = doi,
    title        = title,
    first_author = first_author,
    year         = year,
    valid        = TRUE,
    note         = "OK"))
}

#' Function to extract and validate all IDs
#' @export
#' @examples
#' ## Placeholder Example ##
get_reference_table<-function(text,collection_name,main_cache_path,
                              crossref_email     = NULL,
                              pmid_field_pattern = "\\(\\s*PMID\\s*:[^)]*\\)",
                              pmid_clean_pattern = "\\d{5,9}",
                              doi_field_pattern = "\\(\\s*DOI\\s*:[^()]*(?:\\([^()]*\\)[^()]*)*\\)",
                              doi_clean_pattern = "10\\.\\d{4,9}/(?:[^\\s,:;()]*(?:\\([^()]*\\))?)+(?=\\.?(?:\\s*[,:;)]|\\)$))") {

  #Define cache path and file
  cache_file<-file.path(main_cache_path,collection_name,"zotero_cache.qs2")

  #Define empty ref_table
  empty_ref_table<-tibble::tibble(
    pmid         = character(),
    doi          = character(),
    title        = character(),
    first_author = character(),
    year         = character(),
    valid        = logical(),
    note         = character())

  #Load existing cache or initialise empty
  cache<-if (file.exists(cache_file)) {
    qs2::qs_read(cache_file)
  } else {
    empty_ref_table
  }

  #Extract only new identifiers not already in cache
  #PMIDs
  known_pmids <-cache$pmid[!is.na(cache$pmid)]
  new_pmids   <-extract_ids(
    text           = text,
    input_pattern  = pmid_field_pattern,
    output_pattern = pmid_clean_pattern,
    known_ids      = known_pmids)
  new_pmid_results<- if (length(new_pmids) > 0) {
    purrr::map_dfr(new_pmids,validate_pmid) |>
      dplyr::distinct(pmid,.keep_all=TRUE) #Deduplicate PMIDs
  } else {
    empty_ref_table
  }

  #DOIs
  #Skip DOIs already captured via PMID
  known_dois       <-cache$doi[!is.na(cache$doi)]
  known_dois_pmids <-c(known_dois,new_pmid_results$doi[!is.na(new_pmid_results$doi)])
  new_dois         <-extract_ids(
    text           = text,
    input_pattern  = doi_field_pattern,
    output_pattern = doi_clean_pattern,
    known_ids      = known_dois_pmids)
  new_doi_results<- if (length(new_dois) > 0) {
      purrr::map_dfr(new_dois,validate_doi,crossref_email=crossref_email) |>
      dplyr::distinct(doi,.keep_all=TRUE) #Deduplicate DOIs
    } else {
      empty_ref_table
    }

  #Merge results
  new_results<-dplyr::bind_rows(new_pmid_results,new_doi_results)

  #Merge with cache and overwrite cache file (cache only valid entries)
  if (nrow(new_results) > 0) {
    new_results_valid<-new_results |> dplyr::filter(valid)
    new_cache<-dplyr::bind_rows(cache,new_results_valid) |>
      (\(x) dplyr::bind_rows(
        x |> dplyr::filter(!is.na(doi)) |> dplyr::distinct(doi,.keep_all=TRUE),
        x |> dplyr::filter( is.na(doi)) |> dplyr::distinct(pmid,.keep_all=TRUE)
      ))()
    qs2::qs_save(new_cache,cache_file)
    cat("\nSaved",sum(new_results$valid),"valid new entry(ies) to",cache_file,
        " (", nrow(new_results),"identifier(s) evaluated).\n")
  } else {
    cat("\nNo new identifiers. Using cached results.\n")
    return(list(successful=cache,failed=empty_ref_table,new=empty_ref_table))
  }

  #Report
  if (nrow(new_results) > 0) {
    failed<-new_results |> dplyr::filter(!valid)
    if (nrow(failed) > 0) {
      cat("\nAttention!",nrow(failed),"identifier(s) failed validation:\n")
      print(failed |> dplyr::select(pmid,doi,note))
      cat("\nFix these in the source text and re-run. Valid entries are cached and won't be re-checked.\n")
    } else {
      cat("\nValidated",sum(new_results$valid),"of",nrow(new_results),"new identifier(s).\n")
    }
  }
  return(list(successful=new_cache,failed=failed,new=new_results_valid))
}

#' Function to print IDs for Zotero import
#' @export
#' @examples
#' ## Placeholder Example ##
zotero_import_ids<-function(table) {

  valid_dois<-table |>
    dplyr::filter(valid,!is.na(doi),doi != "") |>
    dplyr::pull(doi)

  if (length(valid_dois)==0) {
    cat("\nValid DOIs already cached and added to Zotero!\n")
  } else {
    cat(valid_dois,sep=" ")
    cat("\n\nPaste the DOIs above into Zotero's Add Item by Identifier (magic wand).\n")
    cat("Make sure all items appear under '",collection_name,"'.\n",sep="")
    cat("Deduplicate references if needed:\n")
    cat("Click on Duplicate Items (left-side panel).\n")
    cat("Then click on duplicates (Zotero selects all copies of that reference).\n")
    cat("Finally click Merge (top right panel). Do this for every duplicated item.\n")
  }

  no_doi_pmids<-table |>
    dplyr::filter(valid,is.na(doi) | doi=="") |>
    dplyr::pull(pmid)

  if (length(no_doi_pmids) > 0) {
    cat("\nDo not proceed just yet! Below are valid refs that have no DOI in NCBI's esummary (old papers)!\n")
    cat(no_doi_pmids,sep=" ")
    cat("\nPaste the PMIDs above into Zotero's Add Item by Identifier (magic wand).\n")
    cat("Make sure all items appear under '",collection_name,"'.\n",sep="")
    cat("Remember to deduplicate references if needed.\n")
  }
}

#' Function to extract citation keys from each reference (with pagination)
#' @export
#' @examples
#' ## Placeholder Example ##
get_citekeys<-function(user_id,api_key,collection_key,zotero_version="3") {

  all_items <-list()
  start     <-0
  limit     <-100

  repeat {
    batch<-request(paste0("https://api.zotero.org/users/",user_id,
                          "/collections/",collection_key,"/items")) |>
      req_method("GET") |>
      req_headers(!!!zotero_headers(api_key,zotero_version)) |>
      req_url_query(format="json",limit=limit,start=start) |>
      req_perform() |>
      resp_body_json()

    all_items<-c(all_items,batch)

    #If fewer results than limit returned, we have everything
    if (length(batch)==0 || length(batch) < limit) break
    start<-start+limit
  }

  purrr::map_dfr(all_items, function(it) {
    data<-it$data
    if (data$itemType=="attachment") return(NULL)
    tibble::tibble(
      doi     = data$DOI         %||% NA_character_,
      pmid    = dplyr::na_if(data$PMID %||% NA_character_, ""),
      citekey = data$citationKey %||% NA_character_)
  }) |>
    dplyr::filter(!is.na(citekey))
}

#' Function to join citekeys properly
#' @export
#' @examples
#' ## Placeholder Example ##
join_citekeys<-function(table,citekeys) {

  #Fix empty string PMIDs in citekeys
  citekeys<-citekeys |> dplyr::mutate(pmid=dplyr::na_if(pmid,""))

  #Join on DOI
  table<-table |> dplyr::left_join(citekeys |> dplyr::filter(!is.na(doi)) |>
        dplyr::select(doi,citekey_doi=citekey),by="doi")

  #Join on PMID
  table<-table |> dplyr::left_join(citekeys |> dplyr::filter(!is.na(pmid)) |>
        dplyr::select(pmid,citekey_pmid=citekey),by="pmid")

  #Prefer DOI match; otherwise use PMID match
  table |> dplyr::mutate(citekey=dplyr::coalesce(citekey_doi,citekey_pmid)) |>
    dplyr::select(-citekey_doi,-citekey_pmid)
}

#' Function to replace text with Zotero citation keys (silently skips already converted citations)
#' @export
#' @examples
#' ## Placeholder Example ##
rewrite_citations<-function(text,table,citekeys) {

  lookup<-table |> dplyr::filter(valid) |>
    join_citekeys(citekeys) |>
    dplyr::filter(!is.na(citekey))

  #Replace identifiers with citekeys
  for ( i in seq_len(nrow(lookup))) {
    row<-lookup[ i, ]

    #Use whichever identifier actually appears in the text
    if (!is.na(row$pmid)) {
      text<-stringr::str_replace_all(text,stringr::fixed(row$pmid),paste0("@",row$citekey))
    }
    if (!is.na(row$doi)) {
      text<-stringr::str_replace_all(text,stringr::fixed(row$doi),paste0("@",row$citekey))
    }
  }

  #Convert (PMID: @key1; @key2) and (DOI: @key) blocks to [@key1;@key2]
  blocks<-stringr::str_extract_all(text,"\\(\\s*(?:PMID|DOI)\\s*:[^)]*\\)")[[1]]
  for (blk in blocks) {

    #Skip blocks that still contain unresolved raw identifiers
    has_raw_pmid <-stringr::str_detect(blk,"(?<!@)\\b\\d{5,9}\\b")
    has_raw_doi  <-stringr::str_detect(blk,"(?<!@)\\b10\\.\\d{4,}/\\b")
    if (has_raw_pmid || has_raw_doi) next

    keys<-stringr::str_extract_all(blk,"@[\\w.:-]+")[[1]]
    if (length(keys)==0) next

    replacement<-paste0("[",paste(keys,collapse=";"),"]")
    text<-stringr::str_replace_all(text,stringr::fixed(blk),replacement)
  }

  #Warning for unresolved citations
  remaining<-stringr::str_extract_all(text,"\\(\\s*(?:PMID|DOI)\\s*:[^)]*\\)")[[1]]
  remaining<-remaining[
    stringr::str_detect(remaining,"(?<!@)\\b\\d{5,9}\\b") |
      stringr::str_detect(remaining,"(?<!@)\\b10\\.\\d{4,}/\\b")]
  if (length(remaining) > 0) {
    warning(length(remaining)," citation block(s) could not be fully converted.")
    cat("\nReview these citations:")
    cat(remaining,sep="\n")
  }
  return(text)
}

#' Function to get references.bib from Zotero
#' @export
#' @examples
#' ## Placeholder Example ##
get_zotero_bib<-function(collection_name,collection_key,main_cache_path) {

  #Fecth bib refs using collection_key
  url<-paste0(
    "http://localhost:23119/better-bibtex/export/collection?/1/",collection_key,".bib")

  #Define path to save
  bib_path<-file.path(main_cache_path,collection_name,"references.bib")

  tryCatch({
    utils::download.file(url=url,destfile=bib_path,quiet=TRUE)
    cat("\nBibliography exported to:\n",bib_path,"\n")
  },error=function(e) {
    stop("\nCould not reach Zotero. Remember to open the software!\n",conditionMessage(e),call.=FALSE)
  })
  invisible(bib_path)
}

