### REMEMBER TO OPEN ZOTERO SOFTWARE BEFORE RUNNING THIS WRAPPER ###

#' Function to wrap Zotero References pipeline
#' @export
#' @examples
#' ## Placeholder Example ##
process_references<-function(text,collection_name,user_id,api_key,
                             main_cache_path,zotero_version="3",
                             crossref_email=NULL) {

  #Create Zotero collection and local cache directory
  cat("\nCreating Zotero collection...\n")
  collection_key<-create_zotero_collection(
    name            = collection_name,
    user_id         = user_id,
    api_key         = api_key,
    main_cache_path = main_cache_path,
    zotero_version  = zotero_version
  )
  cat("Collection key:",collection_key,"\n")

  #Extract and validate identifiers with caching
  cat("\nExtracting and validating identifiers...\n")
  ref_result<-get_reference_table(
    text            = text,
    collection_name = collection_name,
    main_cache_path = main_cache_path,
    crossref_email  = crossref_email
  )

  ref_table <-ref_result$successful
  failed    <-ref_result$failed
  new_refs  <-ref_result$new

  #Stop if any identifiers failed
  if (nrow(failed) > 0) {
    return(list(
      table  = ref_table,   #successful entries that were cached
      failed = failed,      #rows that need fixing: pmid/doi, note explaining why
      new    = new_refs,
      text   = NULL,
      report = NULL
      ))
  }

  #Manual import of references to Zotero
  if (nrow(new_refs) > 0) {
  cat("\nImport references to Zotero\n\n")
  zotero_import_ids(new_refs)
  cat("\nPress ENTER when done...")
  readline()
  }

  #Fetch BBT citekeys and rewrite manuscript text
  cat("\nFetching BBT citekeys from Zotero...\n")
  citekeys<-get_citekeys(
    user_id        = user_id,
    api_key        = api_key,
    collection_key = collection_key,
    zotero_version = zotero_version
  )

  #Warning about missing citekeys
  missing_keys<-ref_table |>
    dplyr::filter(valid) |>
    join_citekeys(citekeys) |>
    dplyr::filter(is.na(citekey))

  if (nrow(missing_keys) > 0) {
    cat("\nAttention!",nrow(missing_keys),"item(s) have no BBT citekey yet.\n")
    cat("\nIs Better BibTeX installed and Zotero open and synced?\n")
    print(missing_keys |> dplyr::select(pmid,doi,title))
    cat("\nFix and re-run.\n")
  }

  cat("\nRewriting citations in manuscript text...\n")
  text_clean<-rewrite_citations(
    text     = text,
    table    = ref_table,
    citekeys = citekeys
  )

  #Export updated .bib from Zotero
  bib_path<-get_zotero_bib(
    collection_name = collection_name,
    collection_key  = collection_key,
    main_cache_path = main_cache_path
  )

  #Compile final report
  report<-ref_table |>
    join_citekeys(citekeys) |>
    dplyr::select(pmid,doi,title,first_author,year,citekey,valid,note)

  cat("\nALL DONE!\n")

  return(list(
    text     = text_clean,
    bib_path = bib_path,
    report   = report,
    table    = ref_table,
    citekeys = citekeys
  ))
}



