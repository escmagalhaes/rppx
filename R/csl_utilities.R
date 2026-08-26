### CSL file configuration utilities ###

#' Function to register CSL entries
#' @export
#' @examples
#' ## Placeholder Example ##
load_csl_registry<-function() {
  c(
    bcj        = "blood-cancer-journal.csl",
    leukemia   = "nature.csl",
    ajh        = "american-medical-association-no-url.csl",
    aacr       = "american-association-for-cancer-research.csl",
    jco        = "journal-of-clinical-oncology.csl",
    lancet     = "the-lancet.csl",

    nature     = "nature.csl",
    science    = "science.csl",
    cell       = "cell.csl",
    plos       = "plos.csl",
    pnas       = "pnas.csl",
    bba        = "biochimica-et-biophysica-acta.csl",

    ama        = "american-medical-association-no-url.csl",
    apa        = "apa.csl",
    ieee       = "ieee.csl",
    chicago    = "chicago-author-date.csl",
    harvard    = "harvard-cite-them-right.csl",

    elsevier   = "elsevier-vancouver.csl",
    springer   = "springer-vancouver.csl"
  )
}

#' Function to sync CSL entries
#' @export
#' @examples
#' ## Placeholder Example ##
sync_csl_files<-function(csl_dir,registry) {

  #Create directory for csl files if necessary
  if (!dir.exists(csl_dir)) {
    dir.create(csl_dir,recursive=TRUE,showWarnings=FALSE)
  }

  #Download data from reliable url
  base_url_root <-"https://raw.githubusercontent.com/citation-style-language/styles/master/"
  base_url_dep  <-paste0(base_url_root,"dependent/")
  files         <-unique(unname(registry))
  failed        <- character(0)
  n_downloaded  <-0L

  for ( f in files) {
    dest<-file.path(csl_dir, f )
    needs_download<-!file.exists(dest) ||
      difftime(Sys.time(),file.mtime(dest),units="days") > 30

    if (!needs_download) next
    tmp     <-tempfile(fileext=".csl")
    try_url <-c(paste0(base_url_root, f ),paste0(base_url_dep, f ))
    ok      <- FALSE
    for (url in try_url) {
      ok<-tryCatch({
        download.file(url=url,destfile=tmp,quiet=TRUE)
        file.copy(tmp,dest,overwrite=TRUE)
        TRUE
      },error=function(e) FALSE)
      if (ok) break
    }
    unlink(tmp)

    if (ok) {
      n_downloaded<-n_downloaded+1L
    } else {
      warning("Failed CSL download: ", f )
      failed<-c(failed, f )
    }
  }

  #Return number of downloaded files
  if (n_downloaded > 0L) {
    cat("Downloaded",n_downloaded,"CSL file(s)\n")
  }

  #Surface missing files early with a clear message
  still_missing<-files[!file.exists(file.path(csl_dir,files))]
  if (length(still_missing) > 0) {
    stop("CSL sync incomplete. Missing files:\n  ",
         paste(still_missing,collapse="\n  "))
  }

  invisible(TRUE)
}

#' Function to get CSL entry
#' @export
#' @examples
#' ## Placeholder Example ##
get_csl<-function(cache_path,registry,journal="blood") {

  if (is.null(journal) || is.na(journal) || journal == "") {
    journal<-"blood"
  }
  file_name<-registry[[journal]]

  if (is.null(file_name)) {
    warning("Unknown journal:",journal,"Falling back to vancouver.csl")
    file_name<-"vancouver.csl"
  }
  path<-file.path(cache_path,file_name)

  #Safety check
  if (!file.exists(path)) {
    stop("CSL file not found: ", path,
         "\nRun sync_csl_files() to download missing styles.")
  }
  return(path)
}

#' Function wrapper to define CSL style
#' @export
#' @examples
#' ## Placeholder Example ##
define_csl_style<-function(journal,cache_path) {
  csl_registry<-load_csl_registry()
  tryCatch(
    sync_csl_files(csl_dir = cache_path, registry = csl_registry),
    error = function(e) stop(conditionMessage(e), call. = FALSE)
  )
  get_csl(journal=journal,cache_path=cache_path,registry=csl_registry)
}
