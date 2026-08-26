### Adjust manuscript for rendering using quarto ###

#' Function to replace text with Quarto keys for Figures and Tables
#' @export
#' @examples
#' ## Placeholder Example ##
rewrite_tag<-function(text,tag,prefix) {
  stringr::str_replace_all(text,paste0("\\(", tag , ":\\s*([^)]+)\\)"),
    function(m) {
      id<-stringr::str_extract(m,paste0("(?<=", tag , ":\\s{0,10})[^)]+")) |>
        stringr::str_trim() |>
        stringr::str_replace_all("\\s+", "-") |>
        tolower()
      paste0(prefix,id,"]")
      })
}

#' Wrapper funtion to raw text with all Figures and Tables Quarto keys
#' @export
#' @examples
#' ## Placeholder Example ##
rewrite_figure_table_keys <- function(
    text,
    mfig_tag = "FIGURE MAIN",
    sfig_tag = "FIGURE SUP",
    mtbl_tag = "TABLE MAIN",
    stbl_tag = "TABLE SUP") {
  text |>
    rewrite_tag(tag = mfig_tag, prefix = "[@mfig-MAIN-FIG-") |>
    rewrite_tag(tag = sfig_tag, prefix = "[@sfig-SUP-FIG-") |>
    rewrite_tag(tag = mtbl_tag, prefix = "[@mtbl-MAIN-TAB-") |>
    rewrite_tag(tag = stbl_tag, prefix = "[@stbl-SUP-TAB-")
}

#' Function to add source_type to mapping table
#' @export
#' @examples
#' ## Placeholder Example ##
add_source_type<-function(lookup_table) {
  lookup_table |> dplyr::mutate(
    source_type = dplyr::case_when(
      stringr::str_ends(key_name,"-dataset") ~ "dataset",
      key_type %in% c("mfig", "sfig")        ~ "figure",
      key_type %in% c("mtbl", "stbl")        ~ "table",
      .default = NA_character_))
}

#' Function to add file_paths to mapping table
#' @export
#' @examples
#' ## Placeholder Example ##
add_file_paths<-function(lookup_table,figure_export_path,
                         table_export_path,dataset_export_path) {

  lookup_table$file_path<-purrr::map(
    seq_len(nrow(lookup_table)),
    function( i ) {
      row<-lookup_table[ i , ]
      root<-switch(
        row$source_type,
        figure  = figure_export_path,
        table   = table_export_path,
        dataset = dataset_export_path,
        stop("Unknown source_type: ",row$source_type))

      folder <-file.path(root,row$key_name)
      files  <-list.files(folder,full.names=TRUE)

      #Do not get paths related to intermediary files
      base<-basename(files)
      keep<-stringr::str_detect(base,paste0(
        "^",stringr::fixed(row$key_name),"(?:\\.[^.]+|_Page_[0-9]+\\.png)$"))
      files<-files[keep]

      if (length(files) == 0) {
        stop("No files found for: ",row$key)
      }
      return(files)
    })
  return(lookup_table)
}

#' Function to add yaml header to .qmd file
#' @export
#' @examples
#' ## Placeholder Example ##
add_yaml_header<-function(bibliography,csl_style,format="docx",
                          include_refs=TRUE,include_figures = FALSE,
                          dataset_export_path=NULL,table_export_path=NULL,
                          figure_export_path=NULL) {

  params_block<-c(
    "params:",
    "  include_figures: false",
    if (!is.null(figure_export_path))
      sprintf("  figure_export_path: '%s'",figure_export_path),
    if (!is.null(table_export_path))
      sprintf("  table_export_path: '%s'", table_export_path),
    if (!is.null(dataset_export_path))
      sprintf("  dataset_export_path: '%s'",dataset_export_path)
  )

  crossref_block <- if (include_refs) {
    c("crossref:",
      "  fig-title: Figure",
      "  tbl-title: Table",
      "  title-delim: ''",
      "  custom:",
      "    - kind: float",
      "      key: mfig",
      "      latex-env: mfig",
      '      reference-prefix: "Figure"',
      "      space-before-numbering: false",
      "    - kind: float",
      "      key: mtbl",
      "      latex-env: mtbl",
      '      reference-prefix: "Table"',
      "      space-before-numbering: false",
      "    - kind: float",
      "      key: sfig",
      "      latex-env: sfig",
      '      reference-prefix: "Supplemental Figure S"',
      "      space-before-numbering: false",
      "    - kind: float",
      "      key: stbl",
      "      latex-env: stbl",
      '      reference-prefix: "Supplemental Table S"',
      "      space-before-numbering: false")
  } else NULL

  c("---",
    sprintf("format: %s", format),
    sprintf("bibliography: '%s'", bibliography),
    sprintf("csl: '%s'", csl_style),
    sprintf("citeproc: %s", tolower(as.character(include_refs))),
    "link-bibliography: false",
    "suppress-bibliography-fields: [doi, url, isbn, issn]",
    params_block,
    crossref_block,
    "---", "")
}

#' Function to extract keys from text in order of appearance
#' @export
#' @examples
#' ## Placeholder Example ##
extract_float_order<-function(manuscript_path,lookup_table) {

  #Extract text from file
  qmd_text<-paste(readLines(manuscript_path),collapse="\n")

  #Extract all keys in order of appearance
  all_keys<-stringr::str_extract_all(qmd_text,"(?<=\\[@)[\\w-]+")[[1]] |> unique()

  #Filter to only keys in lookup_table and preserve order
  all_keys<-all_keys[all_keys %in% lookup_table$key]

  #Reorder lookup by appearance
  lookup_table |> dplyr::filter(key %in% all_keys) |>
    dplyr::arrange(match(key,all_keys))
}

#' Function to generate float chunks to embed within text
#' @export
#' @examples
#' ## Placeholder Example ##
make_float_chunks<-function(lookup_table) {

  placeholder_fig<-file.path(biocache_root,"Quarto","placeholder.png")

  purrr::map_chr(seq_len(nrow(lookup_table)), function( i ) {
    row        <-lookup_table[ i , ]
    key        <-row$key
    legend     <-row$legend
    file_path  <-row$file_path[[1]] #file_path is a list

    #Adjust caption  and legend to make it safely readable in R
    caption_safe <-stringr::str_replace_all(row$caption,"'","\\\\'")

    #Create an object (graphics_call) for passing placeholder_fig correctly
    if (row$source_type=="dataset") {

      #Datasets always use placeholder_fig
      graphics_call<-paste0("knitr::include_graphics('", placeholder_fig ,"')\n")

    } else {

      #For Figures, add the Figures with paths
      graphics_call<-paste0(
        "if (params$include_figures) {\n",
        "  knitr::include_graphics(c(",
        paste0("'", file_path ,"'",collapse=", "),
        "))\n",
        "} else {\n",
        "  knitr::include_graphics('", placeholder_fig ,"')\n",
        "}\n"
      )
    }

    #Make the float chunks and add legend below in plain markdown
    paste0(
      "\n::: {#", key , "}\n\n",

      #Figure (inside float)
      "```{r}\n",
      "#| label: ", key , "\n",
      "#| echo: false\n",
      graphics_call ,
      "```\n\n",

      #Caption (inside float)
      "```{r}\n",
      "#| echo: false\n",
      "#| results: asis\n",
      "knitr::asis_output('", caption_safe, "')\n",
      "```\n\n",

      ":::\n",

      #Legend (outside float)
      if (!is.na(legend) && nzchar(legend)) {
        paste0(legend,"\n\n")
      } else {
        ""
      }
    )
  }) |> paste(collapse="\n")
}

#' Function to select file types in lookup table to export according to rendering mode
#' @export
#' @examples
#' ## Placeholder Example ##
filter_lookup_table<-function(lookup_table,mode=c("draft","review","submission")) {

  mode<-match.arg(mode)

  #Define extensions according to mode
  figure_ext     <-if (mode=="review") "png" else "pdf"
  main_table_ext <-if (mode=="review") "png" else "pdf"
  sup_table_ext  <-if (mode=="review") "png" else "xlsx"

  #Add selected_ext to lookup_table
  lookup_table_ext<-lookup_table |> dplyr::mutate(selected_ext = case_when(
      source_type == "figure" ~ figure_ext,
      source_type == "table" & file_type == "MAIN" ~ main_table_ext,
      source_type == "table" & file_type == "SUP"  ~ sup_table_ext,
      source_type == "dataset" ~ "xlsx",
      .default = NA_character_))

  #Filter lookup_table
  lookup_table_filt<-lookup_table_ext |> dplyr::mutate(
    file_path = purrr::pmap(
      list(file_path,selected_ext,source_type,file_type),
      function(paths,selected_ext,source_type,file_type) {
        #Extract extensions from file_paths
        ext<-tools::file_ext(paths)
        selected<-paths[ext==selected_ext]

        #For gtsummary tables without xlsx files, fallback to pdf
        if (length(selected)==0 && source_type=="table" && file_type=="SUP") {
            selected<-paths[ext=="pdf"]
          }
          selected
        }))

  #Update selected_ext so it reflects the actual selected file
  lookup_table_filt<-lookup_table_filt |> dplyr::mutate(
    selected_ext=purrr::map_chr(file_path,~tools::file_ext(.x[[1]])))

  #Safety check
  if (any(lengths(lookup_table_filt$file_path) == 0)) {
    stop("Some lookup table entries have no files matching the selected extension.")
  }
  return(lookup_table_filt)
}

#' Function to render manuscript file
#' @export
#' @examples
#' ## Placeholder Example ##
render_manuscript<-function(manuscript_path,manuscript_name="manuscript",
                              mode=c("draft","review","submission")) {
  mode <- match.arg(mode)

  include_figures <-mode == "review"
  include_refs    <-mode != "draft"

  output_file<-sprintf("%s_%s_%s.docx",
                       format(Sys.Date(),"%Y.%m.%d"),
                       manuscript_name,mode)

  quarto::quarto_render(
    input          = manuscript_path,
    output_file    = output_file,
    execute_params = list(include_figures = include_figures),
    cache          = FALSE
  )
}

#' Function to prepare submission folder
#' @export
#' @examples
#' ## Placeholder Example ##
prepare_submission<-function(manuscript_path,lookup_table,output_path,
                             sup_naming="Supplemental") {

  #Create output directory if non-existant
  if (!dir.exists(output_path)) { dir.create(output_path,recursive=TRUE) }

  #Copy and rename files
  for (i in seq_len(nrow(lookup_table))) {
    row  <-lookup_table[ i, ]
    src  <-row$file_path[[1]]
    dest <-file.path(output_path,row$label)
    file.copy( src, dest ,overwrite=TRUE)
    if (stringr::str_ends(src,"\\.xlsx")) {

      #Rewrite A1 cells in xlsx sheets
      wb         <-openxlsx2::wb_load(dest)
      sheets     <-openxlsx2::wb_get_sheet_names(wb)
      identifier <-row$key_name
      sub_label  <-paste0(sup_naming," Table S",row$label_n)

      for (s in seq_along(sheets)) {
        #Extract sheet header (cell A1)
        a1<-as.character(openxlsx2::wb_read(wb, sheet = s , col_names=FALSE)[1,1])

        if (!is.na(identifier) && stringr::str_starts(a1,identifier)) {

          #Extract panel and caption if present at end of identifier
          head_noid <-stringr::str_remove(a1,paste0("^",stringr::fixed(identifier)))
          panel     <-stringr::str_extract(head_noid,"^[A-Z](?=\\.)")
          caption   <-stringr::str_remove(head_noid,"^[A-Z]?\\. ?")

        } else {
          #No identifier available
          panel   <-NA_character_
          caption <-a1
        }
        new_a1<-paste0(
          sub_label,
          if (!is.na(panel)) panel else "",
          ". ",
          caption
          )

        #Replace a1 with new_a1
        wb$add_data(sheet=sheets[ s ],x=new_a1,start_col=1,start_row=1,col_names=FALSE)
      }
      #Save workbook
      openxlsx2::wb_save(wb, dest , overwrite=TRUE)
    }
    cat( row$key ,"renamed to", row$label,"\n")
  }

  # Export lookup table #
  #Create directory for exporting lookup_table if non-existant
  lookup_table_path<-file.path(output_path,'Submission_mapping_table')
  if (!dir.exists(lookup_table_path)) { dir.create(lookup_table_path,recursive=TRUE) }
  export_tabs(
    df_list   = list(submission_map = lookup_table),
    filepath  = lookup_table_path,
    filename  = "Submission_mapping_table",
    tab_title = ""
  )

  invisible(lookup_table)
}

#' Function to wrap manuscript exporting
#' @export
#' @examples
#' ## Placeholder Example ##
manuscript_export<-function(manuscript_path,ref_bib_path,csl_cache_path,
                            lookup_table,dataset_export_path=NULL,
                            table_export_path=NULL,figure_export_path=NULL,
                            output_path=NULL,ref_header="**REFERENCES**",
                            figure_table_header="FIGURES AND TABLES",
                            sup_naming="Supplemental",journal_ref="ajh",
                            manuscript_name="manuscript",format="docx",
                            mode=c("draft","review","submission")){

  mode            <-match.arg(mode)
  include_refs    <-mode != "draft"
  include_figures <-mode == "review"

  #Define CSL style
  csl_style<-define_csl_style(journal = journal_ref,cache_path = csl_cache_path)

  #Add yaml header
  yaml_header<-add_yaml_header(
    bibliography        = ref_bib_path,
    csl_style           = csl_style,
    format              = format,
    include_refs        = include_refs,
    include_figures     = include_figures,
    figure_export_path  = figure_export_path,
    table_export_path   = table_export_path,
    dataset_export_path = dataset_export_path
  )

  #Read qmd, replace yaml header, write back
  qmd_lines     <-readLines(manuscript_path)
  yaml_end      <-which(qmd_lines=="---")[2]
  body_start    <-if (is.na(yaml_end)) 1 else yaml_end + 1
  sentinel_line <-which(qmd_lines == "<!-- FLOAT_CHUNKS_START -->")
  body_end      <-if (length(sentinel_line) > 0) sentinel_line - 1 else length(qmd_lines)
  qmd_body      <-qmd_lines[body_start:body_end]

  #Filter file_paths lookup_table according to mode
  lookup_selected_paths<-filter_lookup_table(lookup_table = lookup_table,mode = mode)

  #Generate float chunks
  ordered_lookup<-extract_float_order(manuscript_path = manuscript_path,
                                      lookup_table    = lookup_selected_paths)

  #Add submission labels
  labeled_lookup<-ordered_lookup |>
    dplyr::group_by(key_type) |>
    dplyr::mutate(label_n = dplyr::row_number()) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      selected_file = purrr::map_chr(file_path, \(x) x[[1]]),
      ext   = tools::file_ext(selected_file),
      label = dplyr::case_when(
        key_type == "mfig" ~ paste0("Figure",               label_n, ".", selected_ext),
        key_type == "mtbl" ~ paste0("Table",                label_n, ".", selected_ext),
        key_type == "sfig" ~ paste0(sup_naming,"_Figure_S", label_n, ".", selected_ext),
        key_type == "stbl" ~ paste0(sup_naming, "_Table_S", label_n, ".", selected_ext)
      ))

  #Group by type for chunk output order
  type_order<-c("mfig","mtbl","sfig","stbl")
  labeled_lookup_chunks<-labeled_lookup |>
    dplyr::mutate(type_order = match(key_type,type_order)) |>
    dplyr::arrange(type_order) |> dplyr::select(-type_order)

  #Define float chunks and where it should begin
  float_chunks <-make_float_chunks(lookup_table = labeled_lookup_chunks)
  sentinel     <-"<!-- FLOAT_CHUNKS_START -->" #Prevents accumulating of floats between render modes

  #Create float chunk header
  float_header<-if (mode=="review") {
    paste0("**",figure_table_header,"**")
  } else {
    paste0("**",figure_table_header," LEGENDS**")
  }

  #Write updated qmd file
  writeLines(
    c(yaml_header,qmd_body,sentinel,"",float_header,"",float_chunks,"",ref_header,""),manuscript_path
    )

  #Render manuscript from .qmd file
  render_manuscript(
    mode            = mode,
    manuscript_name = manuscript_name,
    manuscript_path = manuscript_path
    )

  #Prepare submission folder
  if (mode=="submission") {
    prepare_submission(
      manuscript_path = manuscript_path,
      lookup_table    = labeled_lookup_chunks,
      output_path     = output_path,
      sup_naming      = sup_naming
    )

    #Copy rendered manuscript to submission folder
    rendered_docx<-file.path(
      dirname(manuscript_path),
      sprintf("%s_%s_%s.docx",format(Sys.Date(),"%Y.%m.%d"),mode,manuscript_name))

    file.copy(rendered_docx,file.path(
      output_path,
      sprintf("%s_%s_%s.docx",format(Sys.Date(),"%Y.%m.%d"),mode,manuscript_name)),
      overwrite=TRUE)
  }
  invisible(manuscript_path)
}



