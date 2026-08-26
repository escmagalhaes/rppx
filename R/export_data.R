### Functions to export output ###

#' Function to Export tables to excel
#' @export
#' @examples
#' ## Placeholder Example ##
export_tabs<-function(df_list,filename,filepath,tab_title) {

  if (is.data.frame(df_list)) df_list<-list(df_list)
  n<-length(df_list)
  let_id<-if (n==1) "" else LETTERS[seq_len(n)]
  sh_list<-names(df_list)
  if (is.null(sh_list) || all(sh_list=="")) sh_list<-rep("",n)

  safe_sheet<-function(x, i) {
    x<-gsub("[\\[\\]\\*\\?/\\\\:]","",x)
    x<-substr(x,1,31)
    if (x=="" || is.na(x)) x<-paste0("Sheet",i)
    x
  }
  wb<-openxlsx2::wb_workbook()

  for ( i in seq_along(df_list)) {
    df<-df_list[[ i ]]
    sheet_name<-safe_sheet(sh_list[ i ],i)
    wb$add_worksheet(sheet_name)

    tab_title_i<- if (length(tab_title)==1) tab_title else tab_title[ i ]
    title<-paste0(filename,let_id[ i ], if (let_id[ i ] !="") ". " else ". ",tab_title_i)

    #Write title row
    wb$add_data(sheet=sheet_name,x=title,start_row=1,start_col=1,col_names=FALSE)
    wb$add_font(sheet=sheet_name,dims="A1",bold=TRUE)

    #Write data
    wb$add_data(sheet=sheet_name,x=df,start_row=2,col_names=TRUE)

    #Header style: bold + double bottom border
    header_dims<-openxlsx2::wb_dims(rows=2,cols=seq_len(ncol(df)))
    wb$add_font(sheet=sheet_name,dims=header_dims,bold=TRUE)
    wb$add_border(sheet=sheet_name,dims=header_dims
                  ,bottom_color=openxlsx2::wb_color("black"),bottom_border="double")

    #Numeric column formatting
    num_cols<-which(sapply(df,is.numeric) & seq_along(df) !=1)
    data_rows<-seq(3,nrow(df)+2)

    if (length(num_cols) > 0) {
      for (col in num_cols) {
        col_vals<-df[[ col ]]
        has_fraction<-any(col_vals > 0 & col_vals < 1,na.rm=TRUE)

        if (has_fraction) {
          # Batch by value type: whole, small decimal (sci notation), regular decimal
          whole_rows<-which(col_vals == floor(col_vals) & !is.na(col_vals))+2
          small_rows<-which(col_vals < 0.0001 & col_vals !=floor(col_vals) & !is.na(col_vals))+2
          regular_rows<-which(col_vals >= 0.0001 & col_vals !=floor(col_vals) & !is.na(col_vals))+2

          if (length(whole_rows)   > 0)
            wb$add_numfmt(sheet=sheet_name,dims=openxlsx2::wb_dims(rows=whole_rows,cols=col),numfmt="0")
          if (length(small_rows)   > 0)
            wb$add_numfmt(sheet=sheet_name,dims=openxlsx2::wb_dims(rows=small_rows,cols=col),numfmt="0.00E+00")
          if (length(regular_rows) > 0)
            wb$add_numfmt(sheet=sheet_name,dims=openxlsx2::wb_dims(rows=regular_rows,cols=col),numfmt="0.######")

        } else {
          #Count or ratio column: single format for whole column
          is_count<-all(col_vals==floor(col_vals),na.rm=TRUE)
          wb$add_numfmt(sheet=sheet_name
                        ,dims=openxlsx2::wb_dims(rows=data_rows,cols=col)
                        ,numfmt=if (is_count) "0" else "0.######")
        }
      }
    }
    #Center internal columns (cols 2 to ncol-2)
    mid_cols<- if (ncol(df) > 3) seq(2,ncol(df)-2) else NULL
    if (!is.null(mid_cols)) {
      #Data rows: center alignment
      wb$add_cell_style(sheet=sheet_name
                        ,dims=openxlsx2::wb_dims(rows=data_rows,cols=mid_cols)
                        ,horizontal="center")
      #Header row: bold + center alignment
      wb$add_font(sheet=sheet_name
                  ,dims=openxlsx2::wb_dims(rows=2,cols=mid_cols)
                  ,bold=TRUE)
      wb$add_cell_style(sheet=sheet_name
                        ,dims=openxlsx2::wb_dims(rows=2,cols=mid_cols)
                        ,horizontal="center")
    }
    #Column widths
    wb$set_col_widths(sheet=sheet_name,cols=1,widths=25)
    wb$set_col_widths(sheet=sheet_name,cols=2:ncol(df),widths="auto")
  }
  wb$save(file.path(filepath,paste0(filename,".xlsx")))
}

#' Function to export tables to excel spreadsheets
#' @export
#' @examples
#' ## Placeholder Example ##
export_bxplot_tests<-function(results,filename,filepath,tab_title) {

  test_tables <- lapply(names(results), function(tab_names) {
    df <- results[[tab_names]]$test_table
    if (is.null(df)) return(NULL)

    #Remove plotting columns
    df <- df[, !colnames(df) %in% c("y.position", "groups", "xmin", "xmax", "step.increase"), drop = FALSE]

    #Sort by adjusted p-value if present
    if ("p.adj" %in% colnames(df)) {
      df <- df[order(df$p.adj), ]
    }

    df
  })

  names(test_tables) <- names(results)

  #Remove empty entries
  test_tables <- test_tables[!sapply(test_tables, is.null)]

  #Export
  export_tabs(df_list=test_tables,filename=filename,filepath=filepath,tab_title=tab_title)
}

#' Function to export boxplots
#' @export
#' @examples
#' ## Placeholder Example ##
export_bxplots_pdf<-function(results,file=filename,ncol=3,nrow=2){
  plots<-lapply(results,function(x) x$bxplot)
  per_page<-ncol*nrow
  pdf(file,width=11.69,height=8.27)

  for(i in seq(1,length(plots),by=per_page)){
    page<-plots[ i :min( i +per_page-1,length(plots))]
    print(ggpubr::ggarrange(plotlist=page,ncol=ncol,nrow=nrow))
  }
  dev.off()
}

