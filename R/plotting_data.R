### Functions for general plotting ###

#' Function to make density plots
#' @export
#' @examples
#' ## Placeholder Example ##
plot_density<-function(df,protein,grouping_var=NULL,line_var_label="Type",
                       df_sample_type="AML",title_suffix="AML",
                       df_normal=NULL,normal_sample_type="CD34+",
                       covar=NULL,covar_label=NULL,facet=FALSE,n_panels=1,
                       fill_colors=NULL,color_fill_mode=c('list','vector'),
                       line_colors=RColorBrewer::brewer.pal(7,"Dark2"),
                       normal_fill_colors="grey80",default_density_color="steelblue",
                       legend_mode=c("inside","none","right","left","bottom","top"),
                       legend_position=c(0.75,0.75)) {

  color_fill_mode<-match.arg(color_fill_mode)

  #Define plot scale accoding to panels using helper function
  get_scale<-function(n_panels) {
    return(min(1, 1 / sqrt(n_panels) * 1.5))
  }
  scale<-get_scale(n_panels)

  #Safety checks
  if(!protein %in% names(df))
    stop("\nProtein not found in df\n")

  if(!is.null(df_normal) && !protein %in% names(df_normal))
    stop("\nProtein not found in df_normal\n")

  if(!is.null(grouping_var) && !grouping_var %in% names(df))
    stop(paste("\nVariable",grouping_var,"not found in df\n"))

  if(!is.null(covar) && !covar %in% names(df))
    stop(paste("\nVariable",covar,"not found in df\n"))

  if(!is.null(covar) && !is.null(grouping_var))
    stop("\n If covar is passed then grouping_var has to be NULL\n")

  if (!is.null(grouping_var)) {

    #Add labels to data
    df<-df |> group_by(.data[[grouping_var]]) |>
      mutate(n   = dplyr::n(),
             pct = n / nrow(df) * 100,
             label = paste0(.data[[grouping_var]],"; N=",n," (",round(pct,1),"%)")) |> ungroup()

    #Make sure labels variable have the same order as grouping_variable
    label_levels<-df |> distinct(.data[[grouping_var]],label) |>
      arrange(match(.data[[grouping_var]],levels(df[[grouping_var]]))) |> pull(label)
    df$label<-factor(df$label,levels=label_levels)

    #Get density values
    dens<-density(df[[protein]])
    dens_df<-data.frame(x=dens$x,y=dens$y)

    #Compute data cutoffs according to grouping_var to apply to density values
    cutoff_vals<-tapply(df[[protein]],df[[grouping_var]],max)
    cutoff_vals<-cutoff_vals[levels(df[[grouping_var]])]
    n_levels<-length(cutoff_vals)

    #Determine the grouping_var group of each dens_df value according to cutoffs
    idx<-rowSums(outer(dens_df$x,cutoff_vals,`>`))+1
    idx<-pmin(idx,n_levels)
    dens_df$region<-factor(label_levels[idx],levels=label_levels)

    #Create variable to identify the sample_type if non-existing
    dens_df$sample_type<-factor(rep(df_sample_type,nrow(dens_df)))

    #Compute ranges and median
    range_ptn_df<-paste0(round(min(df[[protein]],na.rm=TRUE),2),
                         " to ",round(max(df[[protein]],na.rm=TRUE),2))
    med_ptn_df<-median(df[[protein]],na.rm=TRUE)

    #Get fill colors
    if (color_fill_mode=="list") {
      if (is.null(fill_colors)) fill_colors<-get_mycolors(mode="list")
      if (n_levels > length(fill_colors)) {
        warning(paste0("Not enough colors in for ",n_levels,
                       " levels. Interpolating palette."))
        fill_colors<-colorRampPalette(fill_colors)(n_levels)
      }
      fill_vals<-get_mycolors(n_colors=n_levels,mode="list",color_list=fill_colors)
    } else {
      if (is.null(fill_colors)) fill_colors<-get_mycolors(mode="vector")
      if (n_levels > length(fill_colors)) {
        warning(paste0("Not enough colors in for ",n_levels,
                       " levels. Interpolating palette."))
        fill_colors<-colorRampPalette(fill_colors)(n_levels)
      }
      fill_vals<-get_mycolors(n_colors=n_levels,mode="vector",color_vector=fill_colors)
    }
    fill_vals_named<-setNames(fill_vals,label_levels)

    #Adjust line colors according to level
    line_colors_named<-setNames(line_colors[1],df_sample_type)

    #Density plot
    dplot<-ggplot()+
      geom_ribbon(data=dens_df,aes(x=x,ymin=0,ymax=y,fill=region),alpha=1)+
      geom_line(data=dens_df,aes(x=x,y=y,color=sample_type)
                ,linetype="solid",linewidth=1 * scale)+
      geom_vline(xintercept=med_ptn_df,color=line_colors_named[df_sample_type]
                 ,linetype="dashed",linewidth=1 * scale)+
      labs(title=paste(protein,"levels in",title_suffix),x=paste(protein,"levels"),y=NULL,
           subtitle=paste0(df_sample_type,": ","Median = ",round(med_ptn_df,2)
                           ," | Range = ",range_ptn_df,"\n")
      )

    #Add normal density plot if df_normal is passed
    if (!is.null(df_normal)) {

      #Calculate density for normal data
      dens_normal<-density(df_normal[[protein]])
      dens_normal_df<-data.frame(x=dens_normal$x,y=dens_normal$y)

      #Create variable to identify the sample_type if non-existing
      dens_normal_df$sample_type<-factor(rep(normal_sample_type,nrow(dens_normal_df)))

      #Calculate range and median of normal data
      range_normal<-paste0(round(min(df_normal[[protein]],na.rm=TRUE),2),
                           " to ",round(max(df_normal[[protein]],na.rm=TRUE),2))
      med_normal<-median(df_normal[[protein]],na.rm=TRUE)

      #Adjust line colors according to level
      line_colors_named<-setNames(line_colors[1:2],c(df_sample_type,normal_sample_type))

      #Adjust fill colors
      fill_vals_named<-c(fill_vals_named,setNames(normal_fill_colors,normal_sample_type))

      #Overlay normal on top of the initial plot
      dplot<-dplot+
        geom_ribbon(data=dens_normal_df,aes(x=x,ymin=0,ymax=y),fill=normal_fill_colors,
                    alpha=0.25,inherit.aes=FALSE,show.legend=FALSE)+
        geom_line(data=dens_normal_df,aes(x=x,y=y,color=sample_type)
                  ,linetype="solid",alpha=0.7,linewidth=1 * scale)+
        geom_vline(xintercept=med_normal,color=line_colors_named[normal_sample_type]
                   ,linetype="dashed",alpha=0.7,linewidth=1 * scale)+
        labs(title=paste(protein,"levels in",df_sample_type,"vs",normal_sample_type),
             subtitle=paste0(normal_sample_type,": ","Median = ",round(med_normal,2)
                             ," | Range = ",range_normal,"\n",
                             df_sample_type,": ","Median = ",round(med_ptn_df,2)
                             ," | Range = ",range_ptn_df,"\n"))
    }

    #Apply colors conditionally
    if (!is.null(df_normal)) {
      dplot<-dplot +
        scale_fill_manual(name=paste(protein,"levels"),values=fill_vals_named) +
        scale_color_manual(name=line_var_label,
                           values=setNames(
                             c(line_colors_named[df_sample_type],
                               line_colors_named[normal_sample_type]),
                             c(df_sample_type, normal_sample_type) ))
    } else {
      dplot<-dplot +
        scale_fill_manual(name=paste(protein,"levels"),values=fill_vals_named) +
        scale_color_manual(name=line_var_label,
                           values=line_colors_named[df_sample_type])
    }
  }

  if (!is.null(covar)) {

    #Filter to complete cases only in covar
    df_filt<-df |> filter(!is.na(.data[[covar]]))
    covar_levels<-levels(df_filt[[covar]])
    n_covar_levels<-nlevels(df_filt[[covar]])

    #Get fill colors
    if (color_fill_mode=="list") {
      if (is.null(fill_colors)) fill_colors<-get_mycolors(mode="list")
      if (n_covar_levels > length(fill_colors)) {
        warning(paste0("Not enough colors in for ",n_covar_levels,
                       " levels. Interpolating palette."))
        fill_colors<-colorRampPalette(fill_colors)(n_covar_levels)
      }
      fill_vals<-get_mycolors(n_colors=n_covar_levels,mode="list",color_list=fill_colors)
    } else {
      if (is.null(fill_colors)) fill_colors<-get_mycolors(mode="vector")
      if (n_covar_levels > length(fill_colors)) {
        warning(paste0("Not enough colors in for ",n_covar_levels,
                       " levels. Interpolating palette."))
        fill_colors<-colorRampPalette(fill_colors)(n_covar_levels)
      }
      fill_vals<-get_mycolors(n_colors=n_covar_levels,mode="vector",color_vector=fill_colors)
    }
    fill_vals_named<-setNames(fill_vals,covar_levels)

    #Adjust line colors according to level
    line_colors_named<-fill_vals_named

    #Calculate range and median of normal data
    med_covar<-df_filt |> group_by(.data[[covar]]) |>
      summarise(median=median(.data[[protein]],na.rm=TRUE))

    #Density plot
    dplot<-ggplot(df_filt,aes(x=.data[[protein]],color=.data[[covar]],fill=.data[[covar]]))+
      geom_density(alpha=0.3,linewidth=1 * scale)+
      geom_vline(data=med_covar,aes(xintercept=median,color=!!sym(covar))
                 ,linetype="dashed",linewidth=1 * scale)+
      scale_fill_manual(name=covar_label,values=fill_vals_named)+
      scale_color_manual(name=covar_label,values=line_colors_named)+
      labs(title=paste(protein,"levels in",title_suffix),x=paste(protein,"levels"),y=NULL)

    #If split by covar into separate plots
    if (facet) {
      dplot<-dplot+facet_wrap(as.formula(paste("~", covar)))
    }
  }

  if (is.null(covar) && is.null(grouping_var) ) {

    #Density plot simple
    dplot<-ggplot(df,aes(x=.data[[protein]]))+
      geom_density(alpha=0.4,linewidth=1 * scale
                   ,color=default_density_color
                   ,fill=default_density_color)+
      labs(x=paste(protein,"levels"),y=NULL)+
      labs(title=paste(protein,"levels in",title_suffix),x=NULL,y=NULL)
  }

  #Compute legend theme dinamically
  legend_dplot_theme<-if (isTRUE(facet)) {
    theme(legend.position = "none")
  } else if (legend_mode=="inside") {
    theme(legend.position        = "inside",
          legend.position.inside = legend_position)
  } else {
    theme(legend.position = legend_mode)
  }

  #Compute theme design
  dplot_theme<-
    theme_classic(base_size= 18 * scale
                  ,base_line_size=max(0.2, 0.3 * scale)
                  ,base_rect_size=max(0.2, 0.3 * scale))+
    theme(
      plot.title             = element_text(size = rel(1.2) ),
      axis.title             = element_text(size = rel(1.1) ),
      plot.subtitle          = element_text(size = rel(0.8), face = "italic"),
      legend.title           = element_text(size = rel(0.8)),
      legend.text            = element_text(size = rel(0.7)),
      legend.key.spacing.y   = unit(0.1 * scale, "cm"),
      legend.key.size        = unit(0.6 * scale, "cm"),
      legend.background      = element_rect(fill=NA,color=NA),
      legend.box.background  = element_rect(fill=NA,color=NA),
      legend.key             = element_rect(fill=NA,color=NA),
      plot.margin            = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
    )

 #Apply themes to density plot
 dplot<-dplot+dplot_theme+legend_dplot_theme

  return(dplot)
}

#' Function to make KM plots
#' @export
#' @examples
#' ## Placeholder Example ##
plot_km<-function(data,time_var,event_var,auto_style=FALSE,grouping_var=NULL,covar=NULL
                  ,title=NULL,ylab='Cumulative probability',xlab="Time (years)"
                  ,xlim=NULL,break_x_by=NULL,table_by_time=FALSE,table_times=NULL
                  ,legend_mode=c("inside","none","right","left","bottom","top")
                  ,legend_position=c(0.7,0.7),legend_labs=NULL,ncol_legend=1
                  ,pval_coord=c(3.5,1),risk_table=FALSE,single_color=NULL
                  ,show_pval_table=TRUE,pval_table_adj_method='BH',n_panels=1
                  ,add_overall=TRUE,linetype_overall="solid",color_overall="black"
                  ,line_types=c("solid","twodash","dashed","dotted","dotdash","longdash")
                  ,colors=c('red2','blue3','green4','purple2','darkorange'
                            ,'darkcyan','darkgoldenrod3','yellow2','deeppink2'
                            ,'lightsalmon2','darkturquoise')) {

  legend_mode<-match.arg(legend_mode)
  missing_levels<-character(0)   #initialise here for interaction between grouping_var and covar

  #Safety checks
  if (!is.null(grouping_var) && !is.factor(data[[grouping_var]]))
    stop("grouping_var must be a factor.")

  if (!is.null(covar) && !is.factor(data[[covar]]) )
    stop("covar must be a factor if provided.")

  #Generate model formula according to variables
  if (is.null(grouping_var) && is.null(covar)) {
    #KM curve with no grouping
    surv_formula<-as.formula(paste0("Surv(",time_var,",",event_var,")~1"))
    int_levels<-NULL
    n_colors<-1L
    n_lt<-1L

  } else if (!is.null(grouping_var) && is.null(covar)) {
    #Single grouping variable
    surv_formula<-as.formula(paste0("Surv(",time_var,",",event_var,")~",grouping_var))
    int_levels<-levels(data[[grouping_var]])
    n_colors<-length(int_levels)
    n_lt<-1L

  } else if (!is.null(grouping_var) && !is.null(covar)) {
    #Drop any stale levels before building the grid
    data[[grouping_var]]<-droplevels(data[[grouping_var]])
    data[[covar]]<-droplevels(data[[covar]])

    #Generate Interaction variable if covar exists
    data$int_var<-interaction(data[[grouping_var]],data[[covar]]
                              ,lex.order=TRUE,drop=TRUE,sep='\n')
    surv_formula<-as.formula(paste0("Surv(",time_var,",",event_var,")~int_var"))

    #All possible combinations vs those that actually exist
    all_combinations<-interaction(data[[grouping_var]],data[[covar]],
                                  lex.order=TRUE,drop=FALSE,sep='\n')
    all_levels<-levels(all_combinations)
    exist_levels<-levels(data$int_var)

    #Warning about missing combinations
    missing_levels<-setdiff(all_levels,exist_levels)
    if (length(missing_levels) > 0) {
      warning("The following combinations have no cases and will be skipped:\n  ",
              paste(missing_levels,collapse="\n  "))
    }
    int_levels<-exist_levels
    n_colors<-nlevels(data[[grouping_var]])
    n_lt<-nlevels(data[[covar]])

  } else if (is.null(grouping_var) && !is.null(covar)) {
    #Covar only: single color, linetype varies
    surv_formula<-as.formula(paste0("Surv(", time_var,",",event_var,")~",covar))
    int_levels<-levels(data[[covar]])
    n_colors<-1L
    n_lt<-nlevels(data[[covar]])

    #Stop if single_color argument is not provided
    if (is.null(single_color))
      stop("single_color must be provided when only covar is specified.")

    #Warning/guarding about single color in case grouping_var is not provided
    if (length(single_color) > 1) {
      warning("single_color should be a single color string. Using only the first value.")
      single_color<-single_color[1]
    }
    colors_used<-single_color
  }

  #Auto-compute xlim & break_x_by
  if (is.null(xlim)) {
    x_upper<-ceiling(max(data[[time_var]],na.rm=TRUE))
    if (x_upper %% 2==0) x_upper<-x_upper+1
    xlim<-c(0,x_upper)
  } else {
    x_upper<-xlim[2]
  }
  if (is.null(break_x_by)) {
    break_x_by<-pretty(c(0,x_upper),n=10) |> diff() |> unique()
  }

  #Colors
  if (is.null(single_color)) {
    if (n_colors > length(colors)) {
      warning("grouping_var has ",n_colors," levels but colors only has ",length(colors),
              " entries. Recycling colors...")
    }
    colors_used<-rep(colors,length.out=n_colors) #recycles if n_colors > length(colors)
  }

  #Linetypes
  if (!is.null(covar) && n_lt > length(line_types)) {
    warning("covar has ",n_lt," levels but line_types only has ",length(line_types),
            " entries. Recycling line_types...")
  }
  line_types_used<-rep(line_types,length.out=n_lt)  #recycles if n_lt > length(line_types)

  #Generate color/linetype combo vectors
  if (!is.null(covar) && !is.null(grouping_var)) {
    #Explicit named lookup table: robust to missing combos and factor reordering
    color_levels_all<-levels(data[[grouping_var]])
    lt_levels_all<-levels(data[[covar]])

    full_grid<-expand.grid( gp = color_levels_all, lt = lt_levels_all
                            , stringsAsFactors=FALSE) |>
      dplyr::mutate(
        int_level    = paste(gp, lt , sep = '\n'),
        color_mapped = colors_used[match( gp , color_levels_all)],
        lt_mapped    = line_types_used[match( lt , lt_levels_all)]
      )

    #Filter to existing levels, preserving original identity and level order
    full_grid_exist<-full_grid |>
      dplyr::filter(int_level %in% exist_levels) |>
      dplyr::arrange(match(int_level,exist_levels))

    color_seq<-full_grid_exist$color_mapped
    full_linetype_vec<-full_grid_exist$lt_mapped

  } else if (!is.null(covar) && is.null(grouping_var)) {
    #Covar present, grouping_var absent: single color, linetypes vary
    color_seq<-rep(colors_used, n_lt)
    full_linetype_vec<-line_types_used[seq_len(n_lt)]

  } else {
    #Covar, but grouping_var present: simple color sequence, all same linetype
    color_seq<-colors_used
    full_linetype_vec<-rep(line_types_used[1], n_colors)
  }

  #Add overall parameters if add_overall=TRUE
  pal_vec<-if (add_overall) c(color_overall,color_seq)  else color_seq
  linetype_vec<-if (add_overall) c(linetype_overall,full_linetype_vec) else full_linetype_vec

  #Legend labels
  if(is.null(legend_labs)){
    legend_labs<-if (add_overall) c("Overall",int_levels) else int_levels
  }

  #Fit surv model & plot
  mod_fit<-surv_fit(surv_formula,data=data)

  #Determine plot and pval_table scale based on panel number
  get_scale<-function(n_panels) {
    return(min(1, 1 / sqrt(n_panels) * 1.5))
  }
  scale<-get_scale(n_panels)

  #KM plot
  km_plot<-ggsurvplot(mod_fit,data=data,
                      legend.title        = " ",
                      legend.labs         = if (!auto_style) legend_labs  else NULL,
                      pal                 = if (!auto_style) pal_vec      else NULL,
                      linetype            = if (!auto_style) linetype_vec else 1,
                      xlim                = xlim,
                      break.x.by          = break_x_by,
                      conf.int            = FALSE,
                      title               = title,
                      xlab                = xlab,
                      ylab                = ylab,
                      risk.table.height   = 0.3,
                      tables.col          = "strata",
                      risk.table          = risk_table,
                      risk.table.fontsize = NULL,
                      size                = 1 * scale,
                      fontsize            = 10 * scale,
                      add.all             = add_overall,
                      pval                = TRUE,
                      pval.coord          = pval_coord,
                      pval.size           = 6 * scale,
                      censor.size         = 5 * scale,
                      ggtheme=theme_classic(base_size= 13 * scale
                                            ,base_line_size=max(0.2, 0.3 * scale)
                                            ,base_rect_size=max(0.2, 0.3 * scale))
  )

  #Compute legend theme dinamically
  legend_theme<-if (legend_mode=="inside") {
    theme(legend.position        = legend_mode,
          legend.position.inside = legend_position)
  } else {
    theme(legend.position = legend_mode )
  }

  #Adjust plot theme
  km_plot$plot<-km_plot$plot+theme(
    plot.title             = element_text(size = rel(1.2) ),
    axis.title             = element_text(size = rel(1.1) ),
    axis.text              = element_text(size = rel(0.8) ),
    legend.text            = element_text(size = rel(0.8) ),
    legend.key.size        = unit( 2.5 * scale ,"line"),
    legend.background      = element_rect(fill=NA,color=NA),
    legend.box.background  = element_rect(fill=NA,color=NA),
    legend.key             = element_rect(fill=NA,color=NA),
    plot.margin            = unit(c(0.1, 0.1, 0.1, 0.1), "cm")
  )+legend_theme+guides(colour=guide_legend(ncol=ncol_legend))

  #Adjust risk table theme
  km_plot$table<-km_plot$table +
    theme(
      legend.position        = "none",
      legend.text            = element_text( size = rel(0.8) ),
      legend.key.size        = unit( 2.5 * scale , "line"),
      legend.background      = element_rect(fill=NA,color=NA),
      legend.box.background  = element_rect(fill=NA,color=NA),
      legend.key             = element_rect(fill=NA,color=NA),
      axis.title.x           = element_text(size = rel(0.8) ),
      axis.text              = element_text(size = rel(0.8) ),
      plot.margin            = unit(c( 0.1, 0.1 ,0.1 ,0.1 ), "cm")
    )

  #Compute pairwise tests
  if (isTRUE(show_pval_table)) {
    pairwise_syms<-symnum(
      pairwise_survdiff(surv_formula,data=data,p.adjust.method=pval_table_adj_method)$p.value,
      cutpoints     = c( 0, 0.001, 0.01, 0.05 ,1 ),
      symbols       = c( "p<0.001", "p<0.01", "p<0.05", "ns" ),
      abbr.colnames = FALSE,
      na            = "-"
    )

    #Pval table fill vectors
    #col_fills must have length n_int, row_fills must have length n_int - 1
    #color_seq has length n_colors * n_lt == n_int (after filtering)
    n_int<-length(int_levels)
    col_fills<-color_seq
    row_fills<-c("white",color_seq[-1])

    #Helper function to adjust pval_tab widths
    col_widths<-pmax(apply(pairwise_syms,2,function(col) {
        max(nchar(as.character(col)),na.rm=TRUE)
      }), nchar(colnames(pairwise_syms)))

    #Pairwise p-value table
    pval_tab<-tableGrob(
      pairwise_syms,
      theme=ttheme_minimal(
        core    = list(fg_params = list(cex  = 1.75 * scale, fontface = 2 ),
                       bg_params = list(fill="white",col="black"),
                       padding   = unit(c( 10 * scale,  3 * scale), "mm")
                       ),
        colhead = list(fg_params = list(cex  = 1.25 * scale, col ="white",fontface = 2 ),
                       bg_params = list(fill = col_fills,col=NA),
                       padding   = unit(c( 10 * scale,  4 * scale), "mm")
                       ),
        rowhead = list(fg_params = list(cex  = 1.25 * scale, col  = "white",fontface = 2 ),
                       bg_params = list(fill = row_fills, col=NA),
                       padding   = unit(c( 2 * scale,  3 * scale), "mm")
                       )
      )#,widths  = unit(col_widths, "null")
      )
  } else {
    pval_tab<-NULL
  }

  #Survival table, if computed
  if (table_by_time) {

    #Compute times to be used
    if (is.null(table_times)) {
      table_times<-seq(1,x_upper,by=break_x_by)
    }

    #Median + CI per strata
    med_fit<-as.data.frame(summary(mod_fit)$table) |>
      tibble::rownames_to_column("strata") |>
      dplyr::select(strata,n=records,events,median,ci_lower_med=`0.95LCL`,ci_upper_med=`0.95UCL`)

    #Survival probability at each time point
    fit_by_time<-summary(mod_fit,times=table_times,extend=TRUE)

    #Merge tables
    fit_pivot_wide<-tibble::tibble(
      strata=as.character(fit_by_time$strata)
      ,time=fit_by_time$time
      ,surv=round(fit_by_time$surv,2)) |>
      tidyr::pivot_wider(names_from=time,values_from=surv,names_prefix="t=")
    fit_tab<-dplyr::left_join(med_fit,fit_pivot_wide,by="strata")

  } else {
    fit_tab<-NULL
  }

  #Assemble & return
  invisible(list(
    plot         = km_plot$plot,
    risk_table   = km_plot$table,
    pval_table   = pval_tab,
    time_table   = fit_tab,
    pal_vec      = pal_vec,        # for debugging
    linetype_vec = linetype_vec,   # for debugging
    legend_labs  = legend_labs     # for debugging
  ))
}

#' Function to make upset plots
#' @export
#' @examples
#' ## Placeholder Example ##
plot_upset<-function(df_list,groups=NULL,colors=NULL,
                     only_color="grey70",exclusive_color="gold",
                     highlight_exclusive=FALSE,plot_title=NULL){

  if(length(df_list)==0) stop("df_list is empty")

  upset_mod<-make_comb_mat(df_list)
  comb_names<-comb_name(upset_mod)
  n_comb<-length(comb_names)
  set_names<-set_name(upset_mod)

  #Sanity check
  if(!is.null(groups)){
    missing_groups<-setdiff(groups,set_names)
    if(length(missing_groups) > 0)
      warning("These groups are not in the data: ",paste(missing_groups,collapse=", "))
  }

  #convert combinations to matrix
  comb_mat<-do.call(rbind,strsplit(comb_names,""))
  comb_mat<-apply(comb_mat,2,as.numeric)

  #default color
  comb_colors<-rep("black",n_comb)

  if(!is.null(groups)){
    if(is.null(colors))
      colors<-rep("red",length(groups))
    if(length(colors) != length(groups))
      stop("Length of 'colors' must match length of 'groups'")
    #map group indices
    group_idx<-match(groups,set_names)
    #number of sets per combination
    set_count<-rowSums(comb_mat)
    for(i in seq_along(groups)){
      idx<-group_idx[i]
      if(is.na(idx)) next
      in_group<-comb_mat[,idx]==1
      #priority coloring
      comb_colors[in_group]<-colors[i]
      #group-only intersections
      only_group<-in_group & set_count==1
      comb_colors[only_group]<-only_color
      #optional exclusive intersections
      if(highlight_exclusive){
        exclusive<-in_group & set_count > 1
        comb_colors[exclusive]<-exclusive_color
      }
    }
  }
  upset_plot<-as.ggplot(
    UpSet(upset_mod,pt_size=unit(4,"mm"),lwd=3,comb_col=comb_colors,
          column_title=plot_title,
          column_title_gp=gpar(fontsize=10,fontface="bold"),
          #,left_annotation=upset_left_annotation(upset_mod,add_numbers=TRUE)
          top_annotation=upset_top_annotation(upset_mod,add_numbers=TRUE),
          comb_order=order(-comb_size(upset_mod))))
  return(list(plot=upset_plot,model=upset_mod))
}

#' Function to extract protein names from upset plot
#' @export
#' @examples
#' ## Placeholder Example ##
get_upset_ptns<-function(comb_mat) {

  comb_names<-ComplexHeatmap::comb_name(comb_mat)
  set_names<-ComplexHeatmap::set_name(comb_mat)

  #Helper: convert "110" into labels
  comb_to_label<-function(comb) {
    bits<-as.numeric(strsplit(comb,"")[[1]])
    paste(set_names[bits==1],collapse=" & ")
  }

  #Extract proteins per combination
  comb_list<-lapply(comb_names, function(x) {
    ComplexHeatmap::extract_comb(comb_mat, x)
  })

  #Assign readable names
  names(comb_list)<-sapply(comb_names,comb_to_label)

  return(comb_list)
}

#' Function to make heatmaps
#' @export
#' @examples
#' ## Placeholder Example ##
plot_heatmap<-function(df,proteins,
                       col_anno=TRUE,grouping_var=NULL,covars=NULL,covar_labels=NULL,
                       grouping_var_label=NULL,row_anno=FALSE,row_anno_df=NULL,
                       row_anno_protein_col=NULL,rown_anno_vars=NULL,rown_anno_vars_labels=NULL,
                       colors_ht=NULL,covar_colors_override=NULL,
                       colors_anno_col=c(pal_jco()(2)[2:1],'red2',brewer.pal(7,"Dark2")[c(7,1,2:6)]),
                       colors_anno_row=c(pal_jco()(2)[2:1],'red2',brewer.pal(7,"Dark2")[c(7,1,2:6)]),
                       colors_grouping_var=brewer.pal(7,"Dark2"),colors_NA=pal_jco()(3)[3],
                       plot_title=NULL,cluster_rows=TRUE,cluster_cols=FALSE,
                       legend_breaks=NULL,legend_labels=NULL,breaks=NULL,protein_labels=NULL,
                       show_colnames=FALSE,show_rownames=TRUE,annotation_names_row=FALSE){

  #Input validation
  if (isTRUE(col_anno)) {
    if (is.null(grouping_var) && is.null(covars))
      stop("col_anno=TRUE but both grouping_var and covars are NULL.")
    if (!is.null(grouping_var) && !grouping_var %in% names(df))
      stop(sprintf("grouping_var '%s' not found in input df.",grouping_var))
    if (!is.null(covars) && !all(covars %in% names(df)))
      stop("Some covars not found in input df.")
  }

  if (isTRUE(row_anno)) {
    if (is.null(row_anno_df))
      stop("row_anno=TRUE but row_anno_df is NULL.")
    if (is.null(row_anno_protein_col))
      stop("row_anno=TRUE but row_anno_protein_col is NULL. Specify which column holds protein names.")
    if (!row_anno_protein_col %in% names(row_anno_df))
      stop(sprintf("row_anno_protein_col '%s' not found in row_anno_df.", row_anno_protein_col))
    if (is.null(rown_anno_vars))
      stop("row_anno=TRUE but rown_anno_vars is NULL.")
    if (!all(rown_anno_vars %in% names(row_anno_df)))
      stop("Some rown_anno_vars not found in row_anno_df.")
  }

  #Dataframe subsetting
  df<-as.data.frame(df)
  subset_cols<-c(proteins,grouping_var,covars)
  df<-df[,subset_cols,drop=FALSE]

  #Factor and order grouping_var
  if (!is.null(grouping_var)) {
    df[[grouping_var]]<-factor(df[[grouping_var]])
    df<-df[order(df[[grouping_var]]),]
  }

  #Helper function to set levels and colors
  fix_annotation<-function(x,colors_vec,colors_NA) {

    #Capture factor info before coercion
    is_fct<-is.factor(x)
    existing_lvls<-if (is_fct) levels(x) else NULL

    x<-as.character(x)
    x[is.na(x)]<-"N/A"

    unique_vals<-unique(x[x != "N/A"])
    is_yes_no<-all(unique_vals %in% c("no", "yes"))

    if (is_yes_no) {
      #Always force no, yes, N/A order
      lvls <- c("no", "yes")
      if (any(x=="N/A")) lvls <- c(lvls, "N/A")
    } else if (is_fct) {
      # Preserve original factor level order, append N/A if needed
      lvls <- existing_lvls
      if (any(x=="N/A") && !"N/A" %in% lvls) lvls<-c(lvls,"N/A")
    } else {
      #Use order-as-they-appear, append N/A if needed
      lvls <- unique(x[x != "N/A"])
      if (any(x=="N/A")) lvls<-c(lvls,"N/A")
    }
    x<-factor(x,levels= lvls)
    color_vec<-colors_vec[seq_along(lvls)]

    if ("N/A" %in% lvls) color_vec[lvls=="N/A"]<-colors_NA
    names(color_vec)<-lvls

    return(list(factor=x,colors=color_vec))
  }

  #Column annotation
  if (isTRUE(col_anno)) {
    anno_cols_to_use<-c(grouping_var,covars)
    col_ann<-df[,anno_cols_to_use,drop=FALSE]

    #Labels: fall back to original names if labels not provided
    col_labels<-c(
      if (!is.null(grouping_var)) if (!is.null(grouping_var_label)) grouping_var_label else grouping_var,
      if (!is.null(covars)) if (!is.null(covar_labels)) covar_labels else covars
    )
    colnames(col_ann)<-col_labels
    row.names(col_ann)<-seq_len(nrow(df))

    col_ann_colors<-list()
    for ( i in seq_along(col_ann)) {
      if (colnames(col_ann)[ i ]==col_labels[1] && !is.null(grouping_var)) {
        color_vec<-colors_grouping_var
      } else {
        color_vec<-colors_anno_col
      }
      adj_col_ann_colors<-fix_annotation(col_ann[[ i ]],color_vec,colors_NA)

      #Override colors if provided
      var_name<-colnames(col_ann)[ i ]
      if (!is.null(covar_colors_override) && var_name %in% names(covar_colors_override)) {
        adj_col_ann_colors$colors<-covar_colors_override[[var_name]]

        #Auto-fill N/A color if not specified in override
        if ("N/A" %in% levels(adj_col_ann_colors$factor) && !"N/A" %in% names(adj_col_ann_colors$colors)) {
          adj_col_ann_colors$colors["N/A"]<-colors_NA
        }

      }

      col_ann[[ i ]]<-adj_col_ann_colors$factor
      col_ann_colors[[colnames(col_ann)[ i ]]]<-adj_col_ann_colors$colors
    }
  } else {
    col_ann<- NA
    col_ann_colors<-list()
  }

  #Row annotation
  if (isTRUE(row_anno)) {
    #Subset row_anno dataframe
    row_ann<-row_anno_df[,c(row_anno_protein_col,rown_anno_vars),drop=FALSE]
    row_ann<-row_ann[order(row_ann[[row_anno_protein_col]]),]

    #Set rownames as protein column and then drop it
    rownames(row_ann)<-row_ann[[row_anno_protein_col]]
    row_ann[[row_anno_protein_col]]<-NULL

    #Apply labels
    if (!is.null(rown_anno_vars_labels)) colnames(row_ann)<-rown_anno_vars_labels

    row_ann_colors<-list()
    for ( i in seq_along(row_ann)) {
      adj_row_ann_colors<-fix_annotation(row_ann[[ i ]],colors_anno_row,colors_NA)
      row_ann[[ i ]]<-adj_row_ann_colors$factor
      row_ann_colors[[colnames(row_ann)[ i ]]]<-adj_row_ann_colors$colors
    }
  } else {
    row_ann<- NA
    row_ann_colors<-list()
  }

  #Merge annotation colors lists
  ann_colors<-c(col_ann_colors,row_ann_colors)
  if (length(ann_colors)==0) ann_colors<-NULL

  #Clean heatmap matrix
  ht_mtx<-df[,sort(proteins)]
  ht_mtx<-as.data.frame(t(ht_mtx))
  ht_mtx<-as.matrix(sapply(ht_mtx,function(x) as.numeric(as.character(x))))
  colnames(ht_mtx)<-seq_len(ncol(ht_mtx))
  rownames(ht_mtx)<-sort(proteins)

  #Make sure row_anno is aligned with matrix
  if (isTRUE(row_anno)) { row_ann<-row_ann[rownames(ht_mtx),,drop=FALSE] }

  #Apply protein labels if provided
  if (!is.null(protein_labels)) {

    #Safety check
    if (length(protein_labels) != length(proteins))
      stop("protein_labels must be the same length as proteins.")

    #Map labels onto the sorted protein order
    names(protein_labels)<-proteins
    rownames(ht_mtx)<-protein_labels[sort(proteins)]

    if (isTRUE(row_anno)) { rownames(row_ann)<-protein_labels[sort(proteins)] }
  }

  #Default colors and breaks
  if (is.null(colors_ht)) {
    #Normalize sample signal and create color object
    if(min(ht_mtx)>=-2 & max(ht_mtx)<=2){
      breaks=seq(-2.4,2.4,0.4)
    }else if(min(ht_mtx)>=-2){
      breaks=c(seq(-2.4,2,0.4),max(ht_mtx))
    }else if(max(ht_mtx)<=2){
      breaks=c(min(ht_mtx),seq(-2,2.4,0.4))
    }else{
      breaks=c(min(ht_mtx),seq(-2,2,0.4),max(ht_mtx))
    }
    colors_ht<-matlab::jet.colors(length(breaks)-1)
  }

  ht<-ggplotify::as.ggplot(pheatmap::pheatmap(ht_mtx
                                              ,annotation_colors=ann_colors
                                              ,main=plot_title
                                              ,annotation_col=col_ann
                                              ,annotation_row=row_ann
                                              ,cluster_rows=cluster_rows
                                              ,cluster_cols=cluster_cols
                                              ,clustering_method='ward.D2'
                                              ,fontsize=8,border_color=NA
                                              ,col=colors_ht
                                              ,breaks=breaks
                                              ,legend_breaks=legend_breaks
                                              ,legend_labels=legend_labels
                                              ,treeheight_row=0
                                              ,scale="none"
                                              ,show_colnames=show_colnames
                                              ,show_rownames=show_rownames
                                              ,annotation_names_row=annotation_names_row
                                              ,fontsize_row=8
                                              ,plot=FALSE))
  return(ht)
}


#' Function to create boxplots
#' @export
#' @examples
#' ## Placeholder Example ##
create_bxplot<-function(df, continuous_var, discrete_var, fill_var = NULL,
                          padj = "BH", colors = NULL,
                          title = NULL, show_pairwise = FALSE,
                          show_median_line = TRUE,
                          label_y_test_pos_pct = 0.10,
                          angle_x_var = 0, legend_bxplt=FALSE,
                          legend_jitter=FALSE, fill_label=NULL,
                          ylab=NULL,xlab=NULL) {

  # Set fill label if not provided
  if(!is.null(fill_var) && is.null(fill_label)) fill_label <- fill_var

  # Ensure discrete variables are factors
  df[[discrete_var]] <- factor(df[[discrete_var]])
  if(!is.null(fill_var)) df[[fill_var]] <- factor(df[[fill_var]])

  # Subset and drop unused levels
  df_sub <- df[!is.na(df[[discrete_var]]) & !is.na(df[[continuous_var]]), ]
  df_sub[[discrete_var]] <- droplevels(df_sub[[discrete_var]])
  if(!is.null(fill_var)) df_sub[[fill_var]] <- droplevels(df_sub[[fill_var]])

  # Colors
  fill_for_colors <- if(!is.null(fill_var)) fill_var else discrete_var
  n_fill_levels <- length(levels(df_sub[[fill_for_colors]]))
  if(is.null(colors)) colors <- rainbow(n_fill_levels)
  colors <- rep(colors, length.out = n_fill_levels)

  # Number of groups for global test (x-axis)
  groups <- levels(df[[discrete_var]])
  n_groups <- length(groups)

  # Initialize
  test_table <- NULL
  global_method <- NULL

  #Statistical tests
  if(!is.null(fill_var)) {
    if(show_pairwise) {
      # Pairwise within each discrete group
      test_table <- df_sub |>
        group_by(.data[[discrete_var]]) |>
        wilcox_test(as.formula(paste(continuous_var, "~", fill_var))) |>
        adjust_pvalue(method = padj) |>
        add_significance() |>
        ungroup()

      test_table <- test_table |>
        mutate(
          xmin = map_dbl(group1, ~ which(levels(df_sub[[fill_var]]) == .x)),
          xmax = map_dbl(group2, ~ which(levels(df_sub[[fill_var]]) == .x)),
          y.position = .y.position
        )
    } else {
      # One test per discrete group: compare fill_var levels within each group
      n_fill <- length(levels(df_sub[[fill_var]]))
      test_fn <- if(n_fill == 2) {
        function(d) wilcox_test(d, as.formula(paste(continuous_var, "~", fill_var)))
      } else {
        function(d) dunn_test(d, as.formula(paste(continuous_var, "~", fill_var)), p.adjust.method = padj)
      }

      test_table <- df_sub |>
        group_by(.data[[discrete_var]]) |>
        group_modify(~ as.data.frame(test_fn(.x))) |>
        ungroup() |>
        add_y_position(formula = as.formula(paste(continuous_var, "~", fill_var)), data = df_sub)

      global_method <- if(n_fill == 2) "wilcox.test" else "kruskal.test"
    }
  } else {
    # No fill_var: global test across discrete groups
    if(n_groups == 2){
      test_table <- wilcox_test(df_sub, as.formula(paste(continuous_var, "~", discrete_var)))
      test_table <- as.data.frame(test_table)
      global_method <- "wilcox.test"
    } else if(n_groups > 2){
      test_table <- dunn_test(df_sub, as.formula(paste(continuous_var, "~", discrete_var)), p.adjust.method = padj)
      test_table <- add_y_position(test_table)
      test_table <- as.data.frame(test_table)
      global_method <- "kruskal.test"
    }
  }

  # Dynamic label.y for global test
  y_max <- max(df_sub[[continuous_var]], na.rm = TRUE)
  y_min <- min(df_sub[[continuous_var]], na.rm = TRUE)
  global_label_y <- y_max + label_y_test_pos_pct * (y_max - y_min)

  # Build base plot
  bxplot <- ggplot(df_sub, aes(x = .data[[discrete_var]],
                               y = .data[[continuous_var]],
                               fill = if(!is.null(fill_var)) .data[[fill_var]] else .data[[discrete_var]])) +
    geom_boxplot(position = position_dodge(ifelse(!is.null(fill_var), 0.8, 0.9)),
                 outlier.shape = NA, show.legend = legend_bxplt) +
    geom_jitter(position = position_jitterdodge(dodge.width = ifelse(!is.null(fill_var), 0.8, 0.9),
                                                jitter.width = 0.5),
                alpha = 0.5, size = 1.5, show.legend = legend_jitter) +
    labs(y = ylab, x = xlab, fill = fill_label) +
    scale_fill_manual(values = colors) +
    theme_classic(base_size = 11) +
    theme(axis.text.x = element_text(angle = angle_x_var,
                                     hjust = ifelse(angle_x_var == 0, 0.5, ifelse(angle_x_var <= 45, 1, 0))),
          axis.text.y = element_text(size = 7),
          plot.title = element_text(size = 10),
          legend.position='bottom',plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"))

  # Add per-group test (fill_var defined, not pairwise)
  if(!is.null(fill_var) && !show_pairwise && !is.null(test_table) && nrow(test_table) > 0){

    disc_levels <- levels(df_sub[[discrete_var]])

    # Summarize to one label per discrete group (use overall p or min p)
    label_df <- test_table |>
      group_by(.data[[discrete_var]]) |>
      summarise(
        p_combined = {
          if(n_fill == 2) {
            # Wilcoxon: only one p value anyway
            formatC(min(p.adj), format = "g", digits = 3)
          } else {
            # Kruskal-Wallis: recompute global p from raw data per group
            group_data <- df_sub[df_sub[[discrete_var]] == cur_group()[[1]], ]
            kt <- kruskal.test(as.formula(paste(continuous_var, "~", fill_var)), data = group_data)
            formatC(kt$p.value, format = "g", digits = 3)
          }
        },
        .groups = "drop"
      ) |>
      mutate(
        label      = paste0(if(n_fill == 2) "Wilcoxon p = " else "Kruskal-Wallis p = ", p_combined),
        x          = as.numeric(factor(.data[[discrete_var]], levels = disc_levels)),
        y.position = y_max + label_y_test_pos_pct * (y_max - y_min)
      )

    bxplot <- bxplot +
      geom_text(data  = label_df,
                aes(x = x, y = y.position, label = label),
                inherit.aes = FALSE,
                size        = 3,
                hjust       = 0.5)
  }

  # Add pairwise comparisons (fill_var + show_pairwise)
  if(!is.null(fill_var) && show_pairwise && !is.null(test_table) && nrow(test_table) > 0){
    bxplot <- bxplot +
      stat_pvalue_manual(test_table,
                         label = "p.adj.signif",
                         xmin = "xmin", xmax = "xmax",
                         y.position = "y.position",
                         hide.ns = TRUE,
                         tip.length = 0.03,
                         step.increase = 0.08,
                         bracket.nudge.y = 0.5)
  }

  # Global test only when no fill_var
  if(is.null(fill_var) && !is.null(global_method)){
    bxplot <- bxplot +
      stat_compare_means(method = global_method,
                         label.x = n_groups/2,
                         label.y = global_label_y,
                         size = 3)
  }

  # Median line
  if(show_median_line){
    bxplot <- bxplot + geom_hline(yintercept = median(df_sub[[continuous_var]], na.rm = TRUE),
                                  linetype = "dashed", color = "black", linewidth = 0.5)
  }

  # Add title
  if(!is.null(title)) bxplot <- bxplot + ggtitle(title)

  return(list(test_table = test_table, bxplot = bxplot))
}

#' Function to wrap generation of multiple boxplots
#' @export
#' @examples
#' ## Placeholder Example ##
run_bxplots<-function(df, continuous_vars, discrete_vars, fill_var = NULL,
                        padj = "BH", colors = NULL,
                        disc_var_labels = NULL, show_pairwise = FALSE,
                        title_suffix = NULL, angle_x_var = 0,
                        label_y_test_pos_pct = 0.10, show_median_line = TRUE,
                        legend_bxplt = FALSE, legend_jitter = FALSE,
                        fill_label = NULL,ylab=NULL,xlab=NULL) {

  # Ensure vectors
  continuous_vars <- as.character(continuous_vars)
  discrete_vars <- as.character(discrete_vars)

  # Default discrete labels
  if(is.null(disc_var_labels)) disc_var_labels <- discrete_vars

  if(length(disc_var_labels) != length(discrete_vars)){
    stop("disc_var_labels must match length of discrete_vars")
  }

  results <- list()
  k <- 1

  for(cont_var in continuous_vars){
    for(i in seq_along(discrete_vars)){
      disc_var <- discrete_vars[i]
      disc_label <- disc_var_labels[i]

      # Safety checks
      if(!cont_var %in% colnames(df)) stop(paste("Column not found:", cont_var))
      if(!disc_var %in% colnames(df)) stop(paste("Column not found:", disc_var))
      if(!is.null(fill_var) && !fill_var %in% colnames(df)) stop(paste("Column not found:", fill_var))

      # Generate plot title
      plot_title <- if(is.null(title_suffix)) {
        paste(cont_var, "by", disc_label)
      } else {
        paste(cont_var, "by", disc_label, title_suffix)
      }

      # Call create_bxplot
      res <- create_bxplot(df = df,
                           continuous_var = cont_var,
                           discrete_var = disc_var,
                           fill_var = fill_var,
                           padj = padj,
                           colors = colors,
                           title = plot_title,
                           show_pairwise = show_pairwise,
                           show_median_line = show_median_line,
                           label_y_test_pos_pct = label_y_test_pos_pct,
                           angle_x_var = angle_x_var,
                           legend_bxplt = legend_bxplt,
                           legend_jitter = legend_jitter,
                           fill_label = fill_label,
                           ylab=ylab,
                           xlab=xlab
      )

      # Store result
      results[[k]] <- res
      names(results)[k] <- paste(cont_var, disc_label, sep = "_")
      k <- k + 1
    }
  }

  return(results)
}

#' Function for flexible interleaving multi-plot results (customize order of tables and plots with optional protein/group naming)
#' @export
#' @examples
#' ## Placeholder Example ##
interleave_results<-function(list_of_lists,protein_names=NULL,
                             group_names=NULL,include_protein=TRUE,
                             include_group=TRUE,labels=NULL,suffixes=NULL) {

  n_lists<-length(list_of_lists)

  # Ensure all elements are lists
  if(!all(sapply(list_of_lists,is.list))){
    stop("All elements of list_of_lists must be lists")
  }

  # Ensure all lists have same length (truncate if necessary)
  list_lengths<-sapply(list_of_lists, length)
  min_len<-min(list_lengths)
  if(length(unique(list_lengths)) !=1){
    warning("Lists have unequal lengths. Interleaving will truncate to shortest length")
  }

  # Default labels/suffixes
  if(is.null(labels)) labels<-rep("",n_lists)
  if(is.null(suffixes)) suffixes<-rep("",n_lists)
  if(length(labels) != n_lists) stop("Length of labels must match number of lists")
  if(length(suffixes) != n_lists) stop("Length of suffixes must match number of lists")

  # Default protein/group names
  if(is.null(protein_names)) protein_names<-rep("",min_len)
  if(is.null(group_names)) group_names<-rep("",min_len)

  combined<-vector("list",length=n_lists * min_len)
  names_combined<-character(n_lists * min_len)

  idx<-1
  for( i in seq_len(min_len)){
    for( j in seq_len(n_lists)){
      combined[[ idx ]]<-list_of_lists[[ j ]][[ i ]]

      parts<-c()
      if(include_protein && protein_names[i ] !="") parts<-c(parts,protein_names[ i ])
      if(include_group && group_names[ i ] != "") parts<-c(parts,group_names[ i ])
      base_name <- paste(parts, collapse="_")

      names_combined[ idx ] <- paste0(labels[ j ], base_name, suffixes[ j ])
      idx<-idx+1
    }
  }

  names(combined)<-names_combined
  return(combined)
}


