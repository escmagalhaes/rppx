### Protein networks functions ###

#' Function to extract networks with StringDB
#' @export
#' @examples
#' ## Placeholder Example ##
get_stringdb_data<-function(input_table,first_neighbors=FALSE
                            ,neighbor_min_connections=string_neighbor_min_connections
                            ,neighbor_score_threshold=string_neighbor_score_threshold
                            ,stringdb_version=string_version
                            ,stringdb_score_threshold=string_score_threshold
                            ,stringdb_cache_path=string_cache_path
                            ,stringdb_species=string_species) {

  #In case of running the same network multiple times
  cache_file<-file.path(stringdb_cache_path
                        ,paste0("string_network_",stringdb_species,"_"
                                ,stringdb_score_threshold,"_"
                                ,ifelse(first_neighbors
                                        ,paste0("neighbors_",neighbor_min_connections
                                                ,"_",neighbor_score_threshold, "_"),"")
                                ,digest::digest(sort(unique(input_table$final_names)))
                                ,".rds"))

  if(file.exists(cache_file)){
    message("Loading cached STRING network...")
    return(readRDS(cache_file))
  }

  #Initialize STRINGdb
  message("Initializing STRINGdb...")
  string_db<-STRINGdb$new(version=string_version,
                          species=string_species,
                          score_threshold=string_score_threshold,
                          input_directory=string_cache_path)

  #Extract unique protein names from input
  proteins<-data.frame(final_names=unique(na.omit(input_table$final_names)))

  #Query STRINGdb with unique final_names
  message("Mapping proteins to STRING IDs...")
  mapped_proteins<-string_db$map(proteins,"final_names",removeUnmappedRows=TRUE)
  mapped<-input_table |> inner_join(mapped_proteins,by="final_names",relationship="many-to-many")

  #Add first neighbors if TRUE
  if (first_neighbors) {
    message("Retrieving first neighbors...")
    neighbor_map<-lapply(mapped_proteins$STRING_id, function(id) {
      edges_sub<-string_db$get_interactions(c(id,string_db$get_neighbors(id))) |>
        dplyr::filter(combined_score >= neighbor_score_threshold) |>
        dplyr::filter(from==id | to==id)
      data.frame(source_STRING_id=id,neighbor_STRING_id=c(edges_sub$from,edges_sub$to)) |>
        dplyr::filter(neighbor_STRING_id !=id)
    }) |> dplyr::bind_rows()

    neighbor_connectivity<-neighbor_map |>
      dplyr::filter(!neighbor_STRING_id %in% mapped_proteins$STRING_id) |>
      dplyr::group_by(neighbor_STRING_id) |>
      dplyr::summarise(n_connections=dplyr::n_distinct(source_STRING_id),.groups="drop") |>
      dplyr::filter(n_connections >= neighbor_min_connections)

    #Build neighbor_info from connectivity results only — no redundant get_neighbors call
    neighbor_info<-string_db$get_proteins() |>
      dplyr::filter(protein_external_id %in% neighbor_connectivity$neighbor_STRING_id) |>
      dplyr::select(protein_external_id,preferred_name) |>
      dplyr::rename(STRING_id=protein_external_id,final_names=preferred_name) |>
      dplyr::left_join(neighbor_connectivity,by=c("STRING_id"="neighbor_STRING_id")) |>
      dplyr::mutate(is_neighbor=TRUE,ptm=NA_character_,ptm_type="total"
                    ,base_name=final_names,protein=final_names)

    #Expand mapped with neighbor rows
    mapped<-dplyr::bind_rows(mapped |> dplyr::mutate(is_neighbor=FALSE),neighbor_info) #NA for all correlation/expression columns

  } else {
    mapped<-mapped |> dplyr::mutate(is_neighbor=FALSE)
  }

  #Create nodes table with PTM duplication
  nodes<-mapped |> dplyr::mutate(
    ptm_group = ifelse(is.na(ptm),"total",ptm),
    id        = ifelse(ptm_group == "total",final_names,paste0(final_names,".",ptm)),
    label     = final_names) |>
    dplyr::filter(!is.na(label)) #Remove nodes related to histone marks

  #Introduce transparency bias for PTM vs total
  if (first_neighbors) {
    nodes$node_alpha<-dplyr::case_when(nodes$is_neighbor~180
                                       ,nodes$ptm_group=="total"~220,TRUE~255)
  } else {
    nodes$node_alpha<-ifelse(nodes$ptm_group=="total",220,255)
    }
  nodes$is_ptm<-grepl("\\.",nodes$id)

  #Reorder dataframe
  nodes<-nodes |> select(id,label,ptm_group,STRING_id,node_alpha,is_ptm,everything())

  #Get interactions from STRINGdb
  message("Retrieving interactions from STRINGdb (filtered by score)...")
  interactions_raw<-string_db$get_subnetwork(unique(na.omit(mapped$STRING_id)))

  #Convert igraph to edge table
  edges<-igraph::as_data_frame(interactions_raw,what="edges") |>
    dplyr::rename(STRING_source=from
           ,STRING_target=to
           ,score=combined_score) |>
    dplyr::filter(score >= string_score_threshold)

  #Helper Function to to duplicated PTM edges and adjust if necessary
  correct_ptm_edges<-function(nodes,edges_original) {

    #Start with total protein edges (STRING_id based)
    edges_total<-edges_original |> transmute(source=STRING_source,target=STRING_target,score=score)

    #Join STRING_id to id
    edges_expanded<-edges_total |>
      dplyr::left_join(nodes |> select(id,STRING_id),by=c("source"="STRING_id"),relationship="many-to-many") |>
      dplyr::rename(source_node=id) |>
      dplyr::left_join(nodes |> select(id,STRING_id),by=c("target"="STRING_id"),relationship="many-to-many") |>
      dplyr::rename(target_node=id) |>
      dplyr::filter(!is.na(source_node) & !is.na(target_node)) |>
      dplyr::transmute(source=source_node,target=target_node,score=score) |> distinct()
    return(edges_expanded)
  }

  #Apply helper function to adjust edges
  edges_clean<-correct_ptm_edges(nodes,edges)

  ### SAFETY CHECKS ###
  #Canonicalize original edges (undirected)
  orig_pairs<-paste0(pmin(edges$STRING_source,edges$STRING_target),"_",
                     pmax(edges$STRING_source,edges$STRING_target))
  unique_orig_pairs<-unique(orig_pairs)

  #Canonicalize expanded edges (map back to STRING_id for total nodes only)
  total_nodes<-nodes |> dplyr::filter(ptm_group=="total") |> dplyr::select(id,STRING_id)

  expanded_total_edges<-edges_clean |>
    dplyr::inner_join(total_nodes,by=c("source"="id"),relationship="many-to-many") |>
    dplyr::rename(source_string=STRING_id) |>
    dplyr::inner_join(total_nodes,by=c("target"="id"),relationship="many-to-many") |>
    dplyr::rename(target_string=STRING_id)

  expanded_pairs<-paste0(pmin(expanded_total_edges$source_string
                              ,expanded_total_edges$target_string),"_"
                         ,pmax(expanded_total_edges$source_string
                               ,expanded_total_edges$target_string))

  #Check that no original edge disappeared
  missing_pairs<-setdiff(unique_orig_pairs,unique(expanded_pairs))
  if(length(missing_pairs) > 0){
    warning(length(missing_pairs)," STRING edges do not have PTM counterparts (expected when PTMs are missing).")
  }
  #Check that all expanded edges connect valid nodes
  invalid_edges<-edges_clean |>
    dplyr::filter(!source %in% nodes$id | !target %in% nodes$id)
  if (nrow(invalid_edges) > 0) {
    stop("ERROR: edges_clean contains nodes not present in nodes table.")
  }
  #Check that expansion did not shrink unexpectedly
  if (nrow(edges_clean) < length(unique_orig_pairs)) {
    stop("ERROR: Expanded edge count smaller than unique original edges.")
  }

  message("All Safety checks performed. Nodes and edges tables ready. Network structure is valid. Network can be created in Cytoscape.\n")
  message("Total nodes: ", sum(nodes$ptm_group=="total"),
          "; PTM nodes: ", sum(nodes$ptm_group!="total"))

  saveRDS(list(nodes=nodes,edges=edges_clean),cache_file)
  return(list(nodes=nodes,edges=edges_clean))

}

#' Function to create base network
#' @export
#' @examples
#' ## Placeholder Example ##
create_base_network<-function(nodes,edges
                              ,title="Protein Network"
                              ,collection="Proteomics"
                              ,hide_singletons=TRUE) {

  message("Connecting to Cytoscape...")

  #Connect Cytoscape
  cytoscapePing()
  cytoscapeVersionInfo()

  #Delete previous networks
  deleteAllNetworks()

  #Set seed for reproducibility
  set.seed(123)

  #(optional) Hide singletons from network
  if(!"name" %in% colnames(nodes)) nodes$name<-nodes$id #safety check, cytoscape require 'name' column
  if(hide_singletons){
    connected_nodes<-unique(c(edges$source,edges$target))
    nodes<-nodes[nodes$id %in% connected_nodes,]
  }
  if(hide_singletons){
    message(length(setdiff(nodes$id,unique(c(edges$source,edges$target)))),
            " singleton nodes removed before network creation.")
  }

  message("Creating network...")

  #Create network in cytoscape
  createNetworkFromDataFrames(nodes=nodes,edges=edges
                              ,title=title,
                              collection=collection
                              ,layout.name="none")

  #Wait for Cytoscape to finish
  while(TRUE){
    res<-tryCatch({
      nets<-getNetworkList()
      title %in% nets
    },error=function(e) FALSE)

    if(res) break
    Sys.sleep(1)
  }

  #Load node table
  #loadTableData(nodes,data.key.column="id",table="node",table.key.column="id")

  clearSelection()
  fitContent()

  message("Network successfully created.")
}

#' Function to ajust nodes according to protein connections (used inside apply_base_visual_style)
#' @export
#' @examples
#' ## Placeholder Example ##
highlight_hubs<-function(style_name="dataStyle",min_size=300,max_size=500,top_percent=0.1){

  message("Calculating node connectivity (degree)...")

  nodes<-getTableColumns("node")
  edges<-getTableColumns("edge")

  if(!all(c("source","target") %in% colnames(edges))){
    stop("Edge table missing source/target columns.")
  }

  #Compute node degree
  degree_table<-table(c(edges$source,edges$target))
  degree_df<-data.frame(id=names(degree_table),degree=as.numeric(degree_table),stringsAsFactors=FALSE)

  #Send degree back to Cytoscape
  loadTableData(degree_df,data.key.column="id",table.key.column="id",table="node")

  message("Degree column added to node table.")

  #Continuous node size mapping
  setNodeSizeMapping("degree",c(min(degree_df$degree),max(degree_df$degree)),c(min_size,max_size),mapping.type="c",style.name=style_name)

  #Detect protein hubs
  cutoff<-quantile(degree_df$degree,1-top_percent)
  hubs<-degree_df$id[degree_df$degree>=cutoff]

  message(length(hubs)," hub proteins detected: ",paste(hubs,collapse=', '))
}

#' Function to map STRING confidence to edge width
#' @export
#' @examples
#' ## Placeholder Example ##
apply_edge_confidence<-function(style_name="dataStyle",min_width=2,max_width=6){

  message("Applying STRING confidence edge widths...")

  edges<-getTableColumns("edge")
  if(!"score" %in% colnames(edges)){
    stop("Edge 'score' column not found.")
  }
  score_vals<-edges$score
  score_vals<-score_vals[!is.na(score_vals)]
  setEdgeLineWidthMapping("score",c(min(score_vals), max(score_vals)),
                          c(min_width, max_width),mapping.type="c",style.name=style_name)

  message("Edge confidence mapping applied.")
}

#' Function to apply base visual style
#' @export
#' @examples
#' ## Placeholder Example ##
apply_base_visual_style<-function(style_name="dataStyle",first_neighbors=FALSE) {

  message("Creating base visual style...")

  createVisualStyle(style_name)
  setVisualStyle(style_name)

  layoutNetwork('force-directed defaultSpringLength=150 defaultSpringCoefficient=0.0000001')

  setNodeShapeDefault("ELLIPSE",style_name)
  setNodeLabelMapping("id", style_name)
  setNodeColorDefault("#AAAAAA",style_name)
  setNodeSizeDefault(300, style_name)
  setNodeFontSizeDefault(200,style_name)
  setEdgeLineWidthDefault(1,style_name)
  setNodeBorderWidthDefault(10,style_name)

  #Apply node shape and transparency conditionally
  if (first_neighbors) {
    setNodeShapeMapping("node_shape_group",c("ptm","neighbor","core")
                        ,c("OCTAGON","DIAMOND","ELLIPSE"),style.name=style_name)
    #setNodeOpacityBypass("node_alpha",c(180,220,255),c(180,220,255),mapping.type="c"
                          #,style.name=style_name)
  } else {
    setNodeShapeMapping("node_shape_group",c("ptm","core"),c("DIAMOND","ELLIPSE")
                        ,style.name=style_name)
    #setNodeOpacityBypass("node_alpha",c(220,255),c(220,255),mapping.type="c"
                          #,style.name=style_name)
  }


  #Highlight hubs by increasing their node size
  highlight_hubs(style_name)

  #Apply STRING confidence edge width
  apply_edge_confidence(style_name)

  message("Base visual style applied.")
}

#' Function to color edges from specific networks
#' @export
#' @examples
#' ## Placeholder Example ##
highlight_edge_sets<-function(net_proteins){

  n_nets<-length(net_proteins)
  message("Highlighting edges for ",n_nets," networks...")

  #Assign default betwork names if the list is unnamed
  if(is.null(names(net_proteins)) || any(names(net_proteins)=="")){
    names(net_proteins)<-paste0("Net", seq_along(net_proteins))
  }

  net_edges<-getTableColumns("edge")
  net_nodes<-getTableColumns("node")
  node_map<-net_nodes[,c("id","label")]

  edges_named<-net_edges |>
    left_join(node_map,by=c("source"="id"))  |> dplyr::rename(source_name=label)  |>
    left_join(node_map,by=c("target"="id"))  |> dplyr::rename(target_name=label)

  edges_within<-list()

  # Compute within-set edges
  for(net_name in names(net_proteins)){
    proteins<-net_proteins[[net_name]]
    idx<-edges_named$source_name %in% proteins & edges_named$target_name %in% proteins
    edges_within[[net_name]]<-edges_named$name[idx]
  }

  # Compute intersection edges only if more than one network
  edges_inter<-character(0)
  if(n_nets > 1){
    idx_inter <- rep(FALSE, nrow(edges_named))
    for(i in 1:(n_nets-1)){
      for(j in (i+1):n_nets){
        idx_i<-edges_named$source_name %in% net_proteins[[ i ]] & edges_named$target_name %in% net_proteins[[ j ]]
        idx_j<-edges_named$source_name %in% net_proteins[[ j ]] & edges_named$target_name %in% net_proteins[[ i ]]
        idx_inter<-idx_inter | idx_i | idx_j
      }
    }
    edges_inter<-edges_named$name[idx_inter]
  }

  #Default edges
  edges_default<-edges_named$name[!edges_named$name %in% c(unlist(edges_within),edges_inter)]

  message("Edges categorized:")
  for(names_edges in names(edges_within)){
    message(names_edges, ": ",length(edges_within[[ names_edges ]]))
  }
  message("Intersection: ",length(edges_inter))
  message("Default: ",length(edges_default))

  return(list(
    edges_default=edges_default,
    edges_within=edges_within,
    edges_inter=edges_inter
  ))
}

#' Function to apply specific visual style according to contrast
#' @export
#' @examples
#' ## Placeholder Example ##
apply_visual_contrast<-function(contrast,style_name="dataStyle"
                                ,border_colors=c("deepskyblue2","firebrick1")
                                ,min_border_width=50,max_border_width=100
                                ,default_edge_color='gray70'
                                ,color_edges=FALSE,edge_sets=NULL
                                ,edge_palette=NULL,edge_width_highlight=20
                                ,expr_col=NULL,width_col=NULL,class_col=NULL){

  message("Applying contrast: ",contrast)

  #Define columns depending on contrast
  expr_col<-if (!is.null(expr_col)) expr_col else paste0("mean_expr_",contrast)
  width_col<-if (!is.null(width_col)) width_col else paste0("abs_diff_",contrast)
  class_col<-if (!is.null(class_col)) class_col else paste0("class_",contrast)

  nodes<-getTableColumns(table="node")

  #Define node colors (expression values)
  expr_values<-nodes[[expr_col]]
  expr_values<-expr_values[is.finite(expr_values)]

  #Safety check
  if (length(expr_values)==0) {
    warning("No valid expression values for ",expr_col,". Skipping color mapping.")
  } else {

  #Helper function to adjust Node Color scaling for visualization
  scale_color_breaks<-function(value){
    if(min(value)>=-2 & max(value)<=2){
      node_breaks=seq(-2.4,2.4,0.4)
    }else if(min(value)>=-2){
      node_breaks=c(seq(-2.4,2,0.4),max(value))
    }else if(max(value)<=2){
      node_breaks=c(min(value),seq(-2,2.4,0.4))
    }else{
      node_breaks=c(min(value),seq(-2,2,0.4),max(value))
    }
    return(node_breaks)
  }
  node_breaks<-scale_color_breaks(expr_values)
  node_colors<-jet.colors(length(node_breaks))

  #Attach color node mapping
  setNodeColorMapping(
    expr_col,
    node_breaks,
    node_colors,
    style.name=style_name,
    mapping.type="c"
    )
  }

  #Define border widths (scale according to mode design)
  width_values<-nodes[[width_col]]
  width_values<-width_values[is.finite(width_values)]

  #Safety check
  if (length(width_values)==0) {
    warning("No valid width values for ",width_col,". Skipping width mapping.")
  } else {

    #Attach border width mapping
    setNodeBorderWidthMapping(
      width_col,
      c(min(width_values),max(width_values)),
      c(min_border_width,max_border_width),
      style.name=style_name,
      mapping.type="c"
      )
  }

  #Define border color
  #Safety check
  if (is.null(nodes[[class_col]]) || length(nodes[[class_col]])==0 || all(is.na(nodes[[class_col]]))) {
    warning("Skipping class_col mapping: ",class_col," is NULL/empty/all NA")
  } else {
    #Convert class_col to factor
    nodes[[class_col]]<-factor(nodes[[class_col]])
    class_levels<-levels(nodes[[class_col]])

    #Ensure enough colors
    if (length(border_colors) < length(class_levels)) {
      border_colors_use<-rep(border_colors,length.out=length(class_levels))
    } else {
      border_colors_use<-border_colors
    }

  #Attach border color mapping
  setNodeBorderColorMapping(
    class_col,
    class_levels,
    border_colors,
    style.name=style_name,
    mapping.type="d"
    )
  }

  #Optional edge coloring
  if(color_edges && !is.null(edge_sets)){
    n_nets<-length(edge_sets$edges_within)
    net_names<-names(edge_sets$edges_within)

    # Use user-provided palette or generate default
    if(is.null(edge_palette)){
      if(n_nets <= 12){
        edge_palette<-RColorBrewer::brewer.pal(max(n_nets,3),"Set3")
      } else {
        edge_palette<-grDevices::rainbow(n_nets)
      }
    }

    # Ensure enough colors for within-network + 1 intersection color
    if(length(edge_palette) < n_nets + 1){
      edge_palette<-rep(edge_palette,length.out=n_nets + 1)
    }

    # Within-network edges: first N colors
    for(i in seq_along(net_names)){
      edges<-edge_sets$edges_within[[net_names[ i ]]]
      setEdgeColorBypass(edge.names=edges,new.colors=edge_palette[ i ])
      setEdgeLineWidthBypass(edge.names=edges,new.width=edge_width_highlight)
    }

    # Intersection edges: next color in the palette
    if(length(edge_sets$edges_inter) > 0){
      intersection_color<-edge_palette[n_nets + 1]
      setEdgeColorBypass(edge.names=edge_sets$edges_inter,new.colors=intersection_color)
      setEdgeLineWidthBypass(edge.names=edge_sets$edges_inter,new.width=edge_width_highlight)
    }

    #Default edges
    setEdgeColorBypass(edge.names=edge_sets$edges_default,new.colors=default_edge_color)
    setEdgeLineWidthBypass(edge.names=edge_sets$edges_default,new.width=edge_width_highlight)

    message("Edges highlighted successfully for ",n_nets," networks.")
  }
  message("Contrast '",contrast,"' applied successfully.")
}

#' Function to get a single protein network
#' @export
#' @examples
#' ## Placeholder Example ##
get_single_ptn_network<-function(protein_name,
                                 first_neighbors=TRUE,
                                 neighbor_min_connections=1,
                                 neighbor_score_threshold=string_neighbor_score_threshold,
                                 stringdb_version=string_version,
                                 stringdb_score_threshold=string_score_threshold,
                                 stringdb_cache_path=string_cache_path,
                                 stringdb_species=string_species) {

  #Clean protein name
  clean<-clean_protein_names(protein_name,strict=FALSE)

  #Build minimal input_table
  input_table<-data.frame(final_names=clean$clean_names,original=protein_name
                          ,base_name=protein_name,ptm=NA_character_,stringsAsFactors=FALSE)

  get_stringdb_data(input_table,
                    first_neighbors=first_neighbors,
                    neighbor_min_connections=neighbor_min_connections,
                    neighbor_score_threshold=neighbor_score_threshold,
                    stringdb_version=stringdb_version,
                    stringdb_score_threshold=stringdb_score_threshold,
                    stringdb_cache_path=stringdb_cache_path,
                    stringdb_species=stringdb_species)

}

#' Function to compare networks against reference
#' @export
#' @examples
#' ## Placeholder Example ##
compare_ref_network<-function(ptn_net_list,ref_ptn,first_neighbors=FALSE
                              ,neighbor_min_connections=1
                              ,string_neighbor_score_threshold=700){

  ##REMEMBER TO OPEN CYTOSCAPE SOFTWARE##
  dataset_networks<-list()
  for (nm in names(ptn_net_map)) {
    message("\nBuilding network: ",nm,"...\n")
    ptn_net<-build_string_network(ptn_net_map[[ nm ]],mode="multi",
                                  first_neighbors=first_neighbors,
                                  neighbor_min_connections=neighbor_min_connections,
                                  neighbor_score_threshold=string_neighbor_score_threshold)

    #Store string data
    dataset_networks[[ nm ]]<-ptn_net
  }

  #Helper function to compute jaccard
  jaccard<-function(a, b) {
    a<-unique(na.omit(a))
    b<-unique(na.omit(b))
    length(intersect(a, b)) / length(union(a, b))
  }

  #Single protein network with first neighbors to get reference network
  single_ptn_data<-get_single_ptn_network(protein_name=ref_ptn
                                          ,first_neighbors=TRUE     #always TRUE for this one
                                          ,neighbor_min_connections=neighbor_min_connections
                                          ,neighbor_score_threshold=string_neighbor_score_threshold)

  #Jaccard per dataset network
  jaccard_results<-lapply(names(dataset_networks), function( nm ) {
    dataset_proteins<-unique(dataset_networks[[ nm ]]$nodes$label)
    reference_proteins<-unique(single_ptn_data$nodes$label)
    data.frame(dataset=nm
               ,jaccard=jaccard(dataset_proteins,reference_proteins)
               ,n_overlap=length(intersect(dataset_proteins,reference_proteins))
               ,n_dataset=length(dataset_proteins)
               ,n_reference=length(reference_proteins))
  }) |> dplyr::bind_rows()

  #Get overlapping proteins
  overlap_ptn<-lapply(names(dataset_networks), function( nm ) {
    dataset_proteins<-unique(dataset_networks[[ nm ]]$nodes$label)
    reference_proteins<-unique(single_ptn_data$nodes$label)
    intersect(dataset_proteins,reference_proteins)
  })
  names(overlap_ptn)<-names(ptn_net_map)

  return(list(jaccard=jaccard_results,overlaps=overlap_ptn))
}


