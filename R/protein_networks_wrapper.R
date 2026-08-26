### REMEMBER TO OPEN CYTOSCAPE SOFTWARE BEFORE RUNNING THIS WRAPPER ###

#' Function to wrap protein networks pipeline
#' @export
#' @examples
#' ## Placeholder Example ##
build_string_network<-function(input_table,mode=c("multi","single"),
                               first_neighbors=FALSE,hide_singletons=TRUE,
                               show_first_neighbors=FALSE,
                               neighbor_min_connections=string_neighbor_min_connections,
                               neighbor_score_threshold=string_neighbor_score_threshold,
                               net_proteins=NULL,wait_for_layout_adj=TRUE){

  mode<-match.arg(mode)

  #Get STRING data based on mode
  string_data<- if (mode=="single") {
    get_single_ptn_network(input_table,first_neighbors=first_neighbors,
                           neighbor_min_connections=neighbor_min_connections,
                           neighbor_score_threshold=neighbor_score_threshold)
  } else {
    get_stringdb_data(input_table,first_neighbors=first_neighbors,
                      neighbor_min_connections=neighbor_min_connections,
                      neighbor_score_threshold=neighbor_score_threshold)
  }

  #Filter out first neighbors from network if requested
  nodes_display<-if (!show_first_neighbors && first_neighbors) {
    string_data$nodes |> dplyr::filter(!is_neighbor)
  } else {
    string_data$nodes
  }
  edges_display<-if (!show_first_neighbors && first_neighbors) {
    #Keep only edges between non-neighbor nodes
    non_neighbor_ids<-nodes_display$id
    string_data$edges |>
      dplyr::filter(source %in% non_neighbor_ids & target %in% non_neighbor_ids)
  } else {
    string_data$edges
  }

  #Create base network
  create_base_network(nodes_display,edges_display,hide_singletons=hide_singletons)

  #Store full node and edge table in case of first neighbors=TRUE
  string_data$nodes_display<-nodes_display
  string_data$edges_display<-edges_display

  #Store net membership in node table
  string_data$nodes$net_group<-NA
  if (!is.null(net_proteins)) {
    #Auto-name networks if unnamed
    if (is.null(names(net_proteins)) || any(names(net_proteins)=="")) {
      names(net_proteins)<-paste0("Net",seq_along(net_proteins))
    }

    #Assign net_group in node table
    for (net_name in names(net_proteins)) {
      string_data$nodes$net_group[string_data$nodes$label %in% net_proteins[[net_name]]]<-net_name
    }
    string_data$auto_color_edges<-TRUE
  } else {
    string_data$auto_color_edges<-FALSE
  }

  #Modify node shape conditionally
  if (first_neighbors) {
    string_data$nodes<-string_data$nodes |>
      dplyr::mutate(node_shape_group=dplyr::case_when(is_ptm~"ptm",is_neighbor~ "neighbor",TRUE~ "core"))
  } else {
    string_data$nodes<-string_data$nodes |>
      dplyr::mutate(node_shape_group=dplyr::case_when(is_ptm~"ptm",TRUE~"core"))
  }

  #Attach node data
  loadTableData(string_data$nodes,table="node",table.key.column="id",data.key.column="id")

  #Initial layout and time to adjust it
  apply_base_visual_style(first_neighbors=first_neighbors)
  layoutNetwork('force-directed')  #recalculates layout for better spread after hub sizing

  if(wait_for_layout_adj){
    message("Adjust network layout manually if needed.")
    readline("Press ENTER when finished...")
  }

  return(string_data)
}

#' Function to set constrasts within the protein networks pipeline
#' @export
#' @examples
#' ## Placeholder Example ##
set_contrast<-function(contrast,filename,edge_palette=NULL,wait_for_legend=TRUE
                       ,expr_col=NULL,width_col=NULL,class_col=NULL){

  if(length(getNetworkList())==0){
    stop("No Cytoscape network detected. Run build_string_network() first.")
  }

  nodes<-getTableColumns("node")

  #Automatically enable edge coloring if any net group is present
  color_edges<-any(!is.na(nodes$net_group))

  edge_sets<-NULL
  if(color_edges){
    net_labels<-na.omit(unique(nodes$net_group))
    net_proteins<-lapply(net_labels, function(label) nodes$label[nodes$net_group==label])
    names(net_proteins)<-net_labels

    if(length(net_proteins) >= 1){
      edge_sets<-highlight_edge_sets(net_proteins)
    } else {
      warning("No protein sets found. Skipping edge highlighting.")
      color_edges<-FALSE
    }
  }

  #Apply visual contrast
  apply_visual_contrast(contrast,color_edges=color_edges
                        ,edge_sets=edge_sets,edge_palette=edge_palette
                        ,expr_col=expr_col,width_col=width_col,class_col=class_col)

  if(wait_for_legend){
    message("Add Legends manually in Cytoscape if needed.")
    readline("Press ENTER when finished...")
  }

  fitContent()
  exportImage(filename=filename,type="SVG")
  message("Exported: ",paste0(contrast,"_network.svg"))
}

