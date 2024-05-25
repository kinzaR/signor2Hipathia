suppressPackageStartupMessages(library("org.Hs.eg.db"))
suppressPackageStartupMessages(library(httr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(dplyr))

getPathwayData<-function(api_url="https://signor.uniroma2.it/getPathwayData.php?description&relations"){
  # Make a GET request to the API
  response <- GET(api_url)
  # Check if the request was successful (status code 200)
  if (httr::status_code(response) == 200) {
    # Read the TSV data from the response content
    tsv_data <- content(response, as = "text", encoding = "UTF-8")
    # Find the index where "pathway_id" occurs
    index <- gregexpr("pathway_id", tsv_data) %>% .[[1]] 
    # Split the TSV data into two parts and Convert it to a dataframe
    paths_info <- substr(tsv_data, 1, index - 1) %>% read_tsv
    paths_relations <- substr(tsv_data, index, nchar(tsv_data)) %>% read_tsv
    return(list(paths_info=paths_info, paths_relations=paths_relations))
  }
  stop("Failed to retrieve pathway data. Check the API endpoint and parameters.")
}
getSignorData<-function(api_url){
  if(nchar(api_url)==0) stop("Empty URL")
  # Make a GET request to the API
  response <- GET(paste0("https://signor.uniroma2.it/",api_url))
  # Check if the request was successful (status code 200)
  if (httr::status_code(response) == 200) {
    # Read the TSV data from the response content and Convert it to a dataframe
    data_df <- content(response, as = "text", encoding = "UTF-8") %>% read_tsv(show_col_types = F,skip_empty_rows = T) %>% select_if(~sum(!is.na(.)) > 0)
    return(data_df)
  }
  stop("Failed to retrieve data. Check the API endpoint and parameters.")
}

write_Sif_Att_FromRow <- function(pathway_tsv, spe="hsa",
                          activation_effect = c("up-regulates activity",
                                                              "up-regulates quantity by expression",
                                                              "up-regulates",
                                                              "up-regulates quantity",
                                                              "up-regulates quantity by stabilization"),
                          inhibition_effect = c("down-regulates",
                                                 "down-regulates quantity by repression",
                                                 "down-regulates quantity",
                                                 "down-regulates activity",
                                                 "down-regulates quantity by destabilization"),
                          avoided_effect = c("form complex", "unknown"),
                          proteinFamilies,
                          fusionProteins,
                          complexes,
                          output_folder=NULL,
                          verbose=F){
  
  log_file <- file.path(output_folder,"log.txt")
  path_id <- pathway_tsv$pathway_id %>% unique
  dir.create(output_folder, showWarnings = FALSE, recursive = T)
  write(x = c("* path_id : ", path_id), file = log_file, append = T, sep = "\t", ncolumns = 2)
  if(verbose) cat("Parsing the ",pathway_tsv$pathway_id %>% unique, "...")
  pathway_tsv <- add_hi_effect(pathway_tsv)# Change the interaction to hipathia standard 
  # Check for phenotypes:
  pheno_as_in <- pathway_tsv %>% filter(typea=="phenotype") %>% dim() %>% .[1] != 0
  stimulus_as_in <- pathway_tsv %>% filter(typea=="stimulus") %>% dim() %>% .[1] != 0
  # NOTE: stimulus_as_in? what doing?
  if(pheno_as_in){
    # stop("!------> Phenotype as In was found in ", pathway_tsv$pathway_id %>% unique)
    message("!------> Phenotype as In was found in ", pathway_tsv$pathway_id %>% unique)
    write(x = c("", "Not parsed ------> Phenotype(s) activates/inhibits other molecules."), file = log_file, append = T, sep = "\t", ncolumns = 2)
    return()
  }
  # browser() # check before taking off the phenotypes
  out_nodes<-setdiff(pathway_tsv$idb, pathway_tsv$ida)
  # # adding here gost ionteractions
  # gost_interactions_a <- pathway_tsv %>% filter(idb %in% out_nodes) %>% mutate(entityb=paste0(entitya,"*"),
  #                                                       typeb=typea,
  #                                                       idb=ida,
  #                                                       databaseb=databasea)
  # gost_interactions_b <- pathway_tsv %>% filter(idb %in% out_nodes) %>% mutate(entitya=paste0(entitya,"*"))
  # pathway_tsv<-pathway_tsv %>% filter(!idb %in% out_nodes) %>% rbind(gost_interactions_a,gost_interactions_b)
  if(length(out_nodes[!grepl("SIGNOR-PH",out_nodes)]) == 0){
    if(verbose) message("All effectors are in nodes! this will be taked into acount in futur versions!")# maybe duplicating the node if it is protein
    write(x = c("", "Not parsed ------> No final node found! All Signor effectors are activating/inhibiting other molecule(s). This cyclic issue will be addressed in future versions."), file = log_file, append = T, sep = "\t", ncolumns = 2)
    return()
  }
  effectors_before_pheno <- pathway_tsv %>% filter(idb %in% out_nodes[grepl("SIGNOR-PH",out_nodes)]) %>% select(ida) %>% pull()
  # check if all effectors are NOT last nodes in the pathway:
  lost_effectors <- pathway_tsv %>% filter(!(idb %in% out_nodes[grepl("SIGNOR-PH",out_nodes)]) & (ida %in% effectors_before_pheno)) %>% select(ida) %>% unique() %>% pull
  if(length(lost_effectors)>0)
    write(x = c("", paste0(lost_effectors, collapse = ","), "These effectors will not be considered as effectors in the parsed pathways because they have other interactions as a source and will not be detected as outnodes while performing HiPathia."),
        file = log_file, append = TRUE, sep = "\t", ncolumns = 3)
  pheno_as_out <- pathway_tsv %>% filter(typeb=="phenotype") %>% select(entitya, ida, entityb)
  stimulus_as_out <- pathway_tsv %>% filter(typeb=="stimulus") %>% select(entitya, ida, entityb)
    # stimulus_as_out ?
  if(length(stimulus_as_out$ida)>0)
    write(x = c("", paste0(stimulus_as_out$entityb, collapse = ","), "These stimuli are detected as targets. Should they be treated as phenotypes? Otherwise, they may introduce unnecessary noise during the mechanistic modeling process in the HiPathia method."),
        file = log_file, append = TRUE, sep = "\t", ncolumns = 3)
  pathway_tsv <- pathway_tsv %>% filter(typeb!="phenotype")
  # Here check again if still having out-nodes or is a cyclic pathway or closed-loop pathway?
  out_nodes_after_curation <- setdiff(pathway_tsv$idb, pathway_tsv$ida)
  if(length(out_nodes_after_curation)==0){
    if(verbose) message("This pathway has become a closed-loop pathway after technical curation.")
    write(x = c("", "Not parsed ------> This pathway has become a closed-loop pathway after technical curation!"), file = log_file, append = T, sep = "\t", ncolumns = 2)
    return()
  }
  
  ## change same label for different entity_id
  pathway_tsv <-  pathway_tsv %>% group_by(entitya) %>% mutate(n1 = n()) %>% group_by(entitya,ida) %>% mutate(n2=n()) %>% mutate(entitya_tmp = ifelse(n1!=n2, paste0(entitya,"(",ida,")"),entitya)) %>%
                                  group_by(entityb) %>% mutate(n1 = n()) %>% group_by(entityb,idb) %>% mutate(n2=n()) %>% mutate(entityb_tmp = ifelse(n1!=n2, paste0(entityb,"(",idb,")"),entityb)) %>%
                                  ungroup() %>% select(!c(n1,n2)) 
  ##
  graph<-list()
  graph$att <- get_hi_att(pathway_tsv,path_id, proteinFamilies, fusionProteins, complexes,spe, verbose, log_file) # get the att file of hipathia
  pathway_tsv <- pathway_tsv %>% rowwise %>% mutate(entitya = ifelse(!entitya %in% graph$att$label & entitya_tmp %in% graph$att$label,entitya_tmp,entitya))  %>% 
                              mutate(entityb = ifelse(!entityb %in% graph$att$label & entityb_tmp %in% graph$att$label,entityb_tmp,entityb))
  graph$sif <- get_hi_sif(pathway_tsv, graph$att ,verbose)
  # detect cyclic sub-paths in the whole pathway:
  unreachable_out_nodes<-get_unreachable_out_nodes(sif=graph$sif, verbose)
  if(length(unreachable_out_nodes) > 0){
    unreachable_out_nodes_labels <- graph$att %>%
      filter(ID %in% unreachable_out_nodes$name) %>%
      select(label) %>%
      pull()
    write(x = c("", paste0("Not parsed ------> Unreachable last nodes (effectors):", paste0(unreachable_out_nodes_labels, collapse = ","),". Hipathia requires all last nodes (effectors) to be reachable from the first nodes (receptors).")),
          file = log_file, append = TRUE, sep = "\t", ncolumns = 2)
    return()
  }
  graph <- add_layout2att(graph, verbose)
  # graph$att %>% filter(label %in% pheno_as_out$entitya) %%
  graph$stimulus_as_out <-  stimulus_as_out #only capturing the stimuli
  # Here: phynotypes has to be prepared as annot file
  graph$pheno_as_out <-  pheno_as_out %>% rowwise() %>% 
    left_join(graph$att, by = join_by("entitya"== "label")) %>% select(colnames(pheno_as_out),genesList) %>% unique() 
    # mutate(geneList = unique(graph$att$genesList[graph$att$label== entitya])) %>% unique()
  if(!is.null(output_folder)) write_sif_att_files(graph, spe, path_id, output_folder,verbose)
  write(x = c("", "Parsed!"), file = log_file, append = T, sep = "\t", ncolumns = 2)
  return(graph)
}
add_hi_effect<- function(pathway_tsv){
  pathway_tsv <- pathway_tsv %>% 
    filter(! effect %in% avoided_effect) %>%
    mutate(sign = case_when(effect %in% activation_effect ~ "activation",
                                          effect %in% inhibition_effect  ~ "inhibition",
                                          # effect %in% avoided_effect & !(idb %in% ida) ~ NA,
                                          .default = NA)
  )
  if(anyNA(pathway_tsv$sign)){
    cat("Warning: NA interacting were generated, please check ",unique(pathway_tsv$pathway_id)," Pathway.")
  }
  return(pathway_tsv)
}
resolve_signor_ids <- function(id, type, proteinFamilies, fusionProteins, complexes, verbose=F,log_file = NULL, na_action=T){
  new_id <- switch(type,
                   "protein" = getEntrezFromUniprot(id,verbose, log_file, na_action),
                   "proteinfamily" = getEntrezFromComposed(id,proteinFamilies,sep=",",verbose, log_file, na_action),
                   "fusion protein" = getEntrezFromComposed(id,fusionProteins,sep=",",verbose, log_file, na_action), 
                   "complex" = getEntrezFromComposed(id,complexes,sep=",/,",verbose, log_file, na_action),
                   "smallmolecule" = ifelse(na_action,NA,id) , # has to be changes fopr metabopathia
                   "mirna" = ifelse(na_action,NA,id),
                   "chemical" = ifelse(na_action,NA,id),# has to be changes fopr metabopathia
                   "antibody" = ifelse(na_action,NA,id),
                   "phenotype" = ifelse(na_action,NA,id),
                   "stimulus" = ifelse(na_action,NA,id)
                   )
  return(new_id)
}
get_hi_att <- function(pathway_tsv, path_id, proteinFamilies, fusionProteins, complexes,spe="hsa", verbose=F, log_file = NULL){
  if(verbose) print("Getting Att...")
  all_nodes <- tibble(
    node_id = c(pathway_tsv$ida,pathway_tsv$idb),
    node_label = c(pathway_tsv$entitya, pathway_tsv$entityb),
    node_label_tmp = c(pathway_tsv$entitya_tmp, pathway_tsv$entityb_tmp),
    node_type = c(pathway_tsv$typea, pathway_tsv$typeb)
  ) %>% unique()
  all_nodes <-all_nodes %>% 
    rowwise() %>% mutate(genesList=resolve_signor_ids(node_id, node_type, proteinFamilies, fusionProteins, complexes, verbose=verbose, log_file = log_file, na_action=T)) %>%
    mutate(id_indicator = resolve_signor_ids(node_id, node_type, proteinFamilies, fusionProteins, complexes, verbose=verbose, log_file = log_file, na_action=F)) 
  # %>% mutate(id_indicator= case_when(is.na(genesList) ~ node_id,.default = genesList))
  hi_ids <- sapply(all_nodes$id_indicator, function(g){
    if(grepl(",/,", x = g)) strsplit(x = g, split = ",/,")
    else g
  }) %>% unlist() %>% unique %>% na.omit %>% as_tibble(.) %>%
    mutate(ids=row.names(.))
  all_nodes<- all_nodes %>% rowwise() %>%mutate(hi_id = ifelse(id_indicator %in% hi_ids$value, hi_ids$ids[hi_ids$value==id_indicator],
                                                                      hi_ids[match(strsplit(x = id_indicator, split = ",/,")[[1]], hi_ids$value),2] %>% pull() %>% paste0(collapse = " ")
                                                                )) %>%
    mutate(hi_id= paste0("N-",spe,path_id,"-",hi_id))
  all_nodes<-all_nodes %>% group_by(node_label,hi_id) %>% mutate(n=n()) %>% mutate(node_label_X= ifelse(node_label != node_label_tmp & n==1,node_label_tmp,node_label)) %>% ungroup()
  # create tooltip:
  # all_nodes<-all_nodes %>% mutate(tooltip= paste0("https://signor.uniroma2.it/relation_result.php?id=", node_id))
  all_nodes<-all_nodes %>% mutate(tooltip= paste0("<a target=_blank href=https://signor.uniroma2.it/relation_result.php?id=",node_id,">",node_id,"</a> (",node_label,")"))  
  att <- tibble(ID = all_nodes$hi_id,
                label	=all_nodes$node_label_X,
                X=0,
                Y=0,
                color= "white",
                shape="rectangle",
                type= all_nodes$node_type,
                label.cex= 0.5,
                label.color ="black",
                width=46,
                height=17,
                genesList=all_nodes$genesList,
                tooltip= NA)%>% unique()
  if(any(duplicated(att$ID))) stop("Duplicated hipathia IDs in ATT!")
  # node type
  att <- att %>% mutate(type= case_when(grepl("protein", type) | type=="complex" ~ "gene",
                                          .default = "coumpound"))
  att <- att %>% mutate(shape= case_when(type=="gene" ~ "rectangle",
                                          .default = "circle"))
  # add tooltips
  att<-att %>% rowwise() %>% mutate(tooltip= all_nodes$tooltip[all_nodes$hi_id==ID][1]) # take the first label because sometmes the label is changed and in origin theres more tha one !
  return(att)
}
get_hi_sif <- function(pathway_tsv, att,verbose=F){
  if(verbose) print("Getting Sif...")
  sif<-pathway_tsv %>% rowwise() %>% mutate(hi_ida = att$ID[att$label==entitya]) %>% 
    mutate(hi_idb=att$ID[att$label==entityb]) %>% select(hi_ida,sign, hi_idb) %>% unique()
  return(sif)
}
getEntrezFromUniprot<- function(uniprot, verbose=F, log_file = NULL, na_action = T){
  if(uniprot=="/") return("/") # is a exception to keep it as a separator of complexes
  trans_table <- org.Hs.egUNIPROT %>% as.data.frame() %>%
    filter(uniprot_id == uniprot)
  # Get the entrez gene IDs that are mapped to a Uniprot ID
  l <- length(trans_table$gene_id)
  if(l==1) return(trans_table$gene_id)
  if(l==0){
    if(verbose) cat(uniprot," has no entrez gene id, NA has been returned!")
    if(!is.null(log_file)) write(x = c("", uniprot, "has no Entrez gene ID, NA has been returned!"), 
                                 file = log_file, append = TRUE, sep = "\t", ncolumns = 3)
    return(ifelse(na_action,NA, uniprot))
  }
  if(l>1){
    if(verbose)
      cat(uniprot, "has more than 1 Entrez gene ID", trans_table$gene_id, ", first one has been chosen! (Future enhancement to be confirmed: consider this as a family?)")
    if(!is.null(log_file)) 
      write(x = c("", uniprot, paste0(trans_table$gene_id, collapse = ","), "has more than 1 Entrez gene ID, first one has been chosen! (Future enhancement to be confirmed: consider this as a family?)"),
            file = log_file, append = TRUE, sep = "\t", ncolumns = 4)
    return(trans_table$gene_id[1])
  }
  stop()
}
getEntrezFromComposed <- function(sigID, sigMapper,sep=",", verbose=F, log_file = NULL, na_action = T){
  uniprots <- sigMapper %>% filter(SIG_ID == sigID) %>% select(COMPONENTS) %>% as.character()
  if(length(uniprots)==1 ){
    if(!grepl("SIGNOR-",uniprots)){
      new_genesList <- sapply(strsplit(x = uniprots, split = ",")[[1]], getEntrezFromUniprot, verbose=verbose, log_file = log_file, na_action=na_action) %>% paste(collapse = sep)
    }else{
      new_genesList <- sapply(strsplit(x = uniprots, split = ",")[[1]],function(n){
        if(grepl("SIGNOR-",n))
          return(getEntrezFromComposed(n,sigMapper=sigMapper,sep=sep, verbose, na_action=na_action))
        return(getEntrezFromUniprot(n,verbose,log_file, na_action))
      }) %>% paste(collapse = sep)
    }
  }
  if(length(uniprots)==0) stop("mixer sigmapper is needed, be causion with separators then!!")
  if(length(uniprots)>1) stop("Not a unique entry for signor IDS")
  return(new_genesList)
}

get_unreachable_out_nodes <- function(sif, verbose=F) {
  ig <- graph_from_data_frame(sif[, c("hi_ida", "hi_idb")], directed = TRUE)
  last_nodes<-V(ig)[degree(ig, mode = "out") == 0]
  first_nodes<-V(ig)[degree(ig, mode = "in") == 0]
  unreachable_out_nodes <- last_nodes[!sapply(last_nodes, function(last_node){
                                          any(first_nodes %in% igraph::dfs(ig, last_node, unreachable = F, mode= "in")$order)
                                        })]
  if(verbose){
    if (length(unreachable_out_nodes) > 0) {
      cat("The following out nodes are not reachable from any in nodes:", paste(unreachable_out_nodes$name, collapse = ", "), "\n")
    }
  }
  return(unreachable_out_nodes)
}
# Function to check if all in_nodes can reach all out_nodes
get_unreachable_out_nodes_loopBased <- function(sif, verbose=F) {
  ig <- graph_from_data_frame(sif[, c("hi_ida", "hi_idb")], directed = TRUE)
  last_nodes<-V(ig)[degree(ig, mode = "out") == 0]
  first_nodes<-V(ig)[degree(ig, mode = "in") == 0]
  
  unreachable_out_nodes <- c()
  
  for (last_node in last_nodes) {
    # Check if out_node is reachable from any of the in_nodes
    reachable <- sapply(first_nodes, function(in_node) {
      lengths(shortest_paths(ig, from = in_node, to = last_node)$vpath) > 0
    })
    
    if (!any(reachable)) {
      unreachable_out_nodes <- c(unreachable_out_nodes, V(ig)$name[last_node])
    }
  }
  if(verbose){
    if (length(unreachable_out_nodes) > 0) {
      cat("The following out nodes are not reachable from any in nodes:", paste(unreachable_out_nodes, collapse = ", "), "\n")
    }
  }
  return(unreachable_out_nodes)
}

## Adding layout to graph 
add_layout2att <-function(graph, verbose=F){
  ig <- graph_from_data_frame(graph$sif[, c("hi_ida", "hi_idb")], directed = TRUE)
  l <- layout_nicely(graph = ig,dim = 2)
  rownames(l) <- vertex_attr(ig,name = "name")
  graph$att <- graph$att %>% rowwise() %>% mutate(X = round(l[ID,1]+abs(min(l[,1])), digits = 0)*100,
                                                  Y = round(l[ID,2]+abs(min(l[,2])), digits = 0)*50)
  if(verbose) message("layout added for ", strsplit(x = graph$att$ID[1], split = "-")[[1]][2])
  return(graph)  
}
## write sif and att files
write_sif_att_files<- function(graph,spe, path_id,output_folder,verbose){
  file_pre <- paste0(spe,path_id)
  dir.create(output_folder, showWarnings = FALSE, recursive = T)
  write.table(graph$att, file = file.path(output_folder,paste0(file_pre,".att")), append = F, quote = F, sep = "\t", row.names = F, col.names = T)
  write.table(graph$sif, file = file.path(output_folder,paste0(file_pre,".sif")), append = F, quote = F, sep = "\t", row.names = F, col.names = F)
}
## MGI from SIF patch
mgi_from_sif_patched <- function(sif.folder, spe, verbose = F, entrez_symbol = NULL, dbannot = NULL){
  message("Loading graphs...")
  pgs <- hipathia:::load_graphs(sif.folder, spe, verbose = verbose)
  if (!is.null(dbannot) & !is.null(entrez_symbol)) {
    message("Adding functions to pathways...")
    pgs <- hipathia:::add_funs_to_pgs(pgs, entrez_symbol, dbannot, 
                           maxiter = 1000)
  }
  message("Creating MGI...")
  metaginfo <- create_metaginfo_object_patched(pgs, spe, by.user = TRUE)
  message("Created MGI with ", length(metaginfo$pathigraphs), 
          " pathway(s)")
  return(metaginfo)
}

create_metaginfo_object_patched <- function(fpgs, species, by.user = FALSE, basal.value = 0.5) {
  pathigraph.genes <- hipathia:::all_needed_genes(fpgs)
  labelids <- lapply(fpgs, function(pg) cbind(pg$label.id, 
                                              path.id = pg$path.id, path.name = pg$path.name))
  labelids <- do.call("rbind", labelids)
  rownames(labelids) <- labelids[, 1]
  genes.vals.05 <- matrix(basal.value, ncol = 2, nrow = length(pathigraph.genes), 
                          dimnames = list(pathigraph.genes, c("1", "2")))
  meta.05 <- NULL
  meta.05$pathigraphs <- fpgs
  meta.05$all.labelids <- labelids
  results.05 <- hipathia_patched(genes.vals.05, meta.05, test = FALSE, 
                         verbose = FALSE)
  results.dec.05 <- hipathia_patched(genes.vals.05, meta.05, decompose = TRUE, 
                             test = FALSE, verbose = FALSE)
  metaginfo <- NULL
  metaginfo$species <- species
  metaginfo$all.genes <- pathigraph.genes
  metaginfo$path.norm <- assay(results.dec.05, "paths")[, 
                                                        1]
  metaginfo$eff.norm <- assay(results.05, "paths")[, 1]
  metaginfo$pathigraphs <- fpgs
  metaginfo$all.labelids <- labelids
  metaginfo$group.by <- "pathways"
  metaginfo$by.user <- by.user
  return(metaginfo)
}
##################################### Patched Hipathia
hipathia_patched <- function (genes_vals, metaginfo, uni.terms = FALSE, GO.terms = FALSE, 
          custom.terms = NA, sel_assay = 1, decompose = FALSE, maxnum = 100, 
          verbose = TRUE, tol = 1e-06, test = TRUE) 
{
  # genes_vals <- hipathia:::normalize_data(genes_vals, by_quantiles = FALSE, 
  #                              by_gene = FALSE, percentil = FALSE)
  if (is(genes_vals, "SummarizedExperiment")) {
    coldata <- colData(genes_vals)
    genes_vals <- assay(genes_vals, sel_assay)
  }
  else {
    cols <- colnames(genes_vals)
    coldata <- data.frame(cols = cols, stringsAsFactors = FALSE)
  }
  if (test == TRUE) {
    if (is.null(genes_vals)) 
      stop("Missing input matrix")
    if (is.null(metaginfo)) 
      stop("Missing pathways object")
    hipathia:::test_matrix(genes_vals)
    hipathia:::test_pathways_object(metaginfo)
    hipathia:::test_tolerance(tol)
  }
  pathigraphs <- metaginfo$pathigraphs
  genes_vals <- hipathia:::add_missing_genes(genes_vals, genes = metaginfo$all.genes)
  results <- list()
  if (verbose == TRUE) 
    cat("Computing pathways...\n")
  results$by.path <- lapply(pathigraphs, function(pathigraph) {
    res <- list()
    suppressWarnings(res$nodes.vals <- hipathia:::nodes_values_from_genes(genes_vals, 
                                                               pathigraph$graph))
    if (decompose == FALSE) {
      respaths <- hipathia:::all_path_values(res$nodes.vals, pathigraph$effector.subgraphs, 
                                  maxnum = maxnum, tol = tol)
    }
    else {
      respaths <- hipathia:::all_path_values(res$nodes.vals, pathigraph$subgraphs, 
                                  maxnum = maxnum, tol = tol)
    }
    res$path.vals <- respaths[[1]]
    res$convergence <- respaths[[2]]
    return(res)
  })
  nodes <- do.call("rbind", lapply(results$by.path, function(x) x$nodes.vals))
  nodes_rd <- DataFrame(metaginfo$all.labelids[rownames(nodes), 
  ], node.name = hipathia:::get_node_names(metaginfo, rownames(nodes)), 
  node.type = hipathia:::get_node_type(metaginfo)$type, node.var = apply(nodes, 
                                                              1, var))
  nodes_se <- SummarizedExperiment(list(nodes = nodes), rowData = nodes_rd, 
                                   colData = coldata)
  paths <- do.call("rbind", lapply(results$by.path, function(x) x$path.vals))
  paths_rd <- DataFrame(path.ID = rownames(paths), path.name = hipathia:::get_path_names(metaginfo, 
                                                                              rownames(paths)), path.nodes = hipathia:::get_path_nodes(metaginfo, 
                                                                                                                            rownames(paths), decompose = decompose), decomposed = decompose)
  paths_se <- SummarizedExperiment(list(paths = paths), rowData = paths_rd, 
                                   colData = coldata)
  se_list <- list(nodes = nodes_se, paths = paths_se)
  if (uni.terms == TRUE) {
    if (verbose == TRUE) 
      cat("\nComputing Uniprot terms...\n")
    unis_se <- hipathia:::quantify_funs(paths_se, metaginfo, "uniprot")
    se_list$uni.terms <- unis_se
  }
  if (GO.terms == TRUE) {
    if (verbose == TRUE) 
      cat("\nComputing GO terms...\n")
    gos_se <- hipathia:::quantify_funs(paths_se, metaginfo, "GO")
    se_list$GO.terms <- gos_se
  }
  if (!is.na(custom.terms)) {
    if (verbose == TRUE) 
      cat("\nComputing custom terms...\n")
    custom_se <- hipathia:::quantify_funs(paths_se, metaginfo, dbannot)
    se_list$custom.terms <- custom_se
  }
  resmae <- MultiAssayExperiment(se_list)
  if (verbose == TRUE) 
    message("DONE")
  return(resmae)
}
### Functiona annotation 
# It not used!
## Checkin common complexes with 
get_annots<- function(signor_annot, spe, db="uniprot"){
  ## Here: entrez to symbol used for annotation has to be revised / check the 
  dbs <- c("uniprot", "GO")
  if(!db %in% dbs) stop("Not allowed db. allowed dbs are: ", toString(dbs))
  hipathia_annot<- hipathia:::load_annots(db = db, spe) %>% as_tibble()
  entrez_hgnc <- hipathia:::load_entrez_hgnc(species = spe)  %>% as_tibble()
#   ##  
#   complexes %>% filter(COMPLEX_NAME%in%entrez_hgnc$V2) %>% rowwise() %>% mutate(geneList= getEntrezFromComposed(SIG_ID,complexes,sep=",/,",verbose)) %>%
#     mutate(ENTREZ_in_Hi= entrez_hgnc$V1[entrez_hgnc$V2 == COMPLEX_NAME])
# ## protein family 
#   proteinFamilies %>% filter(PF_NAME %in% entrez_hgnc$V2) %>% rowwise() %>% mutate(geneList= getEntrezFromComposed(SIG_ID,proteinFamilies,sep=",",verbose)) %>%
#     mutate(ENTREZ_in_Hi= entrez_hgnc$V1[entrez_hgnc$V2 == PF_NAME])
# ## fusion protein 
#   fusionProteins  %>% filter(FP_NAME %in% entrez_hgnc$V2) %>% rowwise() %>% mutate(geneList= getEntrezFromComposed(SIG_ID,fusionProteins,sep=",",verbose)) %>%
#   mutate(ENTREZ_in_Hi= entrez_hgnc$V1[entrez_hgnc$V2 == FP_NAME])
  annotations<-list()
  annotations$signor_entrez_hgnc <- rbind(entrez_hgnc,
                              signor_annot %>% select(genesList,entitya) %>% rename( "V1" = genesList , "V2"=entitya)) %>% unique()
  annotations$signor_annot <- rbind(hipathia_annot,
                        signor_annot %>% select(entitya,entityb) %>% rename( "gene" =entitya, "function"=entityb)) %>% unique()
  return(annotations)
}

# this is an Hipathia example
getHipathia_report <- function(mgi, output_folder, verbose){
  port<- servr::random_port()
  fake_comp<- data.frame(row.names = names(mgi$eff.norm)) %>%
    mutate(name= hipathia::get_path_names(mgi, rownames(.)),
           "UP/DOWN"=0,
           statistic=0,
           p.value=1,
           FDRp.value=1)
  hipathia::create_report(comp = fake_comp, metaginfo = mgi, output_folder = "", path = output_folder, verbose = verbose)
  cat("Open a web browser and go to URL http://127.0.0.1:", port, "\n", sep = "")
  message("Press Ctrl + C to stop serving the report...\n")
  servr::httd(paste0(output_folder, "/pathway-viewer"), port = port, browser = TRUE, daemon = F)
}
### general functions
setdiff_all <- function(vec1, vec2) {
  return(setdiff(union(vec1, vec2), intersect(vec1, vec2)))
}
