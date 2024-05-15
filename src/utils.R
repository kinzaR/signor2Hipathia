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
    data_df <- content(response, as = "text", encoding = "UTF-8") %>% read_tsv(show_col_types = F)
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
                          proteinFamilies_path,
                          fusionProteins_path,
                          complexes_path,
                          output_folder=NULL,
                          verbose=F){
  path_id <- pathway_tsv$pathway_id %>% unique %>% gsub(pattern = "-", replacement = "", x = .)
  if(verbose) cat("Parsing the ",pathway_tsv$pathway_id %>% unique, "...")
  pathway_tsv <- add_hi_effect(pathway_tsv)# Change the interaction to hipathia standard 
  proteinFamilies <- getSignorData(proteinFamilies_path)%>% rowwise() %>%
    mutate(COMPONENTS=stringr::str_replace_all(COMPONENTS,pattern = c(",or,"=",","and"="/")))# some correction in the original mappers
  fusionProteins <- getSignorData(fusionProteins_path) %>% rename(SIG_ID=FP_SIG_ID)%>% rowwise() %>%
    mutate(COMPONENTS=stringr::str_replace_all(COMPONENTS,pattern = c(",or,"=",","and"="/")))# some correction in the original mappers
  complexes <- getSignorData(complexes_path) %>% rowwise() %>%
    mutate(COMPONENTS=stringr::str_replace_all(COMPONENTS,pattern = c(",or,"=",","and"="/")))# some correction in the original mappers
  
  graph<-list()
  graph$att <- get_hi_att(pathway_tsv,path_id, proteinFamilies, fusionProteins, complexes,spe, verbose) # get the att file of hipathia
  graph$sif <- get_hi_sif(pathway_tsv, graph$att ,verbose)
  if(!is.null(output_folder)) write_sif_att_files(graph, spe, path_id, output_folder,verbose)
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
resolve_signor_ids <- function(id, type, proteinFamilies, fusionProteins, complexes, verbose=F){
  new_id <- switch(type,
                   "protein" = getEntrezFromUniprot(id,verbose),
                   "proteinfamily" = getEntrezFromComposed(id,proteinFamilies,sep=",",verbose),
                   "fusion protein" = getEntrezFromComposed(id,fusionProteins,sep=",",verbose), 
                   "complex" = getEntrezFromComposed(id,complexes,sep=",/,",verbose),
                   "smallmolecule" = NA, # has to be changes fopr metabopathia
                   "mirna" = NA,
                   "chemical" = NA,# has to be changes fopr metabopathia
                   "phenotype" = NA
                   )
  return(new_id)
}
get_hi_att <- function(pathway_tsv, path_id, proteinFamilies, fusionProteins, complexes,spe="hsa", verbose=F){
  if(verbose) print("Getting Att...")
  all_nodes <- tibble(
    node_id = c(pathway_tsv$ida,pathway_tsv$idb),
    node_label = c(pathway_tsv$entitya, pathway_tsv$entityb),
    node_type = c(pathway_tsv$typea, pathway_tsv$typeb)
  ) %>% unique()
  # genesList= resolve_signor_ids(node_id, node_type, proteinFamilies, fusionProteins, complexes, verbose=F)
  all_nodes<-all_nodes %>% 
    rowwise() %>% mutate(genesList=resolve_signor_ids(node_id, node_type, proteinFamilies, fusionProteins, complexes, verbose=verbose)) %>%
    mutate(id_indicator= case_when(is.na(genesList) ~ node_id,
                                   .default = genesList))
  hi_ids <- sapply(all_nodes$id_indicator, function(g){
    if(grepl(",/,", x = g)) strsplit(x = g, split = ",/,")
    else g
  }) %>% unlist() %>% unique %>% na.omit %>% as_tibble(.) %>%
    mutate(ids=row.names(.))
  all_nodes<- all_nodes %>% rowwise() %>%mutate(hi_id = ifelse(id_indicator %in% hi_ids$value, hi_ids$ids[hi_ids$value==id_indicator],
                                                                hi_ids%>% filter(value %in% (strsplit(x = id_indicator, split = ",/,")[[1]])) %>% select(ids) %>% pull() %>% paste0(collapse = " ")
                                                                )) %>%
    mutate(hi_id= paste0("N-",spe,path_id,"-",hi_id))
  att <- tibble(ID = all_nodes$hi_id,
                label	=all_nodes$node_label,
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
                tooltip=NA) %>% unique()
  return(att)
}
get_hi_sif <- function(pathway_tsv, att,verbose=F){
  if(verbose) print("Getting Sif...")
  sif<-pathway_tsv %>% rowwise() %>% mutate(hi_ida = att$ID[att$label==entitya]) %>% 
    mutate(hi_idb=att$ID[att$label==entityb]) %>% select(hi_ida,sign, hi_idb) %>% unique()
  return(sif)
}
getEntrezFromUniprot<- function(uniprot, verbose=F){
  if(uniprot=="/") return("/") # is a exception to keep it as a separator of complexes
  trans_table <- org.Hs.egUNIPROT %>% as.data.frame() %>%
    filter(uniprot_id == uniprot)
  # Get the entrez gene IDs that are mapped to a Uniprot ID
  l <- length(trans_table$gene_id)
  if(l==1) return(trans_table$gene_id)
  if(l==0){
    if(verbose) cat(uniprot," has no entrez gene id, NA has been returned!")
    return(NA)
    }
  if(l>1){
    if(verbose) cat(uniprot," has more than 1 entrez gene id ",trans_table$gene_id,", first one has been chosen!")
    return(trans_table$gene_id[1])
  }
  stop()
}
getEntrezFromComposed <- function(sigID, sigMapper,sep=",", verbose=F){
  uniprots <- sigMapper %>% filter(SIG_ID == sigID) %>% select(COMPONENTS) %>% as.character()
  if(length(uniprots)==1 ){
    if(!grepl("SIGNOR-",uniprots)){
      new_genesList <- sapply(strsplit(x = uniprots, split = ",")[[1]], getEntrezFromUniprot, verbose=verbose) %>% paste(collapse = sep)
    }else{
      new_genesList <- sapply(strsplit(x = uniprots, split = ",")[[1]],function(n){
        if(grepl("SIGNOR-",n))
          return(getEntrezFromComposed(n,sigMapper=sigMapper,sep=sep, verbose=verbose))
        return(getEntrezFromUniprot(n,verbose))
      }) %>% paste(collapse = sep)
    }
  }
  if(length(uniprots)==0) stop("mixer sigmapper is needed, be causion with separators then!!")
  if(length(uniprots)>1) stop("Not a unique entry fro signor IDS")
  return(new_genesList)
}

## write sif and att files
write_sif_att_files<- function(graph,spe, path_id,output_folder,verbose){
  file_pre <- paste0(spe,path_id)
  dir.create(output_folder, showWarnings = FALSE)
  write.table(graph$att, file = file.path(output_folder,paste0(file_pre,".att")), append = F, quote = F, sep = "\t", row.names = F, col.names = T)
  write.table(graph$sif, file = file.path(output_folder,paste0(file_pre,".sif")), append = F, quote = F, sep = "\t", row.names = F, col.names = F)
}
## MGI from SIF patch
mgi_from_sif_patched <- function(sif.folder, spe, verbose = F, entrez_symbol = NULL, dbannot = NULL){
  message("Loading graphs...")
  pgs <- hipathia:::load_graphs(sif.folder, spe, verbose = verbose)
  if (!is.null(dbannot) & !is.null(entrez_symbol)) {
    message("Adding functions to pathways...")
    pgs <- add_funs_to_pgs(pgs, entrez_symbol, dbannot, 
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


### general functions
setdiff_all <- function(vec1, vec2) {
  return(setdiff(union(vec1, vec2), intersect(vec1, vec2)))
}
