#hipathia function without normalization
library(SummarizedExperiment)
library(MultiAssayExperiment)
local_hipathia <- function(genes_vals, metaginfo, uni.terms = FALSE, GO.terms = FALSE, 
          custom.terms = NA, sel_assay = 1, decompose = FALSE, maxnum = 100, 
          verbose = TRUE, tol = 1e-06, test = TRUE) 
{
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

local_create_metaginfo_object <- function (fpgs, species, by.user = FALSE, basal.value = 0.5){
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
  
  results.05 <- local_hipathia(genes.vals.05, meta.05, test = FALSE,
                               verbose = T)
  results.dec.05 <- local_hipathia(genes.vals.05, meta.05, decompose = TRUE, 
                                   test = FALSE, verbose = T)
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

### for easy graphs 
creatGraph <- function(sif,att){
  rownames(att) <- att[, 1]
  ig <- igraph::graph.data.frame(sif[, c(1, 3)], directed = TRUE)
  E(ig)$relation <- unlist(sapply(sif[, 2], function(x) {
    ifelse(x == "activation", 1, -1)
  }))
  E(ig)$curved <- FALSE
  V(ig)$nodeX <- 0
  V(ig)$nodeY <- 0
  V(ig)$shape <- "rectangle"
  V(ig)$label.cex <- 45
  V(ig)$label.color <- "white"
  V(ig)$width <- 15
  V(ig)$height <- 15
  V(ig)$tooltip <- ""
  glist <- sapply(as.character(att[, "genesList"]), function(x) {
    unlist(strsplit(x, split = ","))
  })
  names(glist) <- att[, "ID"]
  V(ig)$genesList <- glist[V(ig)$name]
  if (!is.null(att$label)) 
    V(ig)$label <- att[V(ig)$name, "label"]
  return(ig)
}


completGraph<-function(s2p,pathway,pathigraphs){
  pathigraphs[[pathway]]$graph <- s2p
  subs <- hipathia:::create_subgraphs(pathigraphs[[pathway]]$graph)
  pathigraphs[[pathway]]$subgraphs <- subs[[1]]
  pathigraphs[[pathway]]$subgraphs.mean.length <- subs[[2]]
  ces <- hipathia:::create_effector_subgraphs(pathigraphs[[pathway]])
  pathigraphs[[pathway]]$effector.subgraphs <- ces
  pathigraphs[[pathway]]$path.name <- pathway
  pathigraphs[[pathway]]$path.id <- pathway
  labs <- cbind(V(pathigraphs[[pathway]]$graph)$name, 
                V(pathigraphs[[pathway]]$graph)$label)
  pathigraphs[[pathway]]$label.id <- labs
  colnames(pathigraphs[[pathway]]$label.id) <- c("name", 
                                                 "label")
  return(pathigraphs)
}
