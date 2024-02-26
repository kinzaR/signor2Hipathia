library("org.Hs.eg.db")
library(httr)
library(readr)
library(dplyr)

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
    data_df <- content(response, as = "text", encoding = "UTF-8") %>% read_tsv
    return(data_df)
  }
  stop("Failed to retrieve data. Check the API endpoint and parameters.")
}

getSifFromRow <- function(pathway_tsv, 
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
                          avoided_effect = c("form complex", "unknown")){
  pathway_tsv <- add_hi_effect(pathway_tsv)
  
  sif <- data.frame(source = pathway_tsv$entitya, source.type = pathway_tsv$typea,
                    sign = pathway_tsv$effect, 
                    target = pathway_tsv$entityb, target.type = pathway_tsv$typeb,
                    stringsAsFactors = F)
  sif$sign[sif$sign%in% activation_effect] <- "activation"
  sif$sign[sif$sign%in% inhibition_effect] <- "inhibition"
  if(length(sif[sif$sign %in% avoided_effect,"target"])>0 && !sif[sif$sign %in% avoided_effect,"target"] %in% sif$source){
    sif <- sif[!sif$sign %in% avoided_effect,]
  }
  ## annotation
  
  # change the node lable to the id
  # example for cpmplex N-hsa04520-67 72
  all_nodes <-unique(rbind(
    data.frame(uniprotID = pathway_tsv$ida,label=pathway_tsv$entitya, type=pathway_tsv$typea),
    data.frame(uniprotID = pathway_tsv$idb,label=pathway_tsv$entityb, type=pathway_tsv$typeb)))
  return(sif)
}
add_hi_effect<- function(pathway_tsv){
  
}
getEntrezFromName<- function(){
  x <- org.Hs.egUNIPROT
  # Get the entrez gene IDs that are mapped to a Uniprot ID
  mapped_genes <- mappedkeys(x)
  # Convert to a list
  xx <- as.list(x[mapped_genes])
  return(xx)
}

### general functions
setdiff_all <- function(vec1, vec2) {
  return(setdiff(union(vec1, vec2), intersect(vec1, vec2)))
}
