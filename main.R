#!/usr/bin/env /opt/R/4.1.2/bin/Rscript
#  #!/usr/bin/env /opt/R/4.1.2/bin/Rscript
## # !/usr/bin/env /opt/R/4.3.1/bin/Rscript <- new version  is imcompatible with web !
# Load needed lirbrary
print(version)
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(hipathia))
# suppressPackageStartupMessages(library(dplyr))
# Source my utils 
source("src/utils.R")
## here I have declared constants
source("src/configs.R")
# Define command-line options
option_list <- list(
  make_option(c("-s","--spe"), type = "character", default = "hsa",
              help = "Species variable. Allowed choices: 'hsa', 'mmu', 'rno'. (default: hsa)"),
  make_option(c("-v", "--verbose"), action = "store_true", default = T,#FALSE 
              help = "Enable verbose mode."),
  make_option(c("-r", "--readySifs"), action = "store_true", default = FALSE, 
              help = "Read from already created sif and att files."),
  make_option(c("-p", "--pathways_list"), type = "character", default = "SIGNOR-AML,SIGNOR-LBC,SIGNOR-PDAP,SIGNOR-AAAM",# "all",#"SIGNOR-AML,SIGNOR-LBC,SIGNOR-PDAP"
              help = "Vector of the IDs of the pathways to parse fromSignor. By default, all available pathways are loaded. Example: 'SIGNOR-AML,SIGNOR-LBC,SIGNOR-PDAP'."),
  make_option("--score", type = "double", default = 0.1,
              help = "The minimum significance score allowed, the range is from 0.1 to 1. Signor set 0.1 as the minimum score as 0 stands for no evidence of interaction.[ for more information : https://www.google.com/url?sa=t&source=web&rct=j&opi=89978449&url=https://signor.uniroma2.it/documentation/SIGNOR_3_score_Documentation_final.docx&ved=2ahUKEwim0JGkkI2GAxX1YPEDHRT8BKIQFnoECBoQAQ&usg=AOvVaw2y_b2VjYMFJgoA3BilRe95]
              By default 0.1"),
  make_option(c("-o","--output_folder"), type = "character", default = "tmp_fa",
              help = "Output folder")
)
# Parse command-line arguments
opt <- parse_args(OptionParser(option_list = option_list))

if(!opt$readySifs){# get data from API
  pathwayData <- getSignorData(api_url = pathways_relations_path)
  # remove special chars in node names
  pathwayData <-pathwayData %>% mutate(entitya= gsub("'|\"", "", entitya),
                         entityb= gsub("'|\"", "", entityb))
  # Species filter with tax_id
  pathwayData <- pathwayData %>% filter(tax_id == switch (opt$spe,"hsa" = "9606","mmu"= "10090","rno" ="10116"),
                                        score >= opt$score)
  # filter if pathway list avail
  print(opt$pathways_list)
  if(!is.null(opt$pathways_list) && opt$pathways_list != "all"){ #4.3.0>:  will throw an error with && and || instead of simply evaluating the first element if the length of the vectors being compared > 1.
    paths_ids <- strsplit(opt$pathways_list, ",")[[1]]
    if(any(paths_ids %in% pathwayData$pathway_id)){
      pathwayData <- pathwayData %>% filter(pathway_id %in% paths_ids)
      if(!all(paths_ids %in% pathwayData$pathway_id))
        warning(paste(setdiff(paths_ids, pathwayData$pathway_id), collapse = ", "), " not found using Signor API!")
      }else
        stop("No intersection between the specified pathway list and the pathways available in the Signor database.")
    }
  if(any(!pathwayData$effect %in% c(activation_effect,inhibition_effect,avoided_effect))){
    paths_new_effect<- pathwayData %>% filter(!effect %in% c(activation_effect,inhibition_effect,avoided_effect)) %>%
      select(pathway_id, pathway_name) %>% unique
    cat("New effect were find, please check these pathways: ")
    print(paths_new_effect)
    cat(".\nThey have these not included interaction(s): ",
         setdiff(pathwayData$effect, c(activation_effect,inhibition_effect,avoided_effect)))
    stop()
  }
  if(any(duplicated(gsub(pattern = pattern4ids, replacement = "0", x = unique(pathwayData$pathway_id)))))
  if(any(duplicated(gsub(pattern = "N|P", replacement = "X", x = gsub(pattern=pattern4ids,
                                                                      replacement="0",
                                                                      x=unique(pathwayData$pathway_id))))))
    stop("re-check IDs of pathways!")
  #remove N from id because it's provoking problem in the functiona nnotation while is replacing N in P! is a bug in hipathia package
  pathwayData <-pathwayData %>% mutate(pathway_id = gsub(pattern = "N|P", replacement = "X", x = gsub(pattern=pattern4ids,
                                                                                                     replacement="0",
                                                                                                     x=pathway_id)))
  # save pathnames and ids in a separate var to save the pathway_names tsv file later
  path_names <- pathwayData %>% select(pathway_id, pathway_name) %>% unique()
  # select only relevant columns
  pathwayData <- pathwayData %>% select(pathway_id,
                                        entitya, typea, ida, databasea,
                                        entityb, typeb, idb, databaseb,
                                        effect) %>% unique
  # curation : has to be moved to a proper function withclear pipeline for graph curation!
  pathwayData <- pathwayData %>% filter(ida!=idb)
  # Here I have to detect cyclic sub-paths!
  pathways_tsv <- split(pathwayData, pathwayData$pathway_id)
  # here a function for all using lapply
  # read annot files from signor APIs before
  # NOTE: here I have to add the geneslist instead of calculating if evry time, this will reduce /optimize time 
  proteinFamilies <- getSignorData(proteinFamilies_path) %>% rowwise() %>%
    mutate(COMPONENTS=stringr::str_replace_all(COMPONENTS,pattern = c(",or,"=",","and"="/")))# some correction in the original mappers
  fusionProteins <- getSignorData(fusionProteins_path) %>% rename(SIG_ID=FP_SIG_ID)%>% rowwise() %>%
    mutate(COMPONENTS=stringr::str_replace_all(COMPONENTS,pattern = c(",or,"=",","and"="/")))# some correction in the original mappers
  complexes <- getSignorData(complexes_path) %>% rowwise() %>%
    mutate(COMPONENTS=stringr::str_replace_all(COMPONENTS,pattern = c(",or,"=",","and"="/")))# some correction in the original mappers
  
  graphs <- lapply(pathways_tsv, write_Sif_Att_FromRow,
         spe=opt$spe,
         activation_effect = activation_effect,
         inhibition_effect = inhibition_effect,
         avoided_effect = avoided_effect,
         proteinFamilies = proteinFamilies,
         fusionProteins = fusionProteins,
         complexes = complexes,
         output_folder=opt$output_folder, verbose=opt$verbose)
  ## Here check only parced pathways with a good id
  not_parsed_pathways <- names(graphs[sapply(graphs, is.null)])
  if(any(c("SIGNOR0AML","SIGNOR0LBC","SIGNOR0PDAP") %in% not_parsed_pathways))
    warning("Maca's pathways are not all parsed!")
  write.table(x = not_parsed_pathways, file = file.path(opt$output_folder,"not_parsed.tsv"), append = F, quote = F, sep = "\t", row.names = F, col.names = F)
  graphs[sapply(graphs, is.null)] <- NULL
  path_names <- path_names %>% filter(pathway_id %in% names(graphs))
    
  write.table(path_names, file = file.path(opt$output_folder,paste0("name.pathways_",opt$spe,".txt")), append = F, quote = F, sep = "\t", col.names = F,row.names = F)
  ### Here I have to write annotations 
  all_pheno_as_out <- do.call(rbind,lapply(graphs, function(g){g$pheno_as_out})) %>% as_tibble()
  # Clean NA genesList
  all_pheno_as_out<-all_pheno_as_out %>% filter(!is.na(genesList))
  all_stimulus_as_out <- do.call(rbind,lapply(graphs, function(g){g$stimulus_as_out})) %>% as_tibble()
  write.table(all_pheno_as_out, file = file.path(opt$output_folder,paste0("phenotypes_annot_",opt$spe,".tsv")), append = F, quote = F,  sep = "\t", col.names = T,row.names = F)
  write.table(all_stimulus_as_out, file = file.path(opt$output_folder,paste0("stimuli_as_out_",opt$spe,".tsv")), append = F, quote = F,  sep = "\t", col.names = T,row.names = F)
  message("Sifs and Atts ready for parsing to mgi!")
}
### from sif/att files to MGIs
cat("hipathia version is : ",package.version("hipathia"))
if(!exists("all_pheno_as_out"))
  all_pheno_as_out <- read_table(file = file.path(opt$output_folder,paste0("phenotypes_annot_",opt$spe,".tsv")), show_col_types = F)
annotations <- get_annots(signor_annot = all_pheno_as_out, spe=opt$spe, db = "uniprot")
if(package.version("hipathia")=="2.10.0"){
  mgi <- hipathia::mgi_from_sif(sif.folder = opt$output_folder, spe = opt$spe, entrez_symbol = as.data.frame(annotations$signor_entrez_hgnc), dbannot = as.data.frame(annotations$signor_annot))
}else{
  mgi <- mgi_from_sif_patched(sif.folder = opt$output_folder, spe = opt$spe, verbose=opt$verbose, entrez_symbol = as.data.frame(annotations$signor_entrez_hgnc), dbannot = as.data.frame(annotations$signor_annot))
}
now <-format(Sys.time(), "%Y_%d_%b_%a_%H_%M")
saveRDS(object = mgi, file = file.path(opt$output_folder,paste0(now,"Vhi_",package.version("hipathia"),"_mgi.RDS")))
message("DONE! Results are in ",opt$output_folder, " folder.")
### to json 
library(jsonlite)
path_j_list<-lapply(mgi$pathigraphs, function(p){
  c(p$path.name)
})
names(path_j_list)<-gsub("hsa", "", names(path_j_list), ignore.case = F)
capsule <- list()
capsule$'SIGNOR Pathways' <- list()
capsule$`SIGNOR Pathways`$'from apis pathways' <- path_j_list
paths_j <- toJSON(capsule, pretty = T, auto_unbox = T)
write(paths_j, file = file.path(opt$output_folder,paste0(now,"Vhi_",package.version("hipathia"),opt$spe,"_pathway_list_signor")))
