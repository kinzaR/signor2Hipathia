#!/usr/bin/env /opt/R/4.3.1/bin/Rscript
# Load needed lirbrary
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(hipathia))
# Source my utils 
source("src/utils.R")
## here I have declared constants
source("src/configs.R")
# Define command-line options
option_list <- list(
  make_option(c("-s","--spe"), type = "character", default = "hsa",
              help = "Species variable. Allowed choices: 'hsa', 'mmu', 'rno'. (default: hsa)"),
  make_option(c("-v", "--verbose"), action = "store_true", default = FALSE, 
              help = "Enable verbose mode."),
  make_option(c("-p", "--pathways_list"), type = "character", default = "SIGNOR-AML,SIGNOR-LBC,SIGNOR-PDAP",
              help = "Vector of the IDs of the pathways to parse fromSignor. By default, all available pathways are loaded. Example: 'SIGNOR-AML,SIGNOR-LBC,SIGNOR-PDAP'."),
  make_option("--score", type = "double", default = 0.1,
              help = "The minimum significance score allowed, the range is from 0.1 to 1. Signor set 0.1 as the minimum score as 0 stands for no evidence of interaction.[ for more information : https://www.google.com/url?sa=t&source=web&rct=j&opi=89978449&url=https://signor.uniroma2.it/documentation/SIGNOR_3_score_Documentation_final.docx&ved=2ahUKEwim0JGkkI2GAxX1YPEDHRT8BKIQFnoECBoQAQ&usg=AOvVaw2y_b2VjYMFJgoA3BilRe95]
              By default 0.1"),
  make_option(c("-o","--output_folder"), type = "character", default = "tmp",
              help = "Output folder")
)
# Parse command-line arguments
opt <- parse_args(OptionParser(option_list = option_list))

# get data from API
pathwayData <- getSignorData(api_url = pathways_relations_path)
# Species filter with tax_id
pathwayData <- pathwayData %>% filter(tax_id == switch (opt$spe,"hsa" = "9606","mmu"= "10090","rno" ="10116"),
                                      score >= opt$score)
# filter if pathway list avail 
if(!is.null(opt$pathways_list)) pathwayData <- pathwayData %>% filter(pathway_id %in% strsplit(opt$pathways_list, ",")[[1]])
if(any(!pathwayData$effect %in% c(activation_effect,inhibition_effect,avoided_effect))){
  paths_new_effect<- pathwayData %>% filter(!effect %in% c(activation_effect,inhibition_effect,avoided_effect)) %>% 
    select(pathway_id, pathway_name) %>% unique
  cat("New effect were find, please check these pathways: ")
  print(paths_new_effect)
  cat(".\nThey have these not included interaction(s): ",
       setdiff(pathwayData$effect, c(activation_effect,inhibition_effect,avoided_effect)))
  stop()
}

# save pathnames and ids in a separate var to save the pathway_names tsv file later
path_names <- pathwayData %>% select(pathway_id, pathway_name) %>% unique()
# select only relevant columns 
pathwayData <- pathwayData %>% select(pathway_id,
                                      entitya, typea, ida, databasea,
                                      entityb, typeb, idb, databaseb,
                                      effect) %>% unique
pathways_tsv <- split(pathwayData, pathwayData$pathway_id)
# here a function for all using lapply

graphs <- lapply(pathways_tsv, write_Sif_Att_FromRow, 
       spe=opt$spe, 
       activation_effect = activation_effect, 
       inhibition_effect = inhibition_effect, 
       avoided_effect = avoided_effect,
       proteinFamilies_path = proteinFamilies_path,
       fusionProteins_path = fusionProteins_path,
       complexes_path = complexes_path,
       output_folder=opt$output_folder, verbose=opt$verbose)
path_names <- path_names %>% mutate(pathway_id= gsub(pattern = "-", replacement = "", x = pathway_id))
write.table(path_names, file = file.path(opt$output_folder,paste0("name.pathways_",opt$spe,".txt")), append = F, quote = F, sep = "\t", col.names = F,row.names = F)
message("Sifs and Atts ready for parsing to mgi!")
### from sif/att files to MGIs
mgi <- mgi_from_sif_patched(sif.folder = opt$output_folder, spe = opt$spe, verbose=opt$verbose,entrez_symbol = NULL, dbannot = NULL )
saveRDS(object = mgi, file = file.path(opt$output_folder,"mgi.RDS"))
message("DONE! Results are in ",opt$output_folder, " folder.")
 