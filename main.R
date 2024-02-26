# Load needed lirbrary
library(magrittr)

# Source my utils 
source("src/utils.R")
## here I have declared constants
source("src/configs.R")
pathwayData <- getSignorData(api_url = pathways_relations)
pathways_tsv <- split(pathwayData, pathwayData$pathway_id)
pathway_tsv <- pathways_tsv$`SIGNOR-AAAM`

sif <- getSifFromRow(pathway_tsv, activation_effect, inhibition_effect, avoided_effect)
# att <- 