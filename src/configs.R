activation_effect <- c("up-regulates activity",
                       "up-regulates quantity by expression",
                       "up-regulates",
                       "up-regulates quantity",
                       "up-regulates quantity by stabilization")
inhibition_effect <- c("down-regulates",
                       "down-regulates quantity by repression",
                       "down-regulates quantity",
                       "down-regulates activity",
                       "down-regulates quantity by destabilization")
avoided_effect <- c("form complex", "unknown")

#### APIs
# pathways
pathways_relations <- "getPathwayData.php?relations"
pathways_description <- "getPathwayData.php?description"
# nodes info
proteinFamilies <- "getDataInternal.php?proteinFamilies=all"
fusionProteins <- "getDataInternal.php?fusionProteins"
complexes <- "getDataInternal.php?complexes=all"

