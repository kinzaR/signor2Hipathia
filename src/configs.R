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
pathways_relations_path <- "getPathwayData.php?relations"
pathways_description_path <- "getPathwayData.php?description"
# nodes info
proteinFamilies_path <- "getDataInternal.php?proteinFamilies=all"
fusionProteins_path <- "getDataInternal.php?fusionProteins"
complexes_path <- "getDataInternal.php?complexes=all"

# separators between components
easy_sep <- "," #OR
proteinFamilies_sep <- "," #OR

complexes_sep <- ",/," #AND
fusionProteins_sep <- ",/," #AND



