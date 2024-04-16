## get data
source("src/configs.R")
source("src/utils.R")
complexes <- getSignorData(complexes_path) %>% .[,c(1:3)] %>% 
  rename_("SIG_ID"= 1, "LABEL"=2, "COMPONENTS"=3) %>%
  mutate(TYPE = "complexe")
fusionProteins <- getSignorData(fusionProteins_path) %>%.[,c(1:3)] %>%
  rename("SIG_ID"= 1, "LABEL"=2, "COMPONENTS"=3) %>%
  mutate(TYPE = "fusionProtein")
proteinFamilies <- getSignorData(proteinFamilies_path) %>%.[,c(1:3)] %>%
  rename("SIG_ID"= 1, "LABEL"=2, "COMPONENTS"=3) %>%
  mutate(TYPE = "proteinFamilie")
all <- rbind(complexes,fusionProteins, proteinFamilies) %>% as.data.frame()
rownames(all) <- all$SIG_ID
## remove internal ids
# patterns "SIGNOR-C" for complexes, "SIGNOR-FP" for Fusion proteins, and "SIGNOR-PF" for protein family
types <- data.frame(pattern =  c("SIGNOR-C", "SIGNOR-FP", "SIGNOR-PF"),
                    types = c( "complexe", "fusionProtein", "proteinFamily"))
all$node_type<-""
test<-sapply(complexes$SIG_ID, function(x){
  components<- all[x,"COMPONENTS"]
  if(grepl("SIGNOR-", components)){
    print(components)
    res <- strsplit(components, ",")[[1]]
    
    zz<-sapply(res, function(r){
      print(r)
      type <- "uniprot"
      # check the pattern
      if(grepl("SIGNOR-", r)) {
        type <- types[grep(paste(types$pattern, collapse="|"), r),2]
        new_component <- all[r,"COMPONENTS"]
        all[x,]$node_type<-paste(all[x,]$node_type, type, sep = ",")
        return(new_component)
      }
      all[x,]$node_type<-paste(all[x,]$node_type, type, sep = ",")
      return(r)
    })
  }
})
result <- lapply(complexes$COMPONENTS, function(vec) strsplit(vec, ",")[[1]])
