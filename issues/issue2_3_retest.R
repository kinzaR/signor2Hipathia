## get data
source("src/configs.R")
source("src/utils.R")
complexes <- getSignorData(complexes_path) %>% .[,c(1:3)] %>% 
  rename_("SIG_ID"= 1, "LABEL"=2, "COMPONENTS"=3) %>%
  mutate(TYPE = "complexe")
complexes_pub <- getSignorData("API/getComplexData.php")
setdiff(complexes$SIG_ID,complexes_pub$SIGNOR_ID)
setdiff(complexes_pub$SIGNOR_ID,complexes$SIG_ID)
length(intersect(complexes_pub$SIGNOR_ID,complexes$SIG_ID)) == length(complexes$SIG_ID)

compare<-tibble(id=complexes$SIG_ID, 
                internalComp=complexes$COMPONENTS, 
                publComp=complexes_pub$MEMBERS[complexes_pub$SIGNOR_ID == complexes$SIG_ID]) 
compare$publComp <- gsub(";", ",", compare$publComp)
compare[compare$internalComp != compare$publComp,]
View(compare[compare$id == "SIGNOR-C246",])

