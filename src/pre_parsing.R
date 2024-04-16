# install.packages("VennDiagram")
library(VennDiagram)
library(readr)
library(dplyr)
source("src/utils.R")
pathways<-list()
pathways$pathways_from_selectors <- read.table("issues/pathway_names.tsv", sep = "\t") %>% .$V1
paths <-getPathwayData()
pathways$pathways_from_info <- unique(paths$paths_info$sig_id)
pathways$pathways_from_relations <- unique(paths$paths_relations$pathway_id)
# Create the Venn diagram
catNames <-  c("From selector [Web]", "From info [API]", "From relations tab [API]  ")
catNamesDesc <-  c("Pathways from selector [Web]", "Pathways from info [API]", "Pathways from relations tab [API]")
cols <- c("#FFC14D", "#52A677", "#759EB2")
venn.plot <- VennDiagram::venn.diagram(x = pathways,
                       category.names = catNames, 
                       filename = NULL, fill=cols, output = TRUE, category.legend = TRUE ,
                       label.fontface = "bold",  # Make the category names bold
                       cat.cex = 1.9,  # Adjust the size of the category names
                       cat.col = cols
                       )
# Add labels for each category
venn.plot[[7]]$label  <- paste(
  c(paste0("(",venn.plot[[7]]$label, ")"),
  setdiff(pathways$pathways_from_selectors, 
          c(pathways$pathways_from_info,pathways$pathways_from_relations))), collapse="\n")
# venn.plot[[8]]$label <- paste(
#   c(paste0("(",venn.plot[[8]]$label, ")"),
#   setdiff(intersect(pathways$pathways_from_selectors,pathways$pathways_from_info), pathways$pathways_from_relations)), collapse="\n")

# venn.plot[[9]]$label <- paste(
#   c(paste0("(", venn.plot[[9]]$label,")"),
#   intersect(pathways$pathways_from_info, setdiff(pathways$pathways_from_relations ,pathways$pathways_from_selectors))), collapse="\n")
venn.plot[[9]]$label <- paste(
  c(paste0("(", venn.plot[[9]]$label,")"),
    setdiff(pathways$pathways_from_info,
            c(pathways$pathways_from_relations,pathways$pathways_from_selectors))), collapse="\n")

venn.plot[[11]]$label <- paste(
  c( paste0("(", venn.plot[[11]]$label,")") )) #,
  # intersect(pathways$pathways_from_info, intersect(pathways$pathways_from_selectors ,pathways$pathways_from_relations))), collapse="\n")

venn.plot[[12]]$label <- paste(
 c(paste0("(",venn.plot[[12]]$label, ")") ))#,
  # setdiff(intersect(pathways$pathways_from_info,pathways$pathways_from_relations),pathways$pathways_from_selectors)), collapse="\n")
venn.plot[[13]]$label <- paste(
  c(paste0("(",venn.plot[[13]]$label, ")"),
    setdiff(pathways$pathways_from_relations, c(pathways$pathways_from_selectors, pathways$pathways_from_info))), collapse="\n")
# Display the Venn diagram
grid.newpage()
grid.draw(venn.plot)

lg <- legendGrob(labels=catNamesDesc, 
                 pch=rep(19,length(catNames)),
                 gp=gpar(col=cols, fill="gray"), byrow = F)

# gridExtra::grid.arrange(lg,venn.plot, nrow =2)
gridExtra::grid.arrange(venn.plot,lg, nrow=2, heights=c(20,5))


#####

library(gplots)
v.table <- venn(pathways)
print(v.table)
print(venn.plot)
