library(hipathia)
mgi<-readRDS("results_FA_all4web_cleaned2/2024_23_may_jue_10_21Vhi_2.10.0_mgi.RDS")
res <-hipathia(genes_vals = normalize_data(hipathia::brca), metaginfo = mgi, verbose = T)

load("results_FA_all4web_cleaned2/job_data(5).rdata")# from web
load("/home/kinza/Downloads/job_data(7).rdata")# from web
#############################
# Pathway visualization
##############################
output_folder<-"tmp/results_FA_4web_cleaned/web"
dir.create(output_folder)
if(difexp == T){
  # Define node colors
  cat("\nSpecies is : ", pathways$species,"********\n")
  colors_de <- node_color_per_de(results, pathways, sample_group, cond1,cond2)
  #colors_de_hipathia <- node_color_per_de(results, pathways, sample_group,"Tumor", "Normal", colors = "hipathia")
  create_report(wt, pathways, path=output_folder,output_folder="", node_colors = colors_de, verbose = T)
  # 4: 
  this_comp<-wt
  pathway<- "hsaSIGXOR0AML_TF"
  metaginfo<-pathways
  hipathia:::create_node_and_edge_attributes(this_comp, pathway, metaginfo)
}else{
  if(difexp == "difexpuni"){
    # Node colors with Uniprot grouping
    #colors_uni <- node_color_per_de(results, pathways, sample_group, "Tumor", "Normal", group_by = "uniprot")
    #create_report(wt, pathways, path=output_folder,output_folder="", node_colors = colors_uni)
    
    colors_de_uni <- node_color_per_de(results, pathways, sample_group, "Tumor","Normal", group_by = "uniprot")
    create_report(wt, pathways, path=output_folder,output_folder="",node_colors = colors_de_uni, group_by = "uniprot")
    
  }else{
    if(difexp == "difexpgp"){
      # Node colors with genes grouping
      colors_gen  <- node_color_per_de(results, pathways, sample_group, "Tumor", "Normal", group_by = "genes")
      create_report(wt, pathways, path=output_folder,output_folder="", node_colors = colors_gen, group_by = "genes")
    }else{
      if(difexp == "difexpgo"){
        # Node colors with genes grouping
        colors_de_go  <- node_color_per_de(results, pathways, sample_group, "Tumor", "Normal", group_by = "GO")
        create_report(wt, pathways, path=output_folder,output_folder="", node_colors = colors_de_go, group_by = "GO")
      }
    }
  }
}
status("95")

################## REPORT

cat("Creating report...\n")
switch(species,
       hsa = {speciesTitle <- "Human (Homo sapiens)"},
       mmu = {speciesTitle <- "Mouse (Mus musculus)"},
       rno = {speciesTitle <- "Rat (Rattus norvegicus)"})

results <- init.report(analysis)

results <- add.section(results,"Input parameters",0)

if(dataset=="gex"){
  results <- add.param(results,"Expression file",basename(exp_file),1)
}
results <- add.param(results,"Design file",basename(design_file),1)
if( design_type != "continuous"){
  results <- add.param(results,"Comparison",paste0(cond1," vs ",cond2),1)
  
  results <- add.param(results,"Paired analysis",ifelse(paired, "Yes", "No"),1)
} else {
  results <- add.param(results,"Comparison","Correlation with a continuous variable",1)
}
#results <- add.param(results,"Design file",paste0(basename(design_file)))

results <- add.param(results,"Species",speciesTitle,1)
#results <- add.param(results,"Decomposed paths",decompose)

results <- add.section(results,"Pathways",0)
results <- add.pathwayViewer(results, "session", "pathways", 1)
results <- add.section(results,"Circuit values",0)
results <- add.download(results,"Circuit values","paths_vals.txt",1)
if( design_type != "continuous"){
  if(file.exists(paste0(output_folder,"/paths_heatmap.png"))){
    results <- add.image(results,"Heatmap","paths_heatmap.png",1)
  }
  if(file.exists(paste0(output_folder,"/paths_pca.png"))){
    results <- add.image(results,"PCA","paths_pca.png",1)
  }
}
results <- add.table(results,"Pathways summary","pathways_summary_table.txt",1, page_size = 18, "pathways_summary_vals.txt")
results <- add.table(results,"Circuit significance","paths_significance.txt",1)
if(go==T | uniprot==T){
  results <- add.section(results,"Function based analysis",0)
  if(go==T){
    results <- add.download(results,"GO terms values","go_vals.txt",1)
    if( design_type != "continuous"){
      if(file.exists(paste0(output_folder,"/go_heatmap.png"))){
        results <- add.image(results,"Heatmap","go_heatmap.png",1)
      }
      if(file.exists(paste0(output_folder,"/go_pca.png"))){
        results <- add.image(results,"PCA","go_pca.png",1)
      }
    }
    results <- add.table(results,"GO term significance","go_significance.txt",1)
    results <- add.table(results,"Highest significant GO ancestor (HSGOA)","HSGOA.txt",1)
  }
  if(uniprot==T){
    results <- add.download(results,"Uniprot keywords values","uniprot_vals.txt",1)
    if( design_type != "continuous"){
      #paste0(output_folder,"/go","_heatmap.png")
      if(file.exists(paste0(output_folder,"/uniprot_heatmap.png"))){
        results <- add.image(results,"Heatmap","uniprot_heatmap.png",1)            
      }
      if(file.exists(paste0(output_folder,"/uniprot_pca.png"))){
        results <- add.image(results,"PCA","uniprot_pca.png",1)            
      }
    }
    results <- add.table(results,"Uniprot keyword significance","uniprot_significance.txt",1)
  }
}
write(render.xml(results),file=paste0(output_folder,"/report.xml"))
#unlink(paste0(output_folder,"/sifs4CellMaps"),recursive = T)

# if(exists("wt")){
#   wt$status <- wt$"UP/DOWN"
#   wt$has_changed <- T
#   path_json <- create.path.info(wt,metaginfo$pathigraphs)
#   write(path_json,file=paste0(output_folder,"/pathways/path_info.json"))
# }
cp_command2 <- paste0("cp  -r ",output_folder,"/pathway-viewer/pathways ",output_folder)
system(cp_command2)

cat("[Finished]\n")

status("100")
