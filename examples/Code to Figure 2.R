#Code to Figure 2
####Libraries####
library(SummArIzeR)
library(factoextra)
library(enrichR)
library(ggVennDiagram)
#####Load data####
data("example_data")

#####Figure 2B####

Termlist_all<-extractMultipleTerms(cytokine_df, condition_col = c("Comparison_short"), categories = c("GO_Biological_Process_2021"), pval_threshold = 0.05,min_genes_threshold = 5, n = 5, split_by_reg = F)
Genelist_test_cluster<-returnIgraphCluster(Termlist_all, ts = 0.2)
plot<-TRUplotIgraph(Termlist_all, ts  = 0.2)
generateGPTPrompt(Genelist_test_cluster)
#Enter vector from LLM 
cluster_summary <- c(
  '1' = 'Cytokine Signaling Response',
  '2' = 'Regulation of Cell Migration and MAPK Signaling',
  '3' = 'Inflammatory and Chemotactic Response',
  '4' = 'Type I Interferon-Mediated Antiviral Response',
  '5' = 'Extracellular Matrix and Structural Organization',
  '6' = 'Cellular Response to Oxidative Stress',
  '7' = 'Glycosaminoglycan Biosynthesis',
  '8' = 'Kinetochore Assembly and Function',
  '9' = 'Positive Regulation of Cellular Processes'
)


Annotation_list<- annotateClusters(Genelist_test_cluster, cluster_summary)

#Analysis of clustering differences
#Comparing genes enriched in GO:terms and old/new cluster
#Defense response to virus
List_a<-unique(Termlist_all$Genes[Termlist_all$Term == "defense response to virus (GO:0051607)"])
Termlist <- unique(Genelist_test_cluster$Term[Genelist_test_cluster$Cluster == 4])
Termlist <- Termlist[Termlist != "defense response to virus (GO:0051607)"]
List_b<-unique(Termlist_all$Genes[Termlist_all$Term %in% Termlist])

#VennDiagram

x <- list(defense_to_virus= List_a, rest_cluster = List_b)

ggVennDiagram(x)

#glycosaminoglycan biosynthetc process
List_a<-unique(Termlist_all$Genes[Termlist_all$Term == "glycosaminoglycan biosynthetic process (GO:0006024)"])
Termlist <- unique(Genelist_test_cluster$Term[Genelist_test_cluster$Cluster == 5])
List_b<-unique(Termlist_all$Genes[Termlist_all$Term %in% Termlist])

#VennDiagram

x <- list(Glyco= List_a, rest_cluster = List_b)

ggVennDiagram(x)

#Interferone signaling

List_a<-unique(Termlist_all$Genes[Termlist_all$Term == "cellular response to type I interferon (GO:0071357)"])
Termlist <- unique(Genelist_test_cluster$Term[Genelist_test_cluster$Cluster == 4])
Termlist <- Termlist[Termlist != "cellular response to type I interferon (GO:0071357)"]
List_b<-unique(Termlist_all$Genes[Termlist_all$Term %in% Termlist])

#VennDiagram

x <- list(Type_I_Int= List_a, rest_cluster = List_b)

ggVennDiagram(x)



####Code to Figure 2D####

#Enrichment analysis
Termlist_all<-extractMultipleTerms(cytokine_df, condition_col = c("Comparison_short"), categories = c("GO_Biological_Process_2025", "BioPlanet_2019"), pval_threshold = 0.05, n = 5, split_by_reg = T)
#Clustering
Genelist_test_cluster<-returnIgraphCluster(Termlist_all, ts = 0.05)
#Generating promt
generateGPTPrompt(Genelist_test_cluster)
#Integrating LLM annotation
cluster_summary <- c(
  '1' = 'Skeletal and Cartilage Development',
  '2' = 'Integrated Cellular Signaling and Immune Response',
  '3' = 'Cell Cycle Regulation and Chromosome Segregation',
  '4' = 'Antiviral Defense and Interferon Signaling',
  '5' = 'Calcium Regulation in Cardiac Function',
  '6' = 'Extracellular Matrix Organization',
  '7' = 'Ephrin Receptor B Signaling',
  '8' = 'p38 MAPK Signaling Pathway'
)
print(cluster_summary)
Annotation_list<- annotateClusters(Genelist_test_cluster, cluster_summary)
#Vizualisation
plotHeatmap(Annotation_list, annotation_bar = F, column_names_centered = F)


