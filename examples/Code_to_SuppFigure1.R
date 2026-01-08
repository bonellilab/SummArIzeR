### Validation of SummArIzeR clustering 
## dataset used: https://www.nature.com/articles/s41467-021-22544-y, figure 3g

library(dplyr)
library(ggplot2)
library(SummArIzeR)
library(enrichR)
library(ComplexHeatmap)
library(circlize)  

Epidata<-read.csv2("Epig_paper.csv")

#result as data.frame
Epidata$condition<-paste0("Enhancer genes")
Epidata$genes<-Epidata$Enhancer_genes

Termlist_all<-extractMultipleTerms(Epidata, condition_col = c("condition"), categories = c("KEGG_2019_Human"), pval_threshold = 1 ,min_genes_threshold = 2, n = 500, split_by_reg = F)

#filter terms from the paper 

paper_terms <- c(
  "Hippo signaling pathway",
  "Wnt signaling pathway",
  "Cushing syndrome",
  "Pancreatic cancer",
  "Chronic myeloid leukemia",
  "Colorectal cancer",
  "Basal cell carcinoma",
  "Signaling pathways regulating pluripotency of stem cells",
  "Endometrial cancer",
  "Hepatocellular carcinoma",
  "Gastric cancer",
  "Breast cancer",
  "Acute myeloid leukemia",
  "Nucleotide excision repair",
  "Mismatch repair",
  "Central carbon metabolism in cancer",
  "ErbB signaling pathway",
  "Renal cell carcinoma",
  "Glycerophospholipid metabolism",
  "Arginine biosynthesis",
  "Axon guidance",
  "Small cell lung cancer",
  "Human papillomavirus infection",
  "PI3K-Akt signaling pathway",
  "Focal adhesion",
  "ECM-receptor interaction",
  "Rap1 signaling pathway",
  "Ras signaling pathway",
  "MAPK signaling pathway"
  
)

#use only terms from the publication
Termlist_all<-Termlist_all %>% filter(Term %in% paper_terms)
Genelist_test_cluster<-returnIgraphCluster(Termlist_all, ts = 0.23)
SummArIzeR::evaluateThreshold(Termlist_all)
plot<-TRUplotIgraph(Termlist_all, ts  = 0.23)
generateGPTPrompt(Genelist_test_cluster)
#Enter vector from LLM 
cluster_summary <- c(
  '1' = 'Oncogenic signaling and cancer metabolism',
  '2' = 'Cell–ECM adhesion and survival signaling',
  '3' = 'Ras/MAPK signal transduction',
  '4' = 'Developmental and cancer-associated signaling pathways',
  '5' = 'DNA damage repair pathways',
  '6' = 'Amino acid biosynthesis',
  '7' = 'Neuronal guidance signaling',
  '8' = 'Membrane lipid metabolism',
  '9' = 'Renal cancer–associated pathways'
)



Annotation_list<- annotateClusters(Genelist_test_cluster, cluster_summary)
SummArIzeR::plotHeatmap(Annotation_list, column_names_centered = F, cluster_rows = T, cluster_columns = T)
SummArIzeR::plotBubbleplot(Annotation_list)



## Calculate JaccardIndex 

##manual clusters 

cl1<-paper_terms[1:13]
cl2<-paper_terms[14:15]
cl3<-paper_terms[16]
cl4<-paper_terms[17:18]
cl5<-paper_terms[19]
cl6<-paper_terms[20]
cl7<-paper_terms[21]
cl8<-paper_terms[22:26]
cl9<-paper_terms[27:29]

manual_clusters <- list(
  cluster_1 = cl1,
  cluster_2 = cl2,
  cluster_3 = cl3,
  cluster_4 = cl4,
  cluster_5 = cl5,
  cluster_6 = cl6,
  cluster_7 = cl7,
  cluster_8 = cl8,
  cluster_9 = cl9
)


jaccard <- function(a, b) {
  length(intersect(a, b)) / length(union(a, b))
}


term2genes <- Termlist_all %>%
  group_by(Term) %>%
  summarise(genes = list(unique(Genes)), .groups = "drop") %>%
  deframe()

#
#prepare a table with unique terms

df_term<-as.data.frame(matrix(nrow = length(unique(Termlist_all$Term))))
df_term$terms<-unique(Termlist_all$Term)
cluster_lookup <- manual_clusters %>%
  imap_dfr(~ tibble(
    cluster = .y,
    terms   = .x
  )) %>%
  unnest(terms)

df_term_mapped <- df_term %>%
  left_join(
    cluster_lookup,
    by = c("terms" = "terms")
  )

jaccardfun<-function(Termi, clusteri){
  
  a<-unique(Termlist_all$Genes[Termlist_all$Term == Termi])
  Termlist_b<-unique(df_term_mapped$terms[df_term_mapped$cluster == clusteri &df_term_mapped$terms != Termi])
  b<-unique(Termlist_all$Genes[Termlist_all$Term %in% Termlist_b])
  jac<-jaccard(a,b)
  return(jac)
}

df_term_mapped$jaccard<-mapply(jaccardfun, df_term_mapped$terms, df_term_mapped$cluster)

##compare to SummArIzeR

term_gene_long <- Termlist_all %>%
  select(Term, Genes) %>%
  distinct() %>%
  rename(Gene = Genes)
cluster_term_long <- Annotation_list %>%
  select(Cluster_Annotation, unique_terms_per_cluster) %>%
  unnest(unique_terms_per_cluster) %>%
  rename(Term = unique_terms_per_cluster) %>%
  distinct()


cluster_term_gene_long <- cluster_term_long %>%
  inner_join(term_gene_long, by = "Term")


jaccardfun<-function(Termi, clusteri){
  
  a<-unique(Termlist_all$Genes[Termlist_all$Term == Termi])
  Termlist_b<-unique(cluster_term_long$Term[cluster_term_long$Cluster_Annotation == clusteri &cluster_term_long$Term != Termi])
  b<-unique(Termlist_all$Genes[Termlist_all$Term %in% Termlist_b])
  jac<-jaccard(a,b)
  return(jac)
}



cluster_term_long$jaccard<-mapply(jaccardfun, cluster_term_long$Term, cluster_term_long$Cluster_Annotation)
df_term_mapped$Term<-df_term_mapped$terms

cluster_joined <- df_term_mapped %>%
 left_join(cluster_term_long,
    by = "Term",
    suffix = c("_manual", "_auto")
  ) %>% filter(jaccard_auto>0)


mat<-as.matrix(cluster_joined %>% select(jaccard_manual, jaccard_auto))
rownames(mat)<-cluster_joined$Term

row_anno <- rowAnnotation(
  Cluster_Annotation = cluster_joined$Cluster_Annotation,
  cluster = cluster_joined$cluster,
  annotation_name_side = "top"
)



# Example: define manual color scale
my_colors <- c("low" = "white", "mid" = "#BDD7FB", "high" = "#00366C")

# Define breakpoints for the scale
breaks <- c(min(mat, na.rm = TRUE), median(mat, na.rm = TRUE), max(mat, na.rm = TRUE))

Heatmap(
  mat,
  row_split = cluster_joined$Cluster_Annotation,
  left_annotation = row_anno,
  col = colorRamp2(
    breaks,
    my_colors,
  ),
  cluster_rows = FALSE,
  cluster_columns = FALSE
)

### Barchart with mean jaccard values

library(ggplot2)
library(ggthemes)
# calculate means excluding zeros
mean_values <- data.frame(
  method = c("Cytoscape", "Summarizer"),
  mean_jaccard = c(
    mean(cluster_joined$jaccard_manual[cluster_joined$jaccard_manual != 0], na.rm = TRUE),
    mean(cluster_joined$jaccard_auto[cluster_joined$jaccard_auto != 0], na.rm = TRUE)
  )
)
mean_values$method <- factor(
  mean_values$method,
  levels = mean_values$method[order(mean_values$mean_jaccard, decreasing = TRUE)]
)

# bar chart
ggplot(mean_values, aes(x = method, y = mean_jaccard, fill = method)) +
  geom_col(width = 0.6) +
  labs(
    x = "Method",
    y = "Mean Jaccard index (non-zero)",
    title = "Comparison of Mean Jaccard Indices"
  ) + scale_fill_manual(values = c("#6E458D", "#33658A", "#C8C8C8", "#86BBD8")) +
  theme_base()


