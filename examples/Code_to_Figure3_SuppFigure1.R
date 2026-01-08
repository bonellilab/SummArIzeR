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
#####Figure 2C####
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

#Interferon signaling

List_a<-unique(Termlist_all$Genes[Termlist_all$Term == "cellular response to type I interferon (GO:0071357)"])
Termlist <- unique(Genelist_test_cluster$Term[Genelist_test_cluster$Cluster == 4])
Termlist <- Termlist[Termlist != "cellular response to type I interferon (GO:0071357)"]
List_b<-unique(Termlist_all$Genes[Termlist_all$Term %in% Termlist])

#VennDiagram

x <- list(Type_I_Int= List_a, rest_cluster = List_b)

ggVennDiagram(x)

##### Supplementary Figure 2C####
##Adding a similarity score to the manual clustering

library(dplyr)
library(tidyr)
library(purrr)
library(tibble)

Termlist_all<-extractMultipleTerms(cytokine_df, condition_col = c("Comparison_short"), categories = c("GO_Biological_Process_2021"), pval_threshold = 0.05,min_genes_threshold = 5, n = 5, split_by_reg = F)

unique(Termlist_all$Term)

manual_cluster_1<-c("eosinophil chemotaxis (GO:0048245)" ,"eosinophil migration (GO:0072677)" ,"cellular response to chemokine (GO:1990869)","positive regulation of cell migration (GO:0030335)" ,"regulation of cell migration (GO:0030334)"  )
manual_cluster_2<-c("cellular response to cytokine stimulus (GO:0071345)"  ,"cytokine-mediated signaling pathway (GO:0019221)")
manual_cluster_3<-c("external encapsulating structure organization (GO:0045229)" ,"extracellular matrix organization (GO:0030198)" ,"extracellular structure organization (GO:0043062)"  ,"glycosaminoglycan biosynthetic process (GO:0006024)")
manual_cluster_3<-c("external encapsulating structure organization (GO:0045229)" ,"extracellular matrix organization (GO:0030198)" ,"extracellular structure organization (GO:0043062)"  ,"glycosaminoglycan biosynthetic process (GO:0006024)")
manual_cluster_4<-c("inflammatory response (GO:0006954)", "response to interleukin-1 (GO:0070555)" ,"response to tumor necrosis factor (GO:0034612)")
manual_cluster_5<-c("cellular response to interferon-gamma (GO:0071346)" , "type I interferon signaling pathway (GO:0060337)","cellular response to type I interferon (GO:0071357)")

manual_cluster_6<-unique(Termlist_all$Term[Termlist_all$Term %in% c(manual_cluster_1, manual_cluster_2, manual_cluster_3, manual_cluster_4, manual_cluster_5) == F])
manual_clusters <- list(
  cluster_1 = manual_cluster_1,
  cluster_2 = manual_cluster_2,
  cluster_3 = manual_cluster_3,
  cluster_4 = manual_cluster_4,
  cluster_5 = manual_cluster_5,
  cluster_6 = manual_cluster_6
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
df_term_man<-df_term_mapped
#### coparison to SummArIzeR 

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


df_term_mapped <- Annotation_list %>%
  select(Cluster_Annotation, unique_terms_per_cluster) %>%
  unnest(unique_terms_per_cluster) %>%
  rename(Term = unique_terms_per_cluster) %>%
  distinct() %>% dplyr::rename("terms" = "Term", "cluster" = "Cluster_Annotation")


df_term_mapped$jaccard<-mapply(jaccardfun, df_term_mapped$terms, df_term_mapped$cluster)
df_term_mapped_summarizer<-df_term_mapped

### Include revigo results 
rvigo<-read.delim("ToolComparisonsDataFigures/rrvgoResults.tsv")
head(rvigo)

#prepare a table with unique terms and clusters
rvigo$term<-paste0(rvigo$term, " (", rvigo$go, ")")
##change some name missmatches (ids similar)
rvigo$term[rvigo$term == "cellular response to type II interferon (GO:0071346)"]<-"cellular response to interferon-gamma (GO:0071346)"
rvigo$term[rvigo$term == "type I interferon-mediated signaling pathway (GO:0060337)"]<-"type I interferon signaling pathway (GO:0060337)"

df_term_mapped<- rvigo %>% dplyr::select(term, parentTerm) %>% dplyr::rename("terms" = "term" , "cluster" = "parentTerm")
df_term_mapped$jaccard<-mapply(jaccardfun, df_term_mapped$terms, df_term_mapped$cluster)
df_term_mapped_revigo<-df_term_mapped

### Include simplify enrichment results 

senrich<-read.delim("ToolComparisonsDataFigures/simplifyEnrichmentResults.tsv")
senrich$cluster<-as.character(senrich$cluster)
head(senrich)
senrich$term<-paste0(senrich$term, " (", senrich$go_id, ")")
#prepare a table with unique terms and clusters
df_term_mapped <- senrich %>% dplyr::select(term, cluster) %>% dplyr::rename("terms" = "term")
df_term_mapped$jaccard<-mapply(jaccardfun, df_term_mapped$terms, df_term_mapped$cluster)
df_term_mapped_senrich<-df_term_mapped

##Heatmap with the jaccards 
head(df_term_man)
head(df_term_mapped_summarizer)
head(df_term_mapped_revigo)
head(df_term_mapped_senrich)


# Prepare each dataframe
df_man <- df_term_man %>%
  select(terms, cluster, jaccard) %>%
  rename(
    cluster_man  = cluster,
    jaccard_man  = jaccard
  )

df_summarizer <- df_term_mapped_summarizer %>%
  select(terms, cluster, jaccard) %>%
  rename(
    cluster_summarizer = cluster,
    jaccard_summarizer = jaccard
  )

df_revigo <- df_term_mapped_revigo %>%
  select(terms, cluster, jaccard) %>%
  rename(
    cluster_revigo = cluster,
    jaccard_revigo = jaccard
  )

df_senrich <- df_term_mapped_senrich %>%
  select(terms, cluster, jaccard) %>%
  rename(
    cluster_senrich = cluster,
    jaccard_senrich = jaccard
  )

# Join all by terms
df_all <- list(
  df_man,
  df_summarizer,
  df_revigo,
  df_senrich
) %>%
  reduce(full_join, by = "terms")

df_all<-df_all %>% filter(jaccard_summarizer != 0 & is.na(jaccard_summarizer) == F)
unique(df_all$cluster_summarizer)

manual_order <- c(
  "Extracellular Matrix and Structural Organization",
  "Type I Interferon-Mediated Antiviral Response",
  "Regulation of Cell Migration and MAPK Signaling",
  "Inflammatory and Chemotactic Response",
  "Cytokine Signaling Response"
  
)

df_all <- df_all %>%
  mutate(
    cluster_summarizer = factor(
      cluster_summarizer,
      levels = manual_order
    )
  ) 
colnames(df_all)
mat<-as.matrix(df_all %>% select(jaccard_summarizer, jaccard_revigo,jaccard_senrich,jaccard_man)) 

rownames(mat)<-df_all$terms

row_anno <- rowAnnotation(
  cluster_summarizer = df_all$cluster_summarizer,
  cluster_man = df_all$cluster_man,
  cluster_revigo = df_all$cluster_revigo,
  cluster_senrich = df_all$cluster_senrich,
  annotation_name_side = "top"
)

library(ComplexHeatmap)
library(circlize)  # for colorRamp2

# Example: define manual color scale
my_colors <- c("low" = "white", "mid" = "#BDD7FB", "high" = "#00366C")

# Define breakpoints for the scale
breaks <- c(min(mat, na.rm = TRUE), median(mat, na.rm = TRUE), max(mat, na.rm = TRUE))

Heatmap(
  mat,
  row_split = df_all$cluster_summarizer,
  left_annotation = row_anno,
  col = colorRamp2(
    breaks,
    my_colors,
  ),
  cluster_rows = FALSE,
  cluster_columns = FALSE
)
####Code to Figure 2D####

### Barchart with mean jaccard values

library(ggplot2)
library(ggthemes)
# calculate means excluding zeros
mean_values <- data.frame(
  method = c("Summarizer", "Manual", "REVIGO", "simplifyEnrichment"),
  mean_jaccard = c(
    mean(df_all$jaccard_summarizer[df_all$jaccard_summarizer != 0], na.rm = TRUE),
    mean(df_all$jaccard_man[df_all$jaccard_man != 0], na.rm = TRUE),
    mean(df_all$jaccard_revigo[df_all$jaccard_revigo != 0], na.rm = TRUE),
    mean(df_all$jaccard_senrich[df_all$jaccard_senrich != 0], na.rm = TRUE)
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



####Code to Figure 2E and supplementary table 1####

#Enrichment analysis
Termlist_all<-extractMultipleTerms(cytokine_df, condition_col = c("Comparison_short"), categories = c("GO_Biological_Process_2025", "BioPlanet_2019"), pval_threshold = 0.05, n = 5, split_by_reg = T)
#Clustering
Genelist_test_cluster<-returnIgraphCluster(Termlist_all, ts = 0.05)
#Generating promt
generateGPTPrompt(Genelist_test_cluster)
#Integrating LLM annotation
## Summary Terms with ChatGPT 4-turbo
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
##Summary Terms with DeepSeek

cluster_summary <- c(
  '1' = 'Skeletal and Muscular System Development',
  '2' = 'Inflammatory and Immune Signaling Pathways',
  '3' = 'Mitotic Cell Cycle and Chromosome Segregation',
  '4' = 'Antiviral Defense and Interferon Signaling',
  '5' = 'Calcium Ion Regulation and Cardiac Muscle Contraction',
  '6' = 'Extracellular Matrix Organization',
  '7' = 'Ephrin Receptor Signaling',
  '8' = 'p38 MAPK Signaling Pathway'
)

##Summary Terms with Claude

cluster_summary <- c(
  '1' = 'Skeletal and Cartilage Development',
  '2' = 'Cytokine Signaling and Inflammatory Response',
  '3' = 'Mitotic Cell Cycle Regulation',
  '4' = 'Antiviral Immune Response',
  '5' = 'Calcium Ion Regulation',
  '6' = 'Extracellular Matrix Organization',
  '7' = 'Ephrin Signaling',
  '8' = 'p38 MAPK Signaling'
)

##Summary Terms with PerplexityAI

#Perplexity AI
cluster_summary <- c(
  '1' = 'Skeletal and Musculoskeletal System Development',
  '2' = 'Immune Response and Intracellular Signaling Regulation',
  '3' = 'Mitotic Cell Cycle and Spindle Assembly Regulation',
  '4' = 'Antiviral Defense and Interferon-Mediated Immune Signaling',
  '5' = 'Regulation of Calcium Signaling and Cardiac Muscle Contraction',
  '6' = 'Extracellular Matrix Organization and Collagen Pathways',
  '7' = 'Ephrin Receptor Signaling Pathway',
  '8' = 'p38 MAPK Signaling Pathway'
)

#Gemini AI

cluster_summary <- c('1' = 'Skeletal and Muscle Development', 
                     '2' = 'Immune and Inflammatory Responses and Signaling', 
                     '3' = 'Mitotic Cell Cycle and Spindle Organization', 
                     '4' = 'Antiviral Defense and Interferon Signaling', 
                     '5' = 'Regulation of Calcium Ion Release in Muscle Contraction', 
                     '6' = 'Extracellular Matrix Organization', 
                     '7' = 'Ephrin Receptor B Signaling', 
                     '8' = 'p38 MAPK Signaling')
print(cluster_summary)
Annotation_list<- annotateClusters(Genelist_test_cluster, cluster_summary)
#Vizualisation
plotHeatmap(Annotation_list, annotation_bar = F, column_names_centered = F)

