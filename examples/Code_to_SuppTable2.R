### Validation of SummArIzeR clustering and annptation
## dataset used: https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2022.868067/full#h11, figure 2

Exdata2_genes<-read.csv2("Example_data2.csv")
Exdata2_hall<-read.csv2("Example_data2_hallmarkgenes.csv")

Exdata2_genes$genes<-Exdata2_genes$id
Exdata2_genes<-Exdata2_genes %>% filter(genes %in% Exdata2_hall$A2M) # gene set used in the publication
Exdata2_genes$condition<-paste0("intercept_hall")

Termlist_all<-extractMultipleTerms(Exdata2_genes, condition_col = c("condition"), categories = c("GO_Biological_Process_2021","KEGG_2021_Human"), pval_threshold = 1 ,min_genes_threshold = 2, n = 500, split_by_reg = F)



Paper_terms<-c("Proteasome", 
               "Prion disease", 
               "Spinocerebellar ataxia", 
               "Parkinson disease", 
               "Glycolysis / Gluconeogenesis", 
               "Fructose and mannose metabolism", 
               "Mitophagy", 
               "HIF-1 signaling pathway", 
               "Central carbon metabolism in cancer", 
               "DNA replication", 
               "Renal cell carcinoma", 
               "Pyrimidine metabolism", 
               "p53 signaling pathway", 
               "Antigen processing and presentation", 
               "Gap junction", 
               "AGE-RAGE signaling pathway in diabetic complications", 
               "Ferroptosis", 
               "DNA replication", 
               "nucleotide-excision repair (GO:0006289)")

Paper_terms[Paper_terms %in% Termlist_all$Term == F]

Termlist_all<-Termlist_all %>% filter(Term %in% Paper_terms)
SummArIzeR::evaluateThreshold(Termlist_all)
plot<-TRUplotIgraph(Termlist_all, ts  = 0.32)
Genelist_test_cluster<-returnIgraphCluster(Termlist_all, ts = 0.32)

generateGPTPrompt(Genelist_test_cluster)
#Enter vector from LLM 
cluster_summary <- c(
  '1' = 'Exdata2 Metabolism and Hypoxia Signaling',
  '2' = 'Neurodegenerative and Protein Degradation Pathways',
  '3' = 'DNA Repair and Replication Processes',
  '4' = 'Carbohydrate Metabolism',
  '5' = 'Diabetic Complications and AGE-RAGE Signaling',
  '6' = 'Antigen Presentation and Immune Response',
  '7' = 'Regulated Cell Death (Ferroptosis)',
  '8' = 'Cell-Cell Communication (Gap Junctions)',
  '9' = 'Mitochondrial Quality Control (Mitophagy)',
  '10' = 'Nucleotide Metabolism (Pyrimidine Pathways)',
  '11' = 'Tumor Suppressor Signaling (p53 Pathway)'
)




Annotation_list<- annotateClusters(Genelist_test_cluster, cluster_summary)

##Similar clusters obseved, compared to publication


##comparison to word cloud annotation 
annotate_go_cluster <- function(
    go_terms,
    max_words = 5,
    exclude_words = c("via", "protein", "factor", "type", "specific")
) {
  
  # --- 1. Explicitly remove GO IDs (GO:0000000) ---
  remove_go_ids <- function(x) {
    gsub("GO:[0-9]+", " ", x)
  }
  
  # --- 2. Local word-count function (self-contained) ---
  count_words_local <- function(term) {
    
    docs <- tm::VCorpus(tm::VectorSource(term))
    docs <- tm::tm_map(docs, tm::content_transformer(remove_go_ids))
    docs <- tm::tm_map(docs, tm::content_transformer(tolower))
    docs <- tm::tm_map(docs, tm::removeNumbers)
    docs <- tm::tm_map(docs, tm::removeWords, tm::stopwords())
    docs <- tm::tm_map(docs, tm::removePunctuation)
    docs <- tm::tm_map(docs, tm::stripWhitespace)
    docs <- tm::tm_map(docs, tm::removeWords, exclude_words)
    
    tdm <- tm::TermDocumentMatrix(docs)
    v <- sort(slam::row_sums(tdm), decreasing = TRUE)
    
    data.frame(word = names(v), freq = v, stringsAsFactors = FALSE)
  }
  
  # --- 3. Run annotation ---
  df <- count_words_local(go_terms)
  if (nrow(df) == 0) return(NA_character_)
  
  df <- head(df[order(df$freq, decreasing = TRUE), ], max_words)
  
  paste(df$word, collapse = " / ")
}



for(i in unique(Annotationlist_Exdata2paper$Cluster_Annotation)){
  go_terms<-unlist(Annotationlist_Exdata2paper$unique_terms_per_cluster[Annotationlist_Exdata2paper$Cluster_Annotation == i])
  
  print(i)
  print(annotate_go_cluster(go_terms))
  
}

### SummArIzeR full potential 


Exdata2_genes<-read.csv2("Exdata2_paper.csv")
Exdata2_hall<-read.csv2("cancer_paper_hallmarkgenes.csv")

Exdata2_genes$genes<-Exdata2_genes$id
Exdata2_genes$condition[Exdata2_genes$genes %in% Exdata2_hall$A2M]<-paste0("intercept_hall")
Exdata2_genes$condition[Exdata2_genes$genes %in% Exdata2_hall$A2M == F]<-paste0("no_intercept_hall")


Termlist_all<-extractMultipleTerms(Exdata2_genes, condition_col = c("condition"), categories = c("GO_Biological_Process_2021","KEGG_2021_Human"), pval_threshold = 0.05 ,min_genes_threshold = 3, n = 20, split_by_reg = F)

SummArIzeR::evaluateThreshold(Termlist_all)
plot<-TRUplotIgraph(Termlist_all, ts  = 0.19)
Genelist_test_cluster<-returnIgraphCluster(Termlist_all, ts = 0.19)

generateGPTPrompt(Genelist_test_cluster)
#Enter vector from LLM 
cluster_summary <- c(
  '1' = 'Exdata2 Metabolism and Hypoxia Signaling',
  '2' = 'Lipid Metabolism, Immunity, and Autophagy',
  '3' = 'DNA Damage Response and Repair Mechanisms',
  '4' = 'Cell Cycle Regulation and Neurodegeneration Pathways',
  '5' = 'Mitotic Spindle and Chromosome Organization',
  '6' = 'Innate Immune Signaling via dsRNA and STAT Activation',
  '7' = 'Retinoid and Eicosanoid Metabolism',
  '8' = 'DNA Replication and Nucleotide Excision Repair',
  '9' = 'Ferroptosis-Mediated Cell Death',
  '10' = 'Cell–Cell Communication via Gap Junctions',
  '11' = 'p53 Tumor Suppressor Signaling Pathway'
)



Annotation_list<- annotateClusters(Genelist_test_cluster, cluster_summary)
SummArIzeR::plotHeatmap(Annotation_list, column_names_centered = F, cluster_rows = T, cluster_columns = T)
SummArIzeR::plotBubbleplot(Annotation_list, plot_colors = list(up = c("#ececec", "#ffb9b1", "#F87060"), down = c("#ececec",
                                                                                                                 "#3172CC", "#102542"),  default = c("#ececec", "#3172CC", "#102542")))

