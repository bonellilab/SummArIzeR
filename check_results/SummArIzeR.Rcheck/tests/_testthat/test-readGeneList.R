library(mockery)

test_that("readGeneList performs enrichment and filters results correctly", {
  # Mock input
  listpathSig <- data.frame(
    genes = c("Gene1", "Gene2", "Gene3", "Gene4"),
    log2fold = c(2, -3, 0.5, -0.2)
  )
  
  # Create a mock function for enrichr
  enrichr_mock <- mock(
    list(
      `GO_Biological_Process` = data.frame(
        Term = c("Process1", "Process2"),
        Genes = c("Gene1;Gene2", "Gene3"),
        Adjusted.P.value = c(0.01, 0.2)
      )
    )
  )
  
  # Replace enrichr with the mock
  stub(readGeneList, "enrichr", enrichr_mock)
  
  # Run the function and check results
  result <- readGeneList(
    listpathSig,
    name = "Test",
    category = "GO_Biological_Process",
    split_by_reg = FALSE,
    pval_threshold = 0.05,
    min_genes_threshold = 1
  )
  
  # Assertions
  expect_true(is.data.frame(result))
  expect_named(result, c("Term", "Genes", "adj_pval", "dbs", "condition"))
  expect_equal(unique(result$Term), "Process1")
})

