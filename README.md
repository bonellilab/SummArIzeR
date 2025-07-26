
<!-- README.md is generated from README.Rmd. Please edit that file -->

# SummArIzeR <a><img src="man/figures/SummArIzeRlogo.png" align="right" height="138"/></a>

<!-- badges: start -->
<!-- badges: end -->

`SummArIzeR` is an R package, that allows an easy use of EnrichR to
compare enrichment results from multiple databases of multiple
conditons. It results in a clustering of enriched terms and enables the
annotation of these terms by creating a promt for large language models
such as gpt4. Results can be vizualised in a Heatmap.

If you are using `SummArIzeR` for your publication, please cite us:

SummArIzeR: Simplifying cross-database enrichment result clustering and
annotation via large language models Marie Brinkmann, Michael Bonelli,
Anela Tosevska bioRxiv 2025.05.28.656331; doi:
<https://doi.org/10.1101/2025.05.28.656331>

## Features

- Perform enrichment analysis using `enrichR`.
- Allows analysis of multiple conditions
- Analyze up- and down-regulated genes separately.
- Filter terms by p-value and gene thresholds.
- Calculates similarities of results terms based on included genes
- Clusters terms using random walk algorithm
- Generates a prompt for a LLM to summarized cluster annotations
- Allows easy heatmap visualization

## Installation

You can install the development version of SummArIzeR from
[GitHub](https://github.com/) with:

``` r
# Install devtools if not already installed
install.packages("devtools")

# Install SummArIzeR
devtools::install_github("bonellilab/SummArIzeR")

# Additional packages

# Install factoextra from CRAN
install.packages("factoextra")

# Install `enrichR`
devtools::install_github("wjawaid/enrichR")
```

## Vignettes/Tutorials

Visit the SummArIzeR website for Info’s about:

- [getting started] (https://bonellilab.github.io/SummArIzeR/articles/01-getting-started.html)
- [a detailed workflow and customization] (https://bonellilab.github.io/SummArIzeR/articles/02-detailed_workflow.html)
- [how to prepare the input data from a DESeq2 object] (https://bonellilab.github.io/SummArIzeR/articles/03-dataframe-creation.html)

## Bug Reports

If you encounter any errors or issues, or if you have a suggestion
please file an issue
[here](https://github.com/bonellilab/SummArIzeR/issues).
