06_analysis_CARNIVAL_results
================
Rosa Hernansaiz-Ballesteros
14/08/2020

### License Info

This program is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your
option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
Public License for more details.

Please check <http://www.gnu.org/licenses/>.

## Introduction

This is the final part of our series of transcriptomics tutorials. In
previous parts, we normalised RNASeq data for differential analysis, and
used the differential analysis results in **PROGENy** and **DOROTHea**.
In the previous tutorial, we demonstrated **CARNIVAL** with TF
activities from **DOROTHea**, a prior knowledge network from
**Omnipath**, and weights from **PROGENy**.

## CARNIVAL output

CARNIVAL (CAusal Reasoning for Network identification using Integer
VALue programming) returns a list with 4 elements:

-   weightedSIF. Character matrix. Contains a summary of *sifAll*, where
    *Weight* indicates the number of networks that present that
    interaction.

-   nodesAttributes. Character matrix. Contains a summary of
    *attributesAll*, where all columns ending in *Act* indicate the
    number of networks that present that node in each of the states
    (Zero = doesn’t appear in the network, Up = upregulated, Down =
    downregulated, Avg = average regulation). The column *NodeType*
    indicates if the node is a TF’s (T) or a perturbed one (S).

-   sifAll. List of character matrices. All networks that have been
    calculated.

-   attributesAll. List of character matrices.

NOTE: This is valid for current CARNIVAL version (1.0.0)

## CARNIVAL interpretation

The summary files ( *sifAll* and *attributesAll* ) can be directly used
to identify causal interactions between the perturbed nodes and the
selected Transcription Factors. In addition to extracting direct
information from the network, we can run different downstream analysis
based on the necessities of each project.

Here are described some of the downstream approaches that we have used:

-   Pathway enrichment analysis
-   Network comparisons

# Get ready

``` r
library(readr)
library(piano)
library(dplyr)
library(ggplot2)
library(tibble)
library(tidyr)
library(dplyr)
library(scales)
library(plyr)
library(GSEABase)
library(network)
library(reshape2)
library(cowplot)
library(pheatmap)
library(ggraph)
library(tidygraph)
if (!require(RNAscripts)) {
  devtools::install("../scripts/RNAscripts", upgrade = "never", quiet = TRUE)
}
library(RNAscripts)
# set working directory
# setwd("~/Projects/transcriptutorial/scripts")
# source("./support_networks.r")
## We also load the support functions
if (exists("snakemake")) {
  carnival_sample_result <- snakemake@input[["carnival_sample_result"]]
  tt_basepath <- snakemake@params[["tutorial_source_path"]]

  names(carnival_sample_result) <- snakemake@params[["sample_names"]]
  stopifnot(length(carnival_sample_result) > 0)
  workflow_config <- snakemake@config

  sample_table <- readr::read_tsv(workflow_config$samples)
  organism <- workflow_config$organism
} else {
  BASE_DIR <- "/omics/odcf/analysis/OE0228_projects/VascularAging/rna_sequencing/tec_aplnrKO/results/carnival/sample_resolution"
  workflow_config <- yaml::read_yaml("../../../configs/tec_aplnrKO.yaml")
  sample_table <- readr::read_tsv("../../../data/tumor_vs_ec/Aplnr_KO.tsv")
  sample_ids <- sample_table %>% pull(sample)
  carnival_sample_result <- file.path(BASE_DIR, paste0(sample_ids, "_carn_res.RDS.gz"))
  names(carnival_sample_result) <- sample_ids
  base_path <- "/desktop-home/heyer/projects/Vascular_Aging"
  tt_basepath <- file.path(base_path, "RNAseq/rna-seq-star-deseq2/workflow/scripts/transcriptutorial")
  organism <- "Mus musculus"
}
source(file.path(tt_basepath, "support_enrichment.r"))
source(file.path(tt_basepath, "support_networks.r"))

carnival_sample_resolution <- purrr::map(carnival_sample_result, readRDS)

# Check if any are 0

length_index <- purrr::map(carnival_sample_resolution, length)
# Check if any are null
lgl_index <- purrr::map_int(carnival_sample_resolution, length) == 0
carnival_sample_resolution <- carnival_sample_resolution[!lgl_index]
stopifnot(length(carnival_sample_resolution) > 0)
process_carnival <- function(carnival_res) {
  carnival_res$weightedSIF <- carnival_res$weightedSIF %>% dplyr::filter(Weight != 0)
  #  carnival_res$nodesAttributes <- carnival_res$nodesAttributes %>% dplyr::filter(AvgAct != 0)
  carnival_res
}

carnival_sample_resolution <- purrr::map(carnival_sample_resolution, process_carnival)
```

WARNING: The following samples could not create a network.
VascAge_aplnrtec9, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
NA, NA

``` r
# read CARNIVAL results
omnipath_id <- RNAscripts::get_organism_omnipath_id(organism)
omniR <- OmnipathR::import_omnipath_interactions(organism = omnipath_id)

# signed and directed
omnipath_sd <- omniR %>% dplyr::filter(consensus_direction == 1 &
  (consensus_stimulation == 1 |
    consensus_inhibition == 1
  ))

# changing 0/1 criteria in consensus_stimulation/inhibition to -1/1
omnipath_sd$consensus_stimulation[which(omnipath_sd$consensus_stimulation == 0)] <- -1
omnipath_sd$consensus_inhibition[which(omnipath_sd$consensus_inhibition == 1)] <- -1
omnipath_sd$consensus_inhibition[which(omnipath_sd$consensus_inhibition == 0)] <- 1

# check consistency on consensus sign and select only those in a SIF format
sif <- omnipath_sd[, c("source_genesymbol", "consensus_stimulation", "consensus_inhibition", "target_genesymbol")] %>%
  dplyr::filter(consensus_stimulation == consensus_inhibition) %>%
  unique.data.frame()

sif$consensus_stimulation <- NULL
colnames(sif) <- c("source", "interaction", "target")

# remove complexes
sif$source <- gsub(":", "_", sif$source)
sif$target <- gsub(":", "_", sif$target)
```

``` r
# get only summary files from CARNIVAL results
sifts <- lapply(carnival_sample_resolution, function(x) {
  x$weightedSIF
})
nodos <- lapply(carnival_sample_resolution, function(x) {
  x$nodesAttributes
})
write(names(sifts), file = stderr())
# Calculate the number of edges and nodes in the networks and its density
node_edge <- do.call(rbind, lapply(sifts, count_edges_nodes_degree))
```

    ## There are 25 extra interactions with different sign

    ## There are 24 extra interactions with different sign

``` r
# Calculate degree distribution for a sample
yeet <- do.call(rbind, sifts) %>% unique()
count_degree <- yeet %>% degree_count()

# degree distribution
p <- data.frame(table(count_degree$total_count) / nrow(count_degree))
colnames(p) <- c("Var1", "total_degree")
p <- merge.data.frame(p, data.frame(table(count_degree$in_count) / nrow(count_degree)), all = T)
colnames(p) <- c("Var1", "total_degree", "in_degree")
p <- merge.data.frame(p, data.frame(table(count_degree$out_count) / nrow(count_degree)), all = T)
colnames(p) <- c("k", "total_degree", "in_degree", "out_degree")
p <- melt(p, value.name = "p", id.vars = "k")
p$k <- relevel(p$k, "0")

# visualise
ggdat <- as.data.frame(node_edge) %>%
  tibble::rownames_to_column(var = "sample") %>%
  dplyr::mutate(condition = gsub(".Rep[0-9]{1}", "", sample))

# Plotting

# relation between number of edges and nodes
ggplot(ggdat, aes(x = nodes, y = edges, color = as.factor(condition))) +
  geom_point() +
  geom_text(
    label = ggdat$sample,
    check_overlap = TRUE,
    vjust = 0,
    nudge_y = 0.5,
    show.legend = F
  ) +
  theme_bw(base_size = 15) +
  guides(color = guide_legend(title = "Conditions")) +
  ggtitle("Node-edge composition")
```

![](06_analysis_CARNIVAL_results_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
# degree distribution
levels(p$k) <- levels(p$k) %>%
  as.numeric() %>%
  sort()
dd <- ggplot(data = p, aes(x = k, y = p, group = variable, color = variable)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  theme(legend.position = "none") +
  guides(color = guide_legend(title = "degree type")) +
  ggtitle("Degree distribution")

ddp <- ggplot(data = p, aes(x = as.numeric(k), y = p, group = variable, color = variable)) +
  geom_point() +
  geom_line() +
  scale_x_continuous(
    breaks = as.numeric(p$k),
    trans = scales::log_trans()
  ) +
  scale_y_log10(
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x))
  ) +
  annotation_logticks() +
  theme_bw() +
  guides(color = guide_legend(title = "degree type")) +
  ggtitle("Degree distribution (log scale)") +
  xlab("k (ln)") +
  ylab("p (log10)")

plot_grid(dd, ddp, labels = "auto", rel_widths = c(1, 2))
```

    ## Warning: Removed 20 rows containing missing values (`geom_point()`).

    ## Warning: Removed 6 rows containing missing values (`geom_line()`).

    ## Warning: Removed 20 rows containing missing values (`geom_point()`).

    ## Warning: Removed 6 rows containing missing values (`geom_line()`).

![](06_analysis_CARNIVAL_results_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->
\# Enrichment Analysis

We define two different gene sets in order to conduct the enrichment
analysis. The first set contains the nodes that appear in the CARNIVAL
output, and are therefore relevant in the context of our input
transcriptomic data. The second set contains all the genes in our prior
knowledge network which are used as the background.

We also downloaded from **MSigDB <https://www.gsea-msigdb.org/> the
following dataset: c2.cp.v7.1.symbols.gmt**. It contains several
pathways from different resources and the genes that are known to be
involved in those pathways.

Using the *Piano* R package, we run a gene set analysis (GSA) based on a
list of significant genes (CARNIVAL nodes) and a gene set collection
(background). It uses Fisher’s exact test.

# Topological measurments

We can get a sense of the size of the reconstructed networks using
topological parameters, such as number of edges and nodes. We can also
have a look at the network’s density and the *degree distribution*. The
*density* indicates the proportion of interactions that exist in our
network when comparing with all possible interactions that can be
established. The *degree distribution* shows the number of connections
of a node. In a directed network, we can distinguish between incoming
and outgoing connections.

NOTE: Here the density is calculated for a *directed graph*. As CARNIVAL
can report 2 interactions between the same 2 nodes with different sign,
these “doubled” interactions are excluded when calculating the density.

For this section, and the following ones, we use the CARNIVAL networks
obtained when running the *sample_resolution_carnival.R*.

To find out more about this topic:

-   <https://mathinsight.org/degree_distribution>
-   <https://www.networksciencebook.com/chapter/2#degree>

``` r
get_common_interactions <- function(interaction_list, psmpl_per = 95,
                                    s_table = sample_table,
                                    condition_sel) {
  sample_id_list <- s_table %>%
    filter(condition %in% condition_sel) %>%
    pull(sample)
  interaction_list <- interaction_list[, sample_id_list]
  
  if (!is.null(dim(interaction_list))) {
    shared_interactions_WT <- getCoreInteractions(topology = interaction_list, psmpl = psmpl_per)
    
    # Visualise the interactions
    colnames(shared_interactions_WT) <- c("from", "sign", "to")
    labels_edge <- c("-1" = "inhibition", "1" = "activation")
    nodes <- data.frame(union(shared_interactions_WT$from, shared_interactions_WT$to))
    colnames(nodes) <- "nodes"
    nodes$label <- nodes$nodes
    
    gg_graph <- tidygraph::tbl_graph(nodes = nodes, edges = shared_interactions_WT) %>%
      ggraph::ggraph(layout = "auto") +
      geom_node_point(color = "#C0C0C0", size = 8) +
      geom_edge_link(arrow = arrow(), aes(edge_colour = as.factor(sign))) +
      theme_graph() +
      geom_node_text(aes(label = label), vjust = 0.4)
  } else {
    gg_graph <- NULL
  }
  gg_graph
}
```

Generally speaking, biological networks are not dense, so don’t worry if
the density values are low… they usually are!

The degree distribution is an interesting graphic to look at, as we can
get some immediate information:

-   k = 0 for in-degree indicates the proportion of initial nodes, while
    for out-degree indicates the effectors (here the TFs). For
    total-degree, that would indicate the presence of isolated nodes.

-   Most of the nodes do not have many connections (e.g. k = 1 -\> p =
    0.6), but there are a few that are highly connected (e.g. k \> 6).

-   We can find some hubs when k is higher (plot B). The out-degree ends
    at k = 5; this means that the biggest regulatory-hub regulates at
    most 5 other nodes. In a similar way, the in-degree goes up to k =
    9; This means that there are few hubs (k \> 5) that are regulated by
    at most 9 nodes.

# Network comparison

When we have more than one network, we usually would like to know how
(dis)similar these networks are.

We can use the *Jaccard Index* to measure similarities and diversity
between sample sets (<https://en.wikipedia.org/wiki/Jaccard_index>).
This index can be used in different scenarios. Generally, it’s used with
the nodes, but we can also use it with the edges.

When we have groups, in our case WT and FOXA2KO, we can also extract the
interactions that are common.

``` r
# create a matrix of all interactions for all samples
write(dim(sif), file = stderr())
interactions <- getTopology(networks = sifts, scafoldNET = sif)
colnames(interactions) <- names(carnival_sample_resolution)
# FIxes bug in topology function (To lazy to actually fix)
ncol_interact <- ncol(interactions)
# interactions <- interactions[rowSums(interactions) > 0,]
# interactions <- interactions[,-c(1:ncol_interact/2)]
# get the edges per sample
# get the edges per sample
net_int <- apply(interactions, 2, function(x, r) {
  r[which(!is.na(x))]
}, rownames(interactions))

# calculate Jaccard indexes per pair
combined <- expand.grid(1:length(names(sifts)), 1:length(names(sifts)))
jac_index <- matrix(
  data = NA, nrow = length(names(sifts)), ncol = length(names(sifts)),
  dimnames = list(names(sifts), names(sifts))
)

for (i in 1:nrow(combined)) {
  n <- names(sifts)[combined[i, 1]]
  m <- names(sifts)[combined[i, 2]]
  jac_index[n, m] <- length(intersect(net_int[[n]], net_int[[m]])) / length(union(net_int[[n]], net_int[[m]]))
}

# Visualize the indexes in a heatmap

pheatmap::pheatmap(jac_index,
  fontsize = 14,
  fontsize_row = 10, fontsize_col = 10,
  angle_col = 45, treeheight_col = 0
)
```

![](06_analysis_CARNIVAL_results_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
corrplot::corrplot(jac_index, order = "hclust")
```

![](06_analysis_CARNIVAL_results_files/figure-gfm/unnamed-chunk-4-2.png)<!-- -->

## Inferring mixed networks from samples

Here we investiagte which edges are shared between samples to get a
bettter look at the structure of these. networks

We plot Edges which we can find in at least 60% of samples in one
comparison group.

We compare both in one sample group, but we also have a look which edges
are consisten across groups (noted with Combined) aka. combined_young ==
young-ec + young-EC_apeliun

``` r
groups_to_check <- c(
  workflow_config$comp_groups,
  as.list(setNames(
    object = sample_table$condition %>% unique(),
    nm = sample_table$condition %>% unique()
  ))
)


yoted <- purrr::map(groups_to_check, get_common_interactions, interaction_list = interactions, psmpl_per = 60, s_table = sample_table[!lgl_index, ])
```

    ## 256 interactions found in at least 0 samples out of 0

    ## Using `stress` as default layout

    ## 256 interactions found in at least 0 samples out of 0

    ## Using `stress` as default layout

``` r
if(any(purrr::map(yoted, is.null))) {
  null_index <- purrr::map_lgl(yoted, is.null)
  yoted <- yoted[!null_index]
}
```

    ## Warning in any(purrr::map(yoted, is.null)): coercing argument of type 'list' to
    ## logical

``` r
yoted <- purrr::map2(yoted, names(yoted), function(x, y) x + ggtitle(y))

yoted
```

    ## $tumor_control

    ## Warning: Using the `size` aesthetic in this geom was deprecated in ggplot2 3.4.0.
    ## ℹ Please use `linewidth` in the `default_aes` field and elsewhere instead.

![](06_analysis_CARNIVAL_results_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

    ## 
    ## $ec_KO

![](06_analysis_CARNIVAL_results_files/figure-gfm/unnamed-chunk-5-2.png)<!-- -->

``` r
build_matrix_from_nodes <- function(node_list, node_type = "AvgAct") {
  gene_index <- purrr::map(node_list, ~ pull(., Node)) %>%
    unlist() %>%
    unique()
  node_mat <- purrr::map(node_list, ~ dplyr::select(., Node, !!node_type)) %>% purrr::reduce(full_join, by = "Node")
  colnames(node_mat) <- c("Node", names(node_list))
  node_mat
}

avg_mat <- build_matrix_from_nodes(nodos)
rownames(avg_mat) <- avg_mat$Node
```

    ## Warning: Setting row names on a tibble is deprecated.

``` r
avg_mat <- subset(avg_mat, select = -c(Node)) %>% as.matrix()
non_zero_index <- apply(avg_mat, 1, function(x) length(which(x != 0)) >= 2)



ComplexHeatmap::Heatmap(avg_mat[non_zero_index, ], column_names_rot = 45, row_names_gp = grid::gpar(fontsize = 6), cluster_columns = FALSE)
```

![](06_analysis_CARNIVAL_results_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->
\## Session Info Details

    ## R version 4.2.0 (2022-04-22)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: CentOS Linux 7 (Core)
    ## 
    ## Matrix products: default
    ## BLAS/LAPACK: /home/heyer/mambaforge/lib/libopenblasp-r0.3.21.so
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## attached base packages:
    ## [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
    ## [8] methods   base     
    ## 
    ## other attached packages:
    ##  [1] snowfall_1.84-6.2    snow_0.4-4           RNAscripts_0.1.0    
    ##  [4] tidygraph_1.2.1      ggraph_2.0.5         pheatmap_1.0.12     
    ##  [7] cowplot_1.1.1        reshape2_1.4.4       network_1.18.0      
    ## [10] GSEABase_1.58.0      graph_1.74.0         annotate_1.74.0     
    ## [13] XML_3.99-0.12        AnnotationDbi_1.60.0 IRanges_2.30.0      
    ## [16] S4Vectors_0.34.0     Biobase_2.56.0       BiocGenerics_0.42.0 
    ## [19] plyr_1.8.8           scales_1.2.1         tidyr_1.2.1         
    ## [22] tibble_3.1.8         ggplot2_3.4.0        dplyr_1.0.10        
    ## [25] piano_2.12.1         readr_2.1.3         
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] circlize_0.4.15        readxl_1.4.1           backports_1.4.1       
    ##   [4] fastmatch_1.1-3        corrplot_0.92          systemfonts_1.0.4     
    ##   [7] igraph_1.3.5           shinydashboard_0.7.2   BiocParallel_1.30.3   
    ##  [10] GenomeInfoDb_1.32.2    digest_0.6.30          foreach_1.5.2         
    ##  [13] htmltools_0.5.3        viridis_0.6.2          OmnipathR_3.5.25      
    ##  [16] fansi_1.0.3            magrittr_2.0.3         checkmate_2.1.0       
    ##  [19] memoise_2.0.1          doParallel_1.0.17      cluster_2.1.3         
    ##  [22] tzdb_0.3.0             limma_3.52.2           ComplexHeatmap_2.12.1 
    ##  [25] Biostrings_2.64.0      graphlayouts_0.8.0     matrixStats_0.63.0    
    ##  [28] vroom_1.6.0            svglite_2.1.0          prettyunits_1.1.1     
    ##  [31] colorspace_2.0-3       blob_1.2.3             rvest_1.0.3           
    ##  [34] rappdirs_0.3.3         ggrepel_0.9.2          xfun_0.35             
    ##  [37] crayon_1.5.2           RCurl_1.98-1.7         jsonlite_1.8.3        
    ##  [40] iterators_1.0.14       glue_1.6.2             polyclip_1.10-4       
    ##  [43] gtable_0.3.1           zlibbioc_1.42.0        XVector_0.36.0        
    ##  [46] GetoptLong_1.0.5       shape_1.4.6            DBI_1.1.3             
    ##  [49] relations_0.6-12       Rcpp_1.0.9             viridisLite_0.4.1     
    ##  [52] xtable_1.8-4           progress_1.2.2         clue_0.3-62           
    ##  [55] bit_4.0.5              DT_0.26                htmlwidgets_1.5.4     
    ##  [58] httr_1.4.4             fgsea_1.22.0           gplots_3.1.3          
    ##  [61] RColorBrewer_1.1-3     ellipsis_0.3.2         pkgconfig_2.0.3       
    ##  [64] farver_2.1.1           utf8_1.2.2             labeling_0.4.2        
    ##  [67] tidyselect_1.2.0       rlang_1.0.6            later_1.3.0           
    ##  [70] munsell_0.5.0          cellranger_1.1.0       tools_4.2.0           
    ##  [73] visNetwork_2.1.2       cachem_1.0.6           cli_3.4.1             
    ##  [76] generics_0.1.3         RSQLite_2.2.14         statnet.common_4.7.0  
    ##  [79] evaluate_0.18          stringr_1.4.1          fastmap_1.1.0         
    ##  [82] yaml_2.3.6             knitr_1.41             bit64_4.0.5           
    ##  [85] caTools_1.18.2         purrr_0.3.5            KEGGREST_1.36.2       
    ##  [88] mime_0.12              slam_0.1-50            xml2_1.3.3            
    ##  [91] compiler_4.2.0         rstudioapi_0.14        curl_4.3.3            
    ##  [94] png_0.1-7              marray_1.74.0          tweenr_1.0.2          
    ##  [97] stringi_1.7.8          highr_0.9              logger_0.2.2          
    ## [100] lattice_0.20-45        Matrix_1.5-3           shinyjs_2.1.0         
    ## [103] vctrs_0.5.1            pillar_1.8.1           lifecycle_1.0.3       
    ## [106] GlobalOptions_0.1.2    data.table_1.14.6      bitops_1.0-7          
    ## [109] httpuv_1.6.6           R6_2.5.1               promises_1.2.0.1      
    ## [112] KernSmooth_2.23-20     gridExtra_2.3          codetools_0.2-18      
    ## [115] MASS_7.3-57            gtools_3.9.4           assertthat_0.2.1      
    ## [118] rjson_0.2.21           withr_2.5.0            GenomeInfoDbData_1.2.8
    ## [121] hms_1.1.2              grid_4.2.0             coda_0.19-4           
    ## [124] rmarkdown_2.18         Cairo_1.6-0            sets_1.0-21           
    ## [127] ggforce_0.3.3          shiny_1.7.3
