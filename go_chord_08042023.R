## Load libraries
suppressPackageStartupMessages(library(GOplot))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(janitor))
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(here))
#######################
## Updated code#########
## ####################
options <- commandArgs(trailingOnly = TRUE)

df = read.csv(options[1],
    sep = "\t",
    header = T
)

###All functions########
## Clean the data frame for empty gene symbols, identical gene symbols with different p/fc values
clean_df <- function(x) {
    x %>%
        filter(symbol != "") %>%
        separate_longer_delim(symbol, "///") %>%
        group_by(symbol) %>%
        slice_max(order_by = fc, n = 1) %>%
        slice_min(order_by = p_val, n = 1) %>%
        ungroup() %>%
        as.data.frame()
}

## Function of enriching GO terms for all ontologies and format the ouput from enrichment to a data frame
enrich_GO <- function(x){
    enrichGO(
        x$symbol,
        Org,
        keyType = "SYMBOL",
        ont = "ALL",
        pvalueCutoff = 0.05,
        pAdjustMethod = "BH"
        )@result %>%
    filter(ONTOLOGY == "BP") %>%
    arrange(p.adjust) %>%
        # slice_head(n = 10) %>%
    slice_head(n = as.integer(options[3])) %>%
    dplyr::select(1:3,7,9) %>%
    ungroup() %>%
    mutate(across(geneID, ~str_replace_all(.x, "/", ','))) %>%
    dplyr::select(1:3, 5, 4) %>%
    dplyr::rename(
        "Category" = 1,
        "ID" = 2,
        "Term" = 3,
        "Genes" = 4,
        "adj_pval" = 5
    ) %>%
    as.data.frame()
}

## Merge dataframes: GO enrichment data and cleaned data
merge_df <- function(x,y){
    x %>%
        separate_longer_delim(Genes, ",") %>%
        inner_join(y, by = c("Genes" = "symbol")) %>%
        group_by(ID, Category) %>%
        # slice_max(order_by = abs(fc), n = 10) %>%
       slice_max(order_by = abs(fc), n = as.integer(options[4])) %>%
        ungroup() %>%
        as.data.frame()
}

## Final data format for drawing the chord
final_df <- function(x){
    x %>%
        dplyr::select(Genes,fc,p_val) %>%
        dplyr::rename("ID" = 1, "logFC" = 2) %>%
        distinct() %>%
        mutate(ID = toupper(ID)) %>%
        as.data.frame()
}

#### Build the organism of interest#####
# Select the organism to load db
Org_db <- switch(options[2], "mouse" = "Mm", "rat" = "Rn", "human" = "Hs")

# Build the package name
Org <- paste0("org.",Org_db,".eg.db")
# Org <- "org.Mm.eg.db"
#Load the package name
suppressPackageStartupMessages(require(Org, character.only = T))

### all the work ####
## Final data frame contains 10 statistically significant GO BP and 10 most affected (by logFC) genes. Extract genes, logFC and p-value. Change first two columns and remove duplicate rows. Convert gene symbols to upper case. This is very important. Upper case is hard coded in the code

names(df)=c("symbol","fc","p_val")
df <- df[abs(df$fc) > 0.6 & df$p_val < 0.05, ]
df2 <- clean_df(df)
enrich_df2 <- enrich_GO(df2)
merged_df <- merge_df(enrich_df2,df2)
gene_df <- final_df(merged_df)

## draw diagram.
circ <- circle_dat(enrich_df2, gene_df)
chord <- chord_dat(circ, gene_df, enrich_df2$ID)
GOChord(chord)
# ggsave("test.pdf", dpi = 300, width = 1800, height = 1800, units = "px")
# save.image("go_plot.Rdata")
GOChord

