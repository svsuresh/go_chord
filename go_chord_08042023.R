## Load libraries
suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library(GOplot))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(janitor))
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(org.Mm.eg.db))
#######################
## Updated code#########
## ####################

df = read.csv("example_input.tsv",
    sep = "\t",
    header = T
)
names(df)=c("symbol","fc","p_val")

df <- df[abs(df$fc) > 0.6 & df$p_val < 0.05, ]

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

enrich_df <- function(x){
    enrichGO(
        x$symbol,
        "org.Mm.eg.db",
        keyType = "SYMBOL",
        ont = "ALL",
        pvalueCutoff = 0.05,
        pAdjustMethod = "BH"
        )@result %>%
    filter(ONTOLOGY == "BP") %>%
    arrange(p.adjust) %>%
    slice_head(n = 10) %>%
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

merge_df <- function(x,y){
    x %>%
        separate_longer_delim(Genes, ",") %>%
        inner_join(y, by = c("Genes" = "symbol")) %>%
        group_by(ID, Category) %>%
        slice_max(order_by = abs(fc), n = 10) %>%
        ungroup() %>%
        as.data.frame()
}

df2 <- clean_df(df)
enrich_df2 <- enrich_df(df2)
merged_df <- merge_df(enrich_df2,df2)

## Merge two data frames two to extract final expression data frame. This data frame contains 10 statistically significant GO BP and 10 most affected (by logFC) genes. Extract genes, logFC and p-value. Change first two columns and remove duplicate rows. Convert gene symbols to upper case. This is very important. Upper case is hard coded in the code
final_df <- function(x){
    x %>%
        dplyr::select(Genes,fc,p_val) %>%
        dplyr::rename("ID" = 1, "logFC" = 2) %>%
        distinct() %>%
        mutate(ID = toupper(ID)) %>%
        as.data.frame()
}

gene_df <- final_df(merged_df)

## draw diagram.
circ <- circle_dat(enrich_df2, gene_df)
chord <- chord_dat(circ, gene_df, enrich_df2$ID)
GOChord(chord, space=0.02)

# save.image("go_plot.Rdata")


