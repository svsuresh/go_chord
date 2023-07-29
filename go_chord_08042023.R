## Load libraries
library(GOplot)
library(dplyr)
library(tidyr)
library(stringr)
library(janitor)
library(clusterProfiler)
library(org.Mm.eg.db)
#######################
## Updated code#########
## ####################
df = read.csv("example_input.tsv",
    sep = "\t",
    header = T
)
df <- df[abs(df$logFC) > 0.6 & df$adj.P.Val < 0.05, ]
names(df)

clean_df <- function(x) {
    x %>%
        dplyr::select(Gene.symbol, logFC, adj.P.Val) %>%
        janitor::clean_names() %>%
        filter(gene_symbol != "") %>%
        separate_longer_delim(gene_symbol, "///") %>%
        group_by(gene_symbol) %>%
        slice_max(order_by = log_fc, n = 1) %>%
        slice_min(order_by = adj_p_val, n = 1) %>%
        ungroup() %>%
        as.data.frame()
}
df2 <- clean_df(df)

ego_df2 <-
    enrichGO(
        df2$gene_symbol,
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
    mutate(across(geneID, str_replace_all, "/", ',')) %>%
    dplyr::select(1:3, 5, 4) %>%
    dplyr::rename(
        "Category" = 1,
        "ID" = 2,
        "Term" = 3,
        "Genes" = 4,
        "adj_pval" = 5
    ) %>%
    as.data.frame()

## Merge two data frames two to extract final expression data frame. This data frame contains 10 statistically significant GO BP and 10 most affected (by logFC) genes. Extract genes, logFC and p-value. Change first two columns and remove duplicate rows. Convert gene symbols to upper case. This is very important. Upper case is hard coded in the code

merged_df <- ego_df2 %>%
    separate_longer_delim(Genes, ",") %>%
    inner_join(df2, by = c("Genes" = "gene_symbol")) %>%
    group_by(ID, Category) %>%
    slice_max(order_by = abs(log_fc), n = 10) %>%
    ungroup() %>%
    as.data.frame()

gene_df <- merged_df %>%
    dplyr::select(Genes, log_fc, adj_p_val) %>%
    dplyr::rename("ID" = 1, "logFC" = 2) %>%
    distinct() %>%
    mutate(ID = toupper(ID)) %>%
    as.data.frame()

go_bp <- ego_df2$ID

## draw diagram.
circ <- circle_dat(ego_df2, gene_df)
chord <- chord_dat(circ, gene_df, go_bp)
GOChord(chord, space=0.02)


# save.image("go_plot.Rdata")
