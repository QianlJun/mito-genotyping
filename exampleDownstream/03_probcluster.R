library(ClusterR)
library(mgatk)
library(SummarizedExperiment)
library(ggplot2)
library(BuenColors)
library(MultiAssayExperiment)
library(Matrix)
library(dplyr)

set.seed(1234)

"%ni%" <- Negate("%in%")

raw <- readRDS("../data/mgatk_scRNAseq/cc100_mgatk.rds")

covSE <- raw[["coverage"]]
allSE <- raw[["alleles"]]
rm(raw)
covCell <-  Matrix::colMeans(assays(covSE)[["coverage"]])
scrna <- readRDS("../output/scRNAseq_transcripts.rds")

# Allele Frequency
hits <- assays(allSE)[["counts"]][,covCell > 100 & colnames(assays(allSE)[["counts"]]) %in% colnames(scrna)]
total <- assays(covSE)[["coverage"]][start(rowRanges(allSE)),] [,covCell>100]
rownames(hits) <- paste0(data.frame(rowRanges(allSE))[,c(2)], "_", data.frame(rowRanges(allSE))[,c(7)])
rownames(total) <- paste0(data.frame(rowRanges(allSE))[,c(2)], "_", data.frame(rowRanges(allSE))[,c(7)])

varDF <- readRDS("../output/usefulVariants.rds")
vars <- rownames(varDF[varDF$AF_ATAC > 0.0001,])

# Filter Variants
hits <- hits[vars,] 
total <- total[vars,] 

af <- hits/(total + 0.0001)
dim(af)
dim(af)
paf <- data.matrix(af[rowSums(af > 0.2) > 2,])
paf[paf > 0.2] <- 0.2
paf[paf < 0.01] <- 0.00

#fz <- Cluster_Medoids(t(paf), clusters = 11, fuzzy = TRUE)
#saveRDS(as.character(fz$clusters), "../output/cc100_scRNAseq_clusters.rds")

clusters <- readRDS("../output/cc100_scRNAseq_clusters.rds")
annodf <- data.frame(clusters = factor(clusters, levels = as.character(c(1,9, 10, 3, 11, 4:8, 2, 12))),
                     ID = colnames(paf)
)

annodf2 <- annodf %>% arrange(clusters)
cols <- c(jdb_palettes[["brewer_spectra"]], "black", "grey")
names(cols) <- unique(annodf2$clusters)

ha1 <- HeatmapAnnotation(df = annodf2[,c("clusters"), drop = FALSE],
                         col = list(clusters = cols)
)

orderVars <- c("2044_C", "5317_A", "665_G", "2253_G", "13360_A", "14471_C", "9554_A",  "10103_C", "5948_G",
                "5056_C",  "15173_A", "7318_G", "204_C", "14069_C")

orderVars[orderVars %ni% rownames(paf)]
rownames(paf)[rownames(paf) %ni% orderVars]

pdf(file="../output/plots/fuzzy_wuzzy_RNA.pdf", height = 2, width = 5)  
par(cex.main=0.8,mar=c(1,1,1,1))
Heatmap(data.matrix(paf[orderVars,annodf2$ID]), col=as.character(jdb_palette("brewer_red",type="continuous")),
        cluster_rows = FALSE, cluster_columns = FALSE, show_column_names = FALSE,
        row_names_gp = gpar(fontsize =4),
        top_annotation = ha1,
        name = "Allele\nFrequency")
dev.off()
