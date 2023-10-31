# Create a heatmap of differentially expressed genes
library(heatmaply)

results.BS <- decideTests(ebFit.BS, method="global", adjust.method="BH", p.value=0.05, lfc=1)
colnames(v.DEGList.filtered.norm.BS$E) <- sampleLabels.BS
diffGenes.BS <- v.DEGList.filtered.norm.BS$E[results.BS[,1] !=0,]
diffGenes.BS.df <- as_tibble(diffGenes.BS, rownames = "geneID")

hm.BS <- heatmaply(diffGenes.BS.df[1:50,2:7], 
                   dendrogram = "column",
                   xlab = "Samples", ylab = "DEGs", 
                   main = "DEGs in Solanum lycopersicum",
                   scale = "none", #column / if the differential gene is only important, not the extent
                   margins = c(60,100,40,20),
                   grid_color = "white",
                   grid_width = 0.0000001,
                   titleX = T,
                   titleY = T,
                   hide_colorbar = FALSE,
                   branches_lwd = 0.1,
                   label_names = c("Gene", "Sample:", "Value"),
                   fontsize_row = 5, fontsize_col = 5,
                   labCol = colnames(diffGenes.BS.df)[2:7],
                   labRow = diffGenes.BS.df$geneID[1:50], # scaled
                   heatmap_layers = theme(axis.line=element_blank())
)

results.B3 <- decideTests(ebFit.B3, method="global", adjust.method="BH", p.value=0.05, lfc=5)
colnames(v.DEGList.filtered.norm.B3$E) <- sampleLabels.B3
diffGenes.B3 <- v.DEGList.filtered.norm.B3$E[results.B3[,1] !=0,]
diffGenes.B3.df <- as_tibble(diffGenes.B3, rownames = "geneID")

hm.B3 <- heatmaply(diffGenes.B3.df[1:50,2:7], 
                   dendrogram = "column",
                   xlab = "Samples", ylab = "DEGs", 
                   main = "DEGs in Solanum lycopersicum",
                   scale = "none", #column / if the differential gene is only important, not the extent
                   margins = c(60,100,40,20),
                   grid_color = "white",
                   grid_width = 0.0000001,
                   titleX = T,
                   titleY = T,
                   hide_colorbar = FALSE,
                   branches_lwd = 0.1,
                   label_names = c("Gene", "Sample:", "Value"),
                   fontsize_row = 5, fontsize_col = 5,
                   labCol = colnames(diffGenes.B3.df)[2:7],
                   labRow = diffGenes.B3.df$geneID[1:50], # scaled
                   heatmap_layers = theme(axis.line=element_blank())
)

hm.BS
hm.B3