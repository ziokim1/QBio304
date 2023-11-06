# TxImport ----
library(tidyverse) # provides access to Hadley Wickham's collection of R packages for data science
library(tximport) # package for getting Kallisto results into R
library(biomaRt) # retrieving annotations for tomato

setwd("~/assignment/read_map") # set the directory to where the mapped reads are
targets <- read_tsv("studydesign_assignment.txt") # read in your study design
targets$group2 <- paste(targets$group, targets$stage)

path <- file.path(targets$sample, "abundance.tsv") # set file paths to your mapped data

lyc.anno <- useMart(biomart="plants_mart", dataset = "slycopersicum_eg_gene", host="https://plants.ensembl.org") # annotations for tomato
Tx.lyc <- getBM(attributes=c('ensembl_transcript_id',
                             'ensembl_gene_id', 'description'),
                mart = lyc.anno) %>%
          as_tibble() %>%
          dplyr::rename(target_id = ensembl_transcript_id, gene_name = ensembl_gene_id)  %>%
          dplyr::select("target_id", "gene_name")
Txi_gene <- tximport(path, # annotates target ID to gene ID
                     type = "kallisto", 
                     tx2gene = Tx.lyc, 
                     txOut = T, # data is represented as transcripts
                     countsFromAbundance = "lengthScaledTPM",
                     ignoreTxVersion = T)

# data wrangling ----

library(edgeR) # normalize read counts
library(matrixStats) # filters raw data
library(cowplot) # plot multigraphs

sampleLabels <- targets$sample
myDGEList <- DGEList(Txi_gene$counts)
log2.cpm <- cpm(myDGEList, log=TRUE)

log2.cpm.df <- as_tibble(log2.cpm, rownames = "geneID")
colnames(log2.cpm.df) <- c("geneID", sampleLabels)
log2.cpm.df.pivot <- pivot_longer(log2.cpm.df, # dataframe to be pivoted
                                  cols = CS01:HT12, # column names to be stored as a SINGLE variable
                                  names_to = "samples", # name of that new variable (column)
                                  values_to = "expression") # name of new variable (column) storing all the values (data)

p1 <- ggplot(log2.cpm.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 12, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="unfiltered, non-normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

cpm <- cpm(myDGEList)
keepers <- rowSums(cpm>1)>=3 # the experiment has three replicates for each control/heat and Breaker/Breaker+3 treatments
myDGEList.filtered <- myDGEList[keepers,]

log2.cpm.filtered <- cpm(myDGEList.filtered, log=TRUE)
log2.cpm.filtered.df <- as_tibble(log2.cpm.filtered, rownames = "geneID")
colnames(log2.cpm.filtered.df) <- c("geneID", sampleLabels)
log2.cpm.filtered.df.pivot <- pivot_longer(log2.cpm.filtered.df, # dataframe to be pivoted
                                           cols = CS01:HT12, # column names to be stored as a SINGLE variable
                                           names_to = "samples", # name of that new variable (column)
                                           values_to = "expression") # name of new variable (column) storing all the values (data)

p2 <- ggplot(log2.cpm.filtered.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 12, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="filtered, non-normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

myDGEList.filtered.norm <- calcNormFactors(myDGEList.filtered, method = "TMM")
log2.cpm.filtered.norm <- cpm(myDGEList.filtered.norm, log=TRUE)
log2.cpm.filtered.norm.df <- as_tibble(log2.cpm.filtered.norm, rownames = "geneID")
colnames(log2.cpm.filtered.norm.df) <- c("geneID", sampleLabels)
log2.cpm.filtered.norm.df.pivot <- pivot_longer(log2.cpm.filtered.norm.df, # dataframe to be pivoted
                                                cols = CS01:HT12, # column names to be stored as a SINGLE variable
                                                names_to = "samples", # name of that new variable (column)
                                                values_to = "expression") # name of new variable (column) storing all the values (data)

p3 <- ggplot(log2.cpm.filtered.norm.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 12, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="filtered, TMM normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

plot_grid(p1, p2, p3, labels = c('A', 'B', 'C'), label_size = 18, nrow = 3)

# multivariate ----
library(plotly) # visualize pca plot

pca.res <- prcomp(t(log2.cpm.filtered.norm), scale.=F, retx=T)
screeplot(pca.res) # A screeplot is a standard way to view eigenvalues for each PCA

pc.var <- pca.res$sdev^2 # sdev^2 captures these eigenvalues from the PCA result
pc.per <- round(pc.var/sum(pc.var)*100, 1) 
pca.res.df <- as_tibble(pca.res$x)

group2 <- factor(targets$group2)

pca.plot <- ggplot(pca.res.df) +
  aes(x=PC1, y=PC2, label=sampleLabels, color = group2) +
  geom_point(size=4) +
  xlab(paste0("PC1 (",pc.per[1],"%",")")) + 
  ylab(paste0("PC2 (",pc.per[2],"%",")")) +
  labs(title="PCA plot",
       color = "groups",
       caption=paste0("produced on ", Sys.time())) +
  coord_fixed() +
  theme_bw()

ggplotly(pca.plot)

# differential genes ----
library(limma) # create linear models and decide tests
library(DT) # format tables

# volcano plot for BS+B3
group <- factor(targets$group) # creating a design matrix with independent variables selected for the experiment.
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

v.DEGList.filtered.norm <- voom(myDGEList.filtered.norm, design, plot = FALSE)

fit <- lmFit(v.DEGList.filtered.norm, design)
contrast.matrix <- makeContrasts(difference = control - heat,
                                    levels=design)

fits <- contrasts.fit(fit, contrast.matrix)
ebFit <- eBayes(fits)

myTopHits <- topTable(ebFit, adjust ="BH", coef=1, number=17150, sort.by="logFC")
myTopHits.df <- myTopHits %>%
  as_tibble(rownames = "geneID")

vplot <- ggplot(myTopHits.df) + # plotting the volcano plot
  aes(y=-log10(adj.P.Val), x=logFC, text = paste("Symbol:", geneID)) +
  geom_point(size=2) +
  geom_hline(yintercept = -log10(0.05), linetype="longdash", colour="grey", linewidth=1) +
  geom_vline(xintercept = 3, linetype="longdash", colour="#BE684D", linewidth=1) +
  geom_vline(xintercept = -3, linetype="longdash", colour="#2C467A", linewidth=1) +
  annotate("rect", xmin = 3, xmax = 12, ymin = -log10(0.05), ymax = 7, alpha=.2, fill="#BE684D") +
  annotate("rect", xmin = -3, xmax = -12, ymin = -log10(0.05), ymax = 7, alpha=.2, fill="#2C467A") +
  labs(title="Volcano plot (BS+B3)",
       subtitle = "Solanum lycopersicum",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

# volcano plot for BS
group.BS <- factor(targets$group[targets$stage == "BS"]) # creating a design matrix with independent variables selected for the experiment.
design.BS <- model.matrix(~0 + group.BS)
colnames(design.BS) <- levels(group.BS)

myDGEList.filtered.norm.BS <- myDGEList.filtered.norm
myDGEList.filtered.norm.BS[[1]] <- myDGEList.filtered.norm[[1]][,targets$stage=="BS"]
myDGEList.filtered.norm.BS[[2]] <- myDGEList.filtered.norm[[2]][targets$stage=="BS",]

v.DEGList.filtered.norm.BS <- voom(myDGEList.filtered.norm.BS, design.BS, plot = FALSE)

fit.BS <- lmFit(v.DEGList.filtered.norm.BS, design.BS)
contrast.matrix.BS <- makeContrasts(difference = control - heat,
                                    levels=design.BS)

fits.BS <- contrasts.fit(fit.BS, contrast.matrix.BS)
ebFit.BS <- eBayes(fits.BS)

myTopHits.BS <- topTable(ebFit.BS, adjust ="BH", coef=1, number=17150, sort.by="logFC")
myTopHits.BS.df <- myTopHits.BS %>%
  as_tibble(rownames = "geneID")

vplot.BS <- ggplot(myTopHits.BS.df) + # plotting the volcano plot
  aes(y=-log10(adj.P.Val), x=logFC, text = paste("Symbol:", geneID)) +
  geom_point(size=2) +
  geom_hline(yintercept = -log10(0.05), linetype="longdash", colour="grey", linewidth=1) +
  geom_vline(xintercept = 3, linetype="longdash", colour="#BE684D", linewidth=1) +
  geom_vline(xintercept = -3, linetype="longdash", colour="#2C467A", linewidth=1) +
  annotate("rect", xmin = 3, xmax = 12, ymin = -log10(0.05), ymax = 5, alpha=.2, fill="#BE684D") +
  annotate("rect", xmin = -3, xmax = -12, ymin = -log10(0.05), ymax = 5, alpha=.2, fill="#2C467A") +
  labs(title="Volcano plot (BS)",
       subtitle = "Solanum lycopersicum",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

# volcano plot for B3
group.B3 <- factor(targets$group[targets$stage == "B3"]) # creating a design matrix with independent variables selected for the experiment.
design.B3 <- model.matrix(~0 + group.B3)
colnames(design.B3) <- levels(group.B3)

myDGEList.filtered.norm.B3 <- myDGEList.filtered.norm
myDGEList.filtered.norm.B3[[1]] <- myDGEList.filtered.norm[[1]][,targets$stage=="B3"]
myDGEList.filtered.norm.B3[[2]] <- myDGEList.filtered.norm[[2]][targets$stage=="B3",]

v.DEGList.filtered.norm.B3 <- voom(myDGEList.filtered.norm.B3, design.B3, plot = FALSE)

fit.B3 <- lmFit(v.DEGList.filtered.norm.B3, design.B3)
contrast.matrix.B3 <- makeContrasts(difference = control - heat,
                                    levels=design.B3)

fits.B3 <- contrasts.fit(fit.B3, contrast.matrix.B3)
ebFit.B3 <- eBayes(fits.B3)

myTopHits.B3 <- topTable(ebFit.B3, adjust ="BH", coef=1, number=17150, sort.by="logFC")
myTopHits.B3.df <- myTopHits.B3 %>%
  as_tibble(rownames = "geneID")

vplot.B3 <- ggplot(myTopHits.B3.df) + # plotting the volcano plot
  aes(y=-log10(adj.P.Val), x=logFC, text = paste("Symbol:", geneID)) +
  geom_point(size=2) +
  geom_hline(yintercept = -log10(0.05), linetype="longdash", colour="grey", linewidth=1) +
  geom_vline(xintercept = 3, linetype="longdash", colour="#BE684D", linewidth=1) +
  geom_vline(xintercept = -3, linetype="longdash", colour="#2C467A", linewidth=1) +
  annotate("rect", xmin = 3, xmax = 12, ymin = -log10(0.05), ymax = 5, alpha=.2, fill="#BE684D") +
  annotate("rect", xmin = -3, xmax = -12, ymin = -log10(0.05), ymax = 5, alpha=.2, fill="#2C467A") +
  labs(title="Volcano plot (B3)",
       subtitle = "Solanum lycopersicum",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

plot_grid(vplot, vplot.BS, vplot.B3, label_size = 18, nrow = 1)

# decideTests to pull out the DEGs and make Venn Diagram

results <- decideTests(ebFit, method="global", adjust.method="BH", p.value=0.05, lfc=3)
results.BS <- decideTests(ebFit.BS, method="global", adjust.method="BH", p.value=0.05, lfc=3)
results.B3 <- decideTests(ebFit.B3, method="global", adjust.method="BH", p.value=0.05, lfc=3)

summary(results) # 235 down-regulated and 291 up-regulated genes for BS+B3
summary(results.BS) # 629 down-regulated and 551 up-regulated genes for BS
summary(results.B3) # 180 down-regulated and 254 up-regulated genes for B3

par(mfrow=c(1,3), mai=c(10,10,2,10), cra = 20,
    vennDiagram(results, names="below threshold", include="both", cex = 2), title("BS+B3"), # total of 526 genes are regulated
    vennDiagram(results.BS, names="below threshold", include="both", cex = 2), title("BS"), # total of 1180 genes are regulated
    vennDiagram(results.B3, names="below threshold", include="both", cex = 2), title("B3")) # total of 434 genes are regulated


# create table for candidate genes that are differentially expressed

results <- decideTests(ebFit, method="global", adjust.method="BH", p.value=0.05, lfc=3)
colnames(v.DEGList.filtered.norm$E) <- sampleLabels
diffGenes <- v.DEGList.filtered.norm$E[results[,1] !=0,]
diffGenes.df <- as_tibble(diffGenes, rownames = "geneID")
datatable(diffGenes.df,
          extensions = c('KeyTable', "FixedHeader"),
          caption = 'DEGs in Solanum lycopersicum',
          options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"))) %>%
  formatRound(columns=c(2:13), digits=2)

write_tsv(diffGenes.df,"DiffGenes.txt")

# functional enrichment ----
library(gprofiler2) #tools for accessing the GO enrichment results using g:Profiler web resources

myTopHits <- topTable(ebFit, adjust ="BH", coef=1, number=526, p.value=0.05, lfc=3, sort.by="logFC")
myTopHits.BS <- topTable(ebFit.BS, adjust ="BH", coef=1, number=1180, p.value=0.05, lfc=3, sort.by="logFC")
myTopHits.B3 <- topTable(ebFit.B3, adjust ="BH", coef=1, number=434, p.value=0.05, lfc=3, sort.by="logFC")

# use the 'gost' function from the gprofiler2 package to run GO enrichment analysis
gost.res <- gost(rownames(myTopHits), organism = "slycopersicum", correction_method = "fdr")
gost.res.BS <- gost(rownames(myTopHits.BS), organism = "slycopersicum", correction_method = "fdr")
gost.res.B3 <- gost(rownames(myTopHits.B3), organism = "slycopersicum", correction_method = "fdr")

# produce an interactive manhattan plot of enriched GO terms
mygostplot <- gostplot(gost.res, interactive = F, capped = F) #set interactive=FALSE to get plot for publications
mygostplot.BS <- gostplot(gost.res.BS, interactive = F, capped = F) #set interactive=FALSE to get plot for publications
mygostplot.B3 <- gostplot(gost.res.B3, interactive = F, capped = F)

plot_grid(mygostplot, mygostplot.BS, mygostplot.B3, labels = c('Manhattan plot (BS+B3)', 'Manhattan plot (BS)', 'Manhattan plot (B3)'), label_size = 20, ncol = 1)

# produce a publication quality static manhattan plot with specific GO terms highlighted
# rerun the above gostplot function with 'interactive=F' and save to an object 'mygostplot'
g1 <- publish_gostplot(
  mygostplot.BS, #your static gostplot from above
  highlight_terms = gost.res.BS[[1]][["term_id"]][order(gost.res.BS[[1]][["p_value"]])][1:10], # 10 ontologies with the lowest p-value
  filename = NULL,
  width = NA,
  height = NA)

g2 <- publish_gostplot(
  mygostplot.BS, #your static gostplot from above
  highlight_terms = gost.res.BS[[1]][["term_id"]][order(gost.res.BS[[1]][["term_size"]], decreasing = T)][1:10], # 10 ontologies with the highest term size
  filename = NULL,
  width = NA,
  height = NA)

plot_grid(g1, g2, labels = c('sort by p value (BS)', 'sort by term size (BS)'), label_size = 18, ncol = 1)

# writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
