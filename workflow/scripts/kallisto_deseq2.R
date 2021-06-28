#Carregar as bibliotecas instaladas. Deve ser executado no inicio de toda a análise
library("IHW")
library("DESeq2")
library("ggplot2")
library("RColorBrewer")
library("pheatmap")
library("tximport")
#library("biomart")

print(snakemake@output[[1]])

#Carrega os nomes dos arquivos de expressão de cada biblioteca
tx2gene <- read.csv(snakemake@input[["tx2gene"]], header = FALSE, sep="\t")
head(tx2gene)
print("\n\n")
head(snakemake@params[["samples_file"]])
print("\n\n")
samples_table <- read.csv(snakemake@params[["samples_file"]], header = TRUE, sep=",")
samples <- samples_table$sample
files <- file.path("results/kallisto",samples, "abundance.tsv")
names(files) <- samples
files
print("\n\n")
#Importa os resultados de expressão de cada biblioteca
txi.kallisto.tsv <- tximport(files, type = "kallisto",  tx2gene = tx2gene, dropInfReps = TRUE)

colnames(txi.kallisto.tsv$counts) <- samples
head(txi.kallisto.tsv$counts)

colDatak <-data.frame(row.names=colnames(txi.kallisto.tsv$counts), type=as.factor(samples_table$condition))
colDatak

ddsk <- DESeqDataSetFromTximport(txi = txi.kallisto.tsv, colData = colDatak, design = ~ type)

#Cria objeto do DESeq2
ddsk <- DESeq(ddsk)
print("Criar gráfico do Heatmap")
#Cria gráfico do Heatmap
vsdk <- rlog(ddsk, blind=FALSE)
sampleDistsk <- dist(t(assay(vsdk)))
sampleDistMatrixk <- as.matrix(sampleDistsk)
rownames(sampleDistMatrixk) <- vsdk$type
colnames(sampleDistMatrixk) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
print("Heatmap")
pheatmap(sampleDistMatrixk,
         clustering_distance_rows=sampleDistsk,
         clustering_distance_cols=sampleDistsk,
         col=colors)

#Cria PCA
pcaDatak <- plotPCA(vsdk, intgroup=c("type"), returnData=TRUE)
percentVark <- round(100 * attr(pcaDatak, "percentVar"))
ggplot(pcaDatak, aes(PC1, PC2, color=type, shape=type)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVark[1],"% variance")) +
  ylab(paste0("PC2: ",percentVark[2],"% variance")) +
  coord_fixed()

print("Results")
#Comparações que serão feitas
res<-results(ddsk, contrast=c("type","treated","untreated"),  filterFun=ihw)

#Adiciona anotação para os genes
#res$ensembl_id <- rownames(res)
#martk <- biomaRt::useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
#genes.tablek <- biomaRt::getBM(attributes = c("ensembl_gene_id_version", "external_gene_name", "description", "gene_biotype"), mart = martk)
#res <- merge(x = as.matrix(res), y = genes.tablek, by.x = "ensembl_id", by.y = "ensembl_gene_id_version", all.x = T, all.y = F )
counts_normalizedk<-counts(ddsk, normalized=TRUE)
counts_rawk<-counts(ddsk)
Tabelao<-cbind(res,counts_rawk,counts_normalizedk)

write.table(as.data.frame(Tabelao),file=snakemake@output[[1]], sep = "\t")
