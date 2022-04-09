# Developing Volcano Plot trial 
## Workflow developed from: https://gist.github.com/stephenturner/f60c1934405c127f09a6 
#### and http://www.bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#countmat

#Export gene expression matrix for each group  ==== 

dgcMatrix <- crhb_only@assays$RNA@data
pc <- as.data.frame(dgcMatrix)
pc<- data.frame(dgcMatrix)
write.csv(pc, "p+c+_DESEQ_Input.csv")

dgcMatrix <- pdyn_only@assays$RNA@data
p <- as.data.frame(dgcMatrix)
p<- data.frame(dgcMatrix)
write.csv(p, "p+c-a-_DESEQ_Input.csv")

dgcMatrix <- crhb_avp@assays$RNA@data
pca <- as.data.frame(dgcMatrix)
pca <- data.frame(dgcMatrix)
write.csv(pca, "p+c+a+_DESEQ_Input.csv")

dgcMatrix <- avp_only@assays$RNA@data
pa <- as.data.frame(dgcMatrix)
pa <- data.frame(dgcMatrix)
write.csv(pa, "p+c-a+_DESEQ_Input.csv")

##Hidden Steps: developing DESEQ input files (gene expression matrix & df identifying groups)=====
#I renamed the cells in these files "pca", "pc", "pa", and "p" based on their respective expression with sequential numbers after their identity. (Ex.p, p.1, p.2, p.3...etc)
#then made a new csv file with all cell gene expression combined with gene names as rows and columns as cell names (pca,pc,pa,p)
#After I copied + transposed the column names (cell names) as rownames into a new csv and made the second column their identity without the numbers (p, pca, pa, pc)

library("pasilla")

pasCts <- read.csv("pdyn_GeneExpMatrix.csv", sep = ",", row.names = NULL)
rownames(pasCts)<- make.names(pasCts$gene_id, unique = TRUE)
pasCts$gene_id <- NULL
cts <- as.matrix(pasCts, sep = ",")
#Make geneIDs unique (some are only different by capitalization) and transfrom into a matrix
pasAnno <- read.csv("Conditions_pdyn.csv", row.names = 1)
coldata <- pasAnno
#read in group identity df 
coldata$Genotype <- factor(coldata$Genotype)
#creates the genotypes (pca, p, pa, pc) as factors 

head(cts,2) #check if everything is outputting correctly
coldata

all(rownames(coldata) %in% colnames(cts)) #TRUE
cts <- cts[, rownames(coldata)]
all(rownames(coldata) == colnames(cts)) #TRUE
#CRITICAL: Makes sure rownames are in same order in both files!!

library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = round(cts),
                              colData = coldata,
                              design = ~ Genotype)
#Stores the input values, calculations, and results from DE analysis
#Uses a formula which expressed how the counts of each gene depends on Genotype 
dds <- DESeq(dds) # takes a while ~10 min
#performs default analysis by: 1. estimating size factors 2. estimation of dispersion factors 3. negative binomial GLM fitting and Wald statistic

# Plot dispersions
png("qc-dispersions.png", 1000, 1000, pointsize=20)
plotDispEsts(dds, main="Dispersion plot")
dev.off() 

# Regularized log transformation 
rld <- varianceStabilizingTransformation(dds) #takes very long 
head(assay(rld))
#calculates a variance stabilizing transformation (VST) from the fitted dispersion-mean relation(s) and then transforms the count data (normalized by division by the size factors or normalization factors)

# Get differential expression results
res <- results(dds) #Used default settings: lfc threshold is 0, pvalue= 0.1
#extracts a result table from a DESeq analysis giving base means across samples, log2 fold changes, standard errors, test statistics, p-values and adjusted p-values; 
res <- res[order(res$padj), ]
## Order by adjusted p-value
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
## Merge with normalized count data
names(resdata)[1] <- "Gene" 
head(resdata)
write.csv(resdata, file="total_diffexpr-results.csv")
## Write results. This is the DE of all groups 

test <-resultsNames(dds)
test #outputs group comparison names 
res <- results(dds, name = "Genotype_pca_vs_p")
write.csv(res, "pca_vs_p_DE.csv") 
#outputs DE results 

EnhancedVolcano(res,title = 'p+c+a+/p+c-a-',
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue', 
                pCutoff = 10e-02,
                FCcutoff = 1)

res <- results(dds, name = "Genotype_pc_vs_p")
write.csv(res,"pc_vs_p_DE.csv")

EnhancedVolcano(res,title = 'p+c+a-/p+c-a-',
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue', 
                pCutoff = 10e-02,
                FCcutoff = 1)

res <- results(dds, name = "Genotype_pa_vs_p")
write.csv(res, "pa_vs_p_DE.csv")

library(EnhancedVolcano)
EnhancedVolcano(res,title = 'p+c-a+/p+c-a-',
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue', 
                pCutoff = 10e-02,
                FCcutoff = 1)

