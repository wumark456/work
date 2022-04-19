dir.create('PDFs')
dir.create('files')

rm(list = ls())
library(affy)
library(limma)
library(stringr)
library(dplyr)
source('d:/git/ProjCodes/source.R')



# 注释文件 ###########
library(stringr)
gpl570 <- read.delim('raw_dat/osteoporosis/GSE7158_family/GPL570-tbl-1.txt', header = F)
gpl570_filtered <- gpl570[which(str_split_fixed(gpl570$V11, ' /// ', 2)[, 2] == ''), ]
gpl570_filtered <- gpl570_filtered[gpl570_filtered$V11 != '', ]

length(gpl570_filtered$V11)
length(unique(gpl570_filtered$V11))


# 1、整理数据集信息 ########
# 1.1、基因类型整理 #########
genetag <- read.delim('e:/public/GeneTag.genecode.v32.txt', 
                      header = T, stringsAsFactors = F)
genetag <- genetag[!duplicated(genetag$ENSGID), ]
rownames(genetag) <- genetag$ENSGID
table(genetag$TYPE)
protein_genes <- genetag[genetag$TYPE == 'protein_coding', ]$SYMBOL
protein_genes

# 2、GSE35959 差异基因的鉴定分析 ###############
# 2.1、GSE35959 数据处理 ###########
GSE35959_cli <- getGEOcli('GSE35959')
GSE35959_cli1 <- GSE35959_cli[, c("Acc", "Title", "age")]
colnames(GSE35959_cli1) <- c('Samples', 'Type', 'Age')
GSE35959_cli1$Type <- str_split_fixed(GSE35959_cli1$Type, '_', 7)[, 1]
table(GSE35959_cli1$Type)
GSE35959_cli1$Type <- ifelse(GSE35959_cli1$Type == 'hMSC-osteopo', 'OP', 'Control')
GSE35959_cli1$GSE <- 'GSE35959'
rownames(GSE35959_cli1) <- GSE35959_cli1$Samples
GSE35959_cli1$Age <- as.numeric(str_split_fixed(GSE35959_cli1$Age, ' ', 3)[, 1])



GSE35959_raw = ReadAffy(celfile.path="raw_dat/osteoporosis/GSE35959_RAW")
GSE35959_eset <- rma(GSE35959_raw)
GSE35959_eset <- exprs(GSE35959_eset)
dim(GSE35959_eset)
boxplot(GSE35959_eset,las=2)
colnames(GSE35959_eset)
colnames(GSE35959_eset) <- substr(colnames(GSE35959_eset), 1, 9)
class(GSE35959_eset)
dim(GSE35959_eset)
GSE35959_eset <- as.data.frame(GSE35959_eset)
boxplot(GSE35959_eset,las=2)
GSE35959_exp <- GSE35959_eset[, rownames(GSE35959_cli1)]
dim(GSE35959_exp)
boxplot(GSE35959_exp,las=2)
intersect(rownames(GSE35959_cli1), colnames(GSE35959_exp))

com_probes <- intersect(gpl570_filtered$V1, rownames(GSE35959_exp))

rownames(gpl570_filtered) <- gpl570_filtered$V1
gpl570_filtered <- gpl570_filtered[com_probes, ]
GSE35959_exp <- GSE35959_exp[com_probes, ]
GSE35959_exp$gene <- gpl570_filtered$V11
GSE35959_exp <- aggregate(.~gene, data = GSE35959_exp, median)
rownames(GSE35959_exp) <- GSE35959_exp$gene
GSE35959_exp <- GSE35959_exp[, -1]
class(GSE35959_exp)
boxplot(GSE35959_exp, las = 2)


# 2.2、差异基因的分析 ###############
table(GSE35959_cli1$Type)
GSE35959_limma <- DEGs_limma(GSE35959_exp[, GSE35959_cli1$Samples],
                             GSE35959_cli1[, c("Samples", "Type")],
                             'OP', 'Control')
dim(GSE35959_limma)
GSE35959_limma_filtered <- GSE35959_limma[GSE35959_limma$adj.P.Val < 0.05 & abs(GSE35959_limma$logFC) > 1, ]
dim(GSE35959_limma_filtered)

write.csv(GSE35959_limma_filtered,
          file = 'results/S1.csv')

GSE35959_limma_filtered <- GSE35959_limma_filtered[order(GSE35959_limma_filtered$logFC), ]

GSE35959_up <- rownames(GSE35959_limma_filtered[GSE35959_limma_filtered$logFC > 0, ])
GSE35959_dn <- rownames(GSE35959_limma_filtered[GSE35959_limma_filtered$logFC < 0, ])




# 2.2、火山图绘制 #################
wb_volcano <- function(gene, p.value, logFC,
                       cut_p = 0.05, cut_logFC = 1,
                       xlab = "log2FC", ylab = "-log10(FDR)",
                       title = 'State') {
  library(ggplot2)
  dat <- data.frame(Gene = gene, p.value = p.value, logFC = logFC)
  dat$type = ifelse(dat$p.value < cut_p & abs(dat$logFC) >= cut_logFC, ifelse(dat$logFC > cut_logFC, 'Up', 'Down'), 'Stable')
  
  p <- ggplot(dat, aes(x = logFC,
                       y = -log10(p.value))) +
    geom_point(aes(color = type),alpha=0.4, size=2) +
    scale_color_manual(values=c("#546de5", "#d2dae2","#ff4757"))+
    # 辅助线
    geom_vline(xintercept=c(-cut_logFC,cut_logFC),
               lty=4,col="black",lwd=0.8) +
    geom_hline(yintercept = -log10(cut_p),lty=4,col="black",lwd=0.8) +
    # 坐标轴
    labs(x=xlab, y=ylab) + xlim(-2, 2) +
    theme_bw()+
    # 图例
    theme(plot.title = element_text(hjust = 0.5),
          legend.position="top")
  p
}
pdf('PDFs/GSE35959_diffgene_volcano.pdf', width = 5, height = 5)
wb_volcano(rownames(GSE35959_limma),
           GSE35959_limma$adj.P.Val,
           GSE35959_limma$logFC,
           ylab = "-log10(FDR)",
           cut_p = 0.05,
           cut_logFC = 1)
dev.off()

# 2.3、热图绘制 ###################
library(pheatmap)
library(dplyr)
colnames(GSE35959_cli1)
GSE35959_cli1 <- arrange(GSE35959_cli1, Type)
annotation_col  <- data.frame(Type = factor(GSE35959_cli1$Type))
rownames(annotation_col) <- GSE35959_cli1$Samples
GSE35959_diff_exp_gene_heatmap <- GSE35959_exp[c(GSE35959_up,GSE35959_dn),
                                               GSE35959_cli1$Samples]
bk=unique(c(seq(-1.7, 1.7, length=100)))

pdf('PDFs/GSE35959_diff_exp_gene_heatmap.pdf',width = 6,height = 4)
pheatmap(GSE35959_diff_exp_gene_heatmap,
         scale = 'row',
         breaks = bk,
         annotation_col = annotation_col,
         cluster_cols = F, cluster_rows = F,
         show_rownames = F, show_colnames = F,
         gaps_col = 14, gaps_row = 152,
         # cellwidth = 0.8, cellheight = 0.35,
         # color = colorRampPalette(c("navy", "white", "firebrick3"))(100)
         # ,annotation_colors = ann_colors
)
dev.off()


# 5.1、差异基因功能分析 ######
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

diff_genes <- rownames(GSE35959_limma_filtered)
enrich_diff_genes <- bitr(diff_genes, fromType = "SYMBOL",
                          toType = c("ENTREZID", "ENSEMBL", "SYMBOL"),
                          OrgDb = org.Hs.eg.db)

enrich_diff_genes_BP <- enrichGO(gene          = enrich_diff_genes$ENTREZID,
                                 OrgDb         = org.Hs.eg.db,
                                 ont           = "BP",
                                 pAdjustMethod = "BH",
                                 pvalueCutoff  = 1,
                                 qvalueCutoff  = 1,
                                 readable      = TRUE)
head(enrich_diff_genes_BP)
diff_BP <- barplot(enrich_diff_genes_BP,
                   showCategory = 10,
                   color = "pvalue",
                   title = "diff Genes BP")
diff_BP

enrich_diff_genes_BP_dat <- as.data.frame(enrich_diff_genes_BP)


enrich_diff_genes_MF <- enrichGO(gene          = enrich_diff_genes$ENTREZID,
                                 OrgDb         = org.Hs.eg.db,
                                 ont           = "MF",
                                 pAdjustMethod = "BH",
                                 pvalueCutoff  = 1,
                                 qvalueCutoff  = 1,
                                 readable      = TRUE)
head(enrich_diff_genes_MF)
diff_MF <- barplot(enrich_diff_genes_MF,
                   showCategory = 10,
                   color = "pvalue",
                   title = "diff Genes MF")
diff_MF

enrich_diff_genes_MF_dat <- as.data.frame(enrich_diff_genes_MF)


enrich_diff_genes_CC <- enrichGO(gene          = enrich_diff_genes$ENTREZID,
                                 OrgDb         = org.Hs.eg.db,
                                 ont           = "CC",
                                 pAdjustMethod = "BH",
                                 pvalueCutoff  = 1,
                                 qvalueCutoff  = 1,
                                 readable      = TRUE)
head(enrich_diff_genes_CC)
diff_CC <- barplot(enrich_diff_genes_CC,
                   showCategory = 10,
                   color = "pvalue",
                   title = "diff Genes CC")
diff_CC

enrich_diff_genes_CC_dat <- as.data.frame(enrich_diff_genes_CC)

enrich_diff_genes_KEGG <- enrichKEGG(gene         = enrich_diff_genes$ENTREZID,
                                     organism     = 'hsa',
                                     pvalueCutoff = 1,
                                     pAdjustMethod = "BH")
head(enrich_diff_genes_KEGG)
diff_kegg <- barplot(enrich_diff_genes_KEGG,
                     showCategory = 10,
                     color = "pvalue",
                     title = "diff Genes KEGG")
diff_kegg
enrich_diff_genes_KEGG_dat <- setReadable(enrich_diff_genes_KEGG, 
                                          OrgDb = org.Hs.eg.db, 
                                          keyType="ENTREZID")
enrich_diff_genes_KEGG_dat <- as.data.frame(enrich_diff_genes_KEGG_dat)

colnames(enrich_diff_genes_BP_dat)
colnames(enrich_diff_genes_MF_dat)
colnames(enrich_diff_genes_CC_dat)
colnames(enrich_diff_genes_KEGG_dat)

enrich_diff_genes_BP_dat$TYPE <- 'BP'
enrich_diff_genes_MF_dat$TYPE <- 'MF'
enrich_diff_genes_CC_dat$TYPE <- 'CC'
enrich_diff_genes_KEGG_dat$TYPE <- 'KEGG'

enrich_diff_genes_dat <- rbind(enrich_diff_genes_BP_dat,
                               enrich_diff_genes_MF_dat,
                               enrich_diff_genes_CC_dat,
                               enrich_diff_genes_KEGG_dat)
enrich_diff_genes_dat_filtered <- enrich_diff_genes_dat[enrich_diff_genes_dat$pvalue < 0.05, ]
table(enrich_diff_genes_dat_filtered$TYPE)

write.table(enrich_diff_genes_dat_filtered,
            file = 'results/S2.txt',
            sep = '\t',
            quote = F,
            row.names = F)


diff_go_kegg <- plot_grid(diff_BP, 
                          diff_MF,
                          diff_CC,
                          diff_kegg, 
                          ncol=2, 
                          labels = LETTERS[1:4],
                          align = 'hv')

diff_go_kegg
library(ggplot2)
ggsave(plot = diff_go_kegg,
       filename = 'PDFs/diff_go_kegg.pdf',
       width = 18.5,height = 10)



# WGCNA 分析 ####################
# 4、WGCNA 共表达分析 #############
library(WGCNA)
GSE35959_WGCNA_cli <- GSE35959_cli1
rownames(GSE35959_WGCNA_cli) <- GSE35959_WGCNA_cli$Samples

dim(GSE35959_exp)
datExpr <- as.data.frame(t(GSE35959_exp))
gsg <- goodSamplesGenes(datExpr, verbose = 3)
gsg$allOK

sampleTree = hclust(dist(datExpr), method = "ward.D2")
pdf('PDFs/WGCNA-samples_cluster-1.pdf',width = 6,height = 6)
plot(sampleTree, main = "Sample clustering to detect outliers"
     , sub="", xlab="")
abline(h=180,col='red')
dev.off()
clust = cutreeStatic(sampleTree, cutHeight = 180, minSize = 3)
table(clust)
keepSamples = (clust==1)
datExpr = datExpr[keepSamples, ]
GSE35959_WGCNA_cli <- GSE35959_WGCNA_cli[rownames(datExpr),]
dim(datExpr)

save(datExpr, file = "datExpr.RData")

powers = c(1:20)
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
cex1 = 0.85

pdf('PDFs/WGCNA-samples_cluster-2.pdf',width = 7,height = 6)
par(mfrow = c(1,2))
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=cex1,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

softPower=sft$powerEstimate



net = blockwiseModules(datExpr, power = softPower, maxBlockSize = 22000,
                       TOMType = "unsigned", minModuleSize = 100,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       # saveTOMs = TRUE,
                       # saveTOMFileBase = "TPM-TOM-40", deepSplit = 3,
                       verbose = 3)
save(net, file = 'wgcna-net.RData')
table(net$colors)
mergedColors = labels2colors(net$colors)
table(mergedColors)
pdf('PDFs/WGCNA-samples_cluster-3.pdf',width = 6,height = 6)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    groupLabels = c("Module colors","GS.weight"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()


table(GSE35959_WGCNA_cli$Type)
GSE35959_WGCNA_cli$Samples <- rownames(GSE35959_WGCNA_cli)
GSE35959_WGCNA_cli$Control <- ifelse(GSE35959_WGCNA_cli$Type == 'Control', 1, 0)
GSE35959_WGCNA_cli$OP <- ifelse(GSE35959_WGCNA_cli$Type == 'OP', 1, 0)
MEs_col = net$MEs[GSE35959_WGCNA_cli$Samples, ]
colnames(MEs_col) = labels2colors(as.numeric(gsub('ME','',colnames(net$MEs))))



spms <- GSE35959_WGCNA_cli[, c("Control", "OP", "Age")]

# write.table(cbind(GSE35959_WGCNA_cli$Samples,spms),'files/spms.txt',sep = '\t',quote = F,row.names = F)
dim(spms)
modTraitCor = cor(MEs_col, as.data.frame(spms), use = "p")
modTraitP = corPvalueStudent(modTraitCor, dim(spms)[1])
textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)
pdf('PDFs/WGCNA-samples_cluster-4.pdf',width = 6,height = 6)
labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(spms), 
               yLabels = colnames(MEs_col), 
               cex.lab = 0.5, 
               ySymbols = colnames(MEs_col), colorLabels = FALSE, 
               colors = blueWhiteRed(50), 
               textMatrix = textMatrix, setStdMargins = FALSE, 
               cex.text = 0.5, zlim = c(-1,1), xLabelsAngle = 0,
               main = paste("Module-trait relationships"))
dev.off()



MEs = moduleEigengenes(datExpr[GSE35959_WGCNA_cli$Samples,], mergedColors)$eigengenes
MEs = orderMEs(MEs)
colnames(MEs)=gsub('^ME','',colnames(MEs))
me.inds=which(colnames(MEs)!='grey')
if(length(me.inds)==1){
  cns=colnames(MEs)[me.inds]
  MEs=as.data.frame(MEs[,me.inds])
  colnames(MEs)=cns
}else if(length(me.inds)>1){
  MEs=MEs[,me.inds]
}


nGenes = ncol(datExpr)
nSamples = nrow(datExpr[GSE35959_WGCNA_cli$Samples,])
geneTraitSignificance = as.data.frame(cor(datExpr[GSE35959_WGCNA_cli$Samples,],spms, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nrow(MEs)));
names(geneTraitSignificance) = paste("GS.", colnames(spms), sep="");
names(GSPvalue) = paste("p.GS.", colnames(spms), sep="")
geneModuleMembership = as.data.frame(cor(datExpr[GSE35959_WGCNA_cli$Samples,], MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", names(geneModuleMembership), sep="");
names(MMPvalue) = paste("p.MM", names(geneModuleMembership), sep="");
write.table(cbind(Genes=row.names(geneTraitSignificance),geneTraitSignificance,GSPvalue,Module=mergedColors)
            ,file = "GeneTraitSignificance.txt",sep = '\t',quote = F,row.names = F)


# 5、差异基因与 WGCNA 共表达基因交集以及功能分析 #################
# 5.1、green 模块基因的功能富集分析 ########
green_genes <- colnames(datExpr)[which(mergedColors=='green')]
enrich_green_genes <- bitr(green_genes, fromType = "SYMBOL",
                           toType = c("ENTREZID", "ENSEMBL", "SYMBOL"),
                           OrgDb = org.Hs.eg.db)

enrich_green_genes_BP <- enrichGO(gene          = enrich_green_genes$ENTREZID,
                                  OrgDb         = org.Hs.eg.db,
                                  ont           = "BP",
                                  pAdjustMethod = "BH",
                                  pvalueCutoff  = 1,
                                  qvalueCutoff  = 1,
                                  readable      = TRUE)
head(enrich_green_genes_BP)
green_BP <- barplot(enrich_green_genes_BP,
                    showCategory = 10,
                    color = "pvalue",
                    title = "green Genes BP")
green_BP

enrich_green_genes_BP_dat <- as.data.frame(enrich_green_genes_BP)


enrich_green_genes_MF <- enrichGO(gene          = enrich_green_genes$ENTREZID,
                                  OrgDb         = org.Hs.eg.db,
                                  ont           = "MF",
                                  pAdjustMethod = "BH",
                                  pvalueCutoff  = 1,
                                  qvalueCutoff  = 1,
                                  readable      = TRUE)
head(enrich_green_genes_MF)
green_MF <- barplot(enrich_green_genes_MF,
                    showCategory = 10,
                    color = "pvalue",
                    title = "green Genes MF")
green_MF

enrich_green_genes_MF_dat <- as.data.frame(enrich_green_genes_MF)


enrich_green_genes_CC <- enrichGO(gene          = enrich_green_genes$ENTREZID,
                                  OrgDb         = org.Hs.eg.db,
                                  ont           = "CC",
                                  pAdjustMethod = "BH",
                                  pvalueCutoff  = 1,
                                  qvalueCutoff  = 1,
                                  readable      = TRUE)
head(enrich_green_genes_CC)
green_CC <- barplot(enrich_green_genes_CC,
                    showCategory = 10,
                    color = "pvalue",
                    title = "green Genes CC")
green_CC

enrich_green_genes_CC_dat <- as.data.frame(enrich_green_genes_CC)

enrich_green_genes_KEGG <- enrichKEGG(gene         = enrich_green_genes$ENTREZID,
                                      organism     = 'hsa',
                                      pvalueCutoff = 1,
                                      pAdjustMethod = "BH")
head(enrich_green_genes_KEGG)
green_kegg <- barplot(enrich_green_genes_KEGG,
                      showCategory = 10,
                      color = "pvalue",
                      title = "green Genes KEGG")
green_kegg
enrich_green_genes_KEGG_dat <- setReadable(enrich_green_genes_KEGG, 
                                           OrgDb = org.Hs.eg.db, 
                                           keyType="ENTREZID")
enrich_green_genes_KEGG_dat <- as.data.frame(enrich_green_genes_KEGG_dat)
enrich_green_genes_KEGG_dat

colnames(enrich_green_genes_BP_dat)
colnames(enrich_green_genes_MF_dat)
colnames(enrich_green_genes_CC_dat)
colnames(enrich_green_genes_KEGG_dat)

enrich_green_genes_BP_dat$TYPE <- 'BP'
enrich_green_genes_MF_dat$TYPE <- 'MF'
enrich_green_genes_CC_dat$TYPE <- 'CC'
enrich_green_genes_KEGG_dat$TYPE <- 'KEGG'

enrich_green_genes_dat <- rbind(enrich_green_genes_BP_dat,
                                enrich_green_genes_MF_dat,
                                enrich_green_genes_CC_dat,
                                enrich_green_genes_KEGG_dat)
enrich_green_genes_dat_filtered <- enrich_green_genes_dat[enrich_green_genes_dat$pvalue < 0.05, ]
table(enrich_green_genes_dat_filtered$TYPE)

write.table(enrich_green_genes_dat_filtered,
            file = 'results/S4.txt',
            sep = '\t', quote = F, row.names = F)

green_go_kegg <- plot_grid(green_BP, 
                           green_MF,
                           green_CC,
                           green_kegg, 
                           ncol=2, 
                           labels = LETTERS[1:4],
                           align = 'hv')

green_go_kegg
library(ggplot2)
ggsave(plot = green_go_kegg,
       filename = 'PDFs/green_go_kegg.pdf',
       width = 18.5,height = 10)

# 5.2、pink 模块基因的功能富集分析 ########
pink_genes <- colnames(datExpr)[which(mergedColors=='pink')]
enrich_pink_genes <- bitr(pink_genes, fromType = "SYMBOL",
                          toType = c("ENTREZID", "ENSEMBL", "SYMBOL"),
                          OrgDb = org.Hs.eg.db)

enrich_pink_genes_BP <- enrichGO(gene          = enrich_pink_genes$ENTREZID,
                                 OrgDb         = org.Hs.eg.db,
                                 ont           = "BP",
                                 pAdjustMethod = "BH",
                                 pvalueCutoff  = 1,
                                 qvalueCutoff  = 1,
                                 readable      = TRUE)
head(enrich_pink_genes_BP)
pink_BP <- barplot(enrich_pink_genes_BP,
                   showCategory = 10,
                   color = "pvalue",
                   title = "pink Genes BP")
pink_BP

enrich_pink_genes_BP_dat <- as.data.frame(enrich_pink_genes_BP)


enrich_pink_genes_MF <- enrichGO(gene          = enrich_pink_genes$ENTREZID,
                                 OrgDb         = org.Hs.eg.db,
                                 ont           = "MF",
                                 pAdjustMethod = "BH",
                                 pvalueCutoff  = 1,
                                 qvalueCutoff  = 1,
                                 readable      = TRUE)
head(enrich_pink_genes_MF)
pink_MF <- barplot(enrich_pink_genes_MF,
                   showCategory = 10,
                   color = "pvalue",
                   title = "pink Genes MF")
pink_MF

enrich_pink_genes_MF_dat <- as.data.frame(enrich_pink_genes_MF)


enrich_pink_genes_CC <- enrichGO(gene          = enrich_pink_genes$ENTREZID,
                                 OrgDb         = org.Hs.eg.db,
                                 ont           = "CC",
                                 pAdjustMethod = "BH",
                                 pvalueCutoff  = 1,
                                 qvalueCutoff  = 1,
                                 readable      = TRUE)
head(enrich_pink_genes_CC)
pink_CC <- barplot(enrich_pink_genes_CC,
                   showCategory = 10,
                   color = "pvalue",
                   title = "pink Genes CC")
pink_CC

enrich_pink_genes_CC_dat <- as.data.frame(enrich_pink_genes_CC)

enrich_pink_genes_KEGG <- enrichKEGG(gene         = enrich_pink_genes$ENTREZID,
                                     organism     = 'hsa',
                                     pvalueCutoff = 1,
                                     pAdjustMethod = "BH")
head(enrich_pink_genes_KEGG)
pink_kegg <- barplot(enrich_pink_genes_KEGG,
                     showCategory = 10,
                     color = "pvalue",
                     title = "pink Genes KEGG")
pink_kegg
enrich_pink_genes_KEGG_dat <- setReadable(enrich_pink_genes_KEGG, 
                                          OrgDb = org.Hs.eg.db, 
                                          keyType="ENTREZID")
enrich_pink_genes_KEGG_dat <- as.data.frame(enrich_pink_genes_KEGG_dat)
enrich_pink_genes_KEGG_dat

colnames(enrich_pink_genes_BP_dat)
colnames(enrich_pink_genes_MF_dat)
colnames(enrich_pink_genes_CC_dat)
colnames(enrich_pink_genes_KEGG_dat)

enrich_pink_genes_BP_dat$TYPE <- 'BP'
enrich_pink_genes_MF_dat$TYPE <- 'MF'
enrich_pink_genes_CC_dat$TYPE <- 'CC'
enrich_pink_genes_KEGG_dat$TYPE <- 'KEGG'

enrich_pink_genes_dat <- rbind(enrich_pink_genes_BP_dat,
                               enrich_pink_genes_MF_dat,
                               enrich_pink_genes_CC_dat,
                               enrich_pink_genes_KEGG_dat)
enrich_pink_genes_dat_filtered <- enrich_pink_genes_dat[enrich_pink_genes_dat$pvalue < 0.05, ]
table(enrich_pink_genes_dat_filtered$TYPE)

write.table(enrich_pink_genes_dat_filtered,
            file = 'results/S5.txt',
            sep = '\t', quote = F, row.names = F)

pink_go_kegg <- plot_grid(pink_BP, 
                          pink_MF,
                          pink_CC,
                          pink_kegg, 
                          ncol=2, 
                          labels = LETTERS[1:4],
                          align = 'hv')

pink_go_kegg
library(ggplot2)
ggsave(plot = pink_go_kegg,
       filename = 'PDFs/pink_go_kegg.pdf',
       width = 18.5,height = 10)
# 5.3、tan 模块基因的功能富集分析 ########
tan_genes <- colnames(datExpr)[which(mergedColors=='tan')]
enrich_tan_genes <- bitr(tan_genes, fromType = "SYMBOL",
                         toType = c("ENTREZID", "ENSEMBL", "SYMBOL"),
                         OrgDb = org.Hs.eg.db)

enrich_tan_genes_BP <- enrichGO(gene          = enrich_tan_genes$ENTREZID,
                                OrgDb         = org.Hs.eg.db,
                                ont           = "BP",
                                pAdjustMethod = "BH",
                                pvalueCutoff  = 1,
                                qvalueCutoff  = 1,
                                readable      = TRUE)
head(enrich_tan_genes_BP)
tan_BP <- barplot(enrich_tan_genes_BP,
                  showCategory = 10,
                  color = "pvalue",
                  title = "tan Genes BP")
tan_BP

enrich_tan_genes_BP_dat <- as.data.frame(enrich_tan_genes_BP)


enrich_tan_genes_MF <- enrichGO(gene          = enrich_tan_genes$ENTREZID,
                                OrgDb         = org.Hs.eg.db,
                                ont           = "MF",
                                pAdjustMethod = "BH",
                                pvalueCutoff  = 1,
                                qvalueCutoff  = 1,
                                readable      = TRUE)
head(enrich_tan_genes_MF)
tan_MF <- barplot(enrich_tan_genes_MF,
                  showCategory = 10,
                  color = "pvalue",
                  title = "tan Genes MF")
tan_MF

enrich_tan_genes_MF_dat <- as.data.frame(enrich_tan_genes_MF)


enrich_tan_genes_CC <- enrichGO(gene          = enrich_tan_genes$ENTREZID,
                                OrgDb         = org.Hs.eg.db,
                                ont           = "CC",
                                pAdjustMethod = "BH",
                                pvalueCutoff  = 1,
                                qvalueCutoff  = 1,
                                readable      = TRUE)
head(enrich_tan_genes_CC)
tan_CC <- barplot(enrich_tan_genes_CC,
                  showCategory = 10,
                  color = "pvalue",
                  title = "tan Genes CC")
tan_CC

enrich_tan_genes_CC_dat <- as.data.frame(enrich_tan_genes_CC)

enrich_tan_genes_KEGG <- enrichKEGG(gene         = enrich_tan_genes$ENTREZID,
                                    organism     = 'hsa',
                                    pvalueCutoff = 1,
                                    pAdjustMethod = "BH")
head(enrich_tan_genes_KEGG)
tan_kegg <- barplot(enrich_tan_genes_KEGG,
                    showCategory = 10,
                    color = "pvalue",
                    title = "tan Genes KEGG")
tan_kegg
enrich_tan_genes_KEGG_dat <- setReadable(enrich_tan_genes_KEGG, 
                                         OrgDb = org.Hs.eg.db, 
                                         keyType="ENTREZID")
enrich_tan_genes_KEGG_dat <- as.data.frame(enrich_tan_genes_KEGG)

colnames(enrich_tan_genes_BP_dat)
colnames(enrich_tan_genes_MF_dat)
colnames(enrich_tan_genes_CC_dat)
colnames(enrich_tan_genes_KEGG_dat)

enrich_tan_genes_BP_dat$TYPE <- 'BP'
enrich_tan_genes_MF_dat$TYPE <- 'MF'
enrich_tan_genes_CC_dat$TYPE <- 'CC'
enrich_tan_genes_KEGG_dat$TYPE <- 'KEGG'

enrich_tan_genes_dat <- rbind(enrich_tan_genes_BP_dat,
                              enrich_tan_genes_MF_dat,
                              enrich_tan_genes_CC_dat,
                              enrich_tan_genes_KEGG_dat)
enrich_tan_genes_dat_filtered <- enrich_tan_genes_dat[enrich_tan_genes_dat$pvalue < 0.05, ]
table(enrich_tan_genes_dat_filtered$TYPE)

write.table(enrich_tan_genes_dat_filtered,
            file = 'results/S6.txt',
            row.names = F, sep = '\t', quote = F)


tan_go_kegg <- plot_grid(tan_BP, 
                         tan_MF,
                         tan_CC,
                         # tan_kegg, 
                         ncol=2, 
                         labels = LETTERS[1:3],
                         align = 'hv')

tan_go_kegg
library(ggplot2)
ggsave(plot = tan_go_kegg,
       filename = 'PDFs/tan_go_kegg.pdf',
       width = 15,height = 10)


library(VennDiagram)
venn.plot <- venn.diagram(x = list(diffgenes = rownames(GSE35959_limma_filtered),
                                   green_genes = green_genes,
                                   pink_genes = pink_genes,
                                   tan_genes = tan_genes),
                          filename = NULL,
                          fill = color8[1:4], # 可以自定义颜色
                          alpha = 0.5,
                          euler.d = TRUE,
                          cex = 1,
                          cat.cex = 1,
                          cat.pos = 0,
                          fontfamily = "serif",
                          cat.fontfamily = "serif")
grid.draw(venn.plot)
dev.off()

library(venn)
library(ggplot2)
pdf('PDFs/venn.pdf', width = 5, height = 5)
venn(list(Diff_genes = rownames(GSE35959_limma_filtered),
          Green_genes = green_genes,
          Pink_genes = pink_genes,
          Tan_genes = tan_genes),
     ilabels = TRUE,
     # ellipse = TRUE,
     zcolor = color8[1:4],
     box = FALSE)
dev.off()


com_genes <- intersect(rownames(GSE35959_limma_filtered),
                       c(green_genes, pink_genes, tan_genes))

write.table(com_genes,
            file = 'results/S7.txt',
            row.names = F, col.names = F,
            quote = F)



# 5.1、共有基因模块基因功能分析 ######
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
enrich_com_genes <- bitr(com_genes, fromType = "SYMBOL",
                         toType = c("ENTREZID", "ENSEMBL", "SYMBOL"),
                         OrgDb = org.Hs.eg.db)

enrich_com_genes_BP <- enrichGO(gene          = enrich_com_genes$ENTREZID,
                                OrgDb         = org.Hs.eg.db,
                                ont           = "BP",
                                pAdjustMethod = "BH",
                                pvalueCutoff  = 1,
                                qvalueCutoff  = 1,
                                readable      = TRUE)
head(enrich_com_genes_BP)
com_BP <- barplot(enrich_com_genes_BP,
                  showCategory = 10,
                  color = "pvalue",
                  title = "com Genes BP")
com_BP

enrich_com_genes_BP_dat <- as.data.frame(enrich_com_genes_BP)


enrich_com_genes_MF <- enrichGO(gene          = enrich_com_genes$ENTREZID,
                                OrgDb         = org.Hs.eg.db,
                                ont           = "MF",
                                pAdjustMethod = "BH",
                                pvalueCutoff  = 1,
                                qvalueCutoff  = 1,
                                readable      = TRUE)
head(enrich_com_genes_MF)
com_MF <- barplot(enrich_com_genes_MF,
                  showCategory = 10,
                  color = "pvalue",
                  title = "com Genes MF")
com_MF

enrich_com_genes_MF_dat <- as.data.frame(enrich_com_genes_MF)


enrich_com_genes_CC <- enrichGO(gene          = enrich_com_genes$ENTREZID,
                                OrgDb         = org.Hs.eg.db,
                                ont           = "CC",
                                pAdjustMethod = "BH",
                                pvalueCutoff  = 1,
                                qvalueCutoff  = 1,
                                readable      = TRUE)
head(enrich_com_genes_CC)
com_CC <- barplot(enrich_com_genes_CC,
                  showCategory = 10,
                  color = "pvalue",
                  title = "com Genes CC")
com_CC

enrich_com_genes_CC_dat <- as.data.frame(enrich_com_genes_CC)

enrich_com_genes_KEGG <- enrichKEGG(gene         = enrich_com_genes$ENTREZID,
                                    organism     = 'hsa',
                                    pvalueCutoff = 1,
                                    pAdjustMethod = "BH")
head(enrich_com_genes_KEGG)
com_kegg <- barplot(enrich_com_genes_KEGG,
                    showCategory = 10,
                    color = "pvalue",
                    title = "com Genes KEGG")
com_kegg
enrich_com_genes_KEGG_dat <- setReadable(enrich_com_genes_KEGG, 
                                         OrgDb = org.Hs.eg.db, 
                                         keyType="ENTREZID")
enrich_com_genes_KEGG_dat <- as.data.frame(enrich_com_genes_KEGG_dat)
enrich_com_genes_KEGG_dat

colnames(enrich_com_genes_BP_dat)
colnames(enrich_com_genes_MF_dat)
colnames(enrich_com_genes_CC_dat)
colnames(enrich_com_genes_KEGG_dat)

enrich_com_genes_BP_dat$TYPE <- 'BP'
enrich_com_genes_MF_dat$TYPE <- 'MF'
enrich_com_genes_CC_dat$TYPE <- 'CC'
enrich_com_genes_KEGG_dat$TYPE <- 'KEGG'

enrich_com_genes_dat <- rbind(enrich_com_genes_BP_dat,
                              enrich_com_genes_MF_dat,
                              enrich_com_genes_CC_dat,
                              enrich_com_genes_KEGG_dat)
enrich_com_genes_dat_filtered <- enrich_com_genes_dat[enrich_com_genes_dat$pvalue < 0.05, ]
table(enrich_com_genes_dat_filtered$TYPE)

write.table(enrich_com_genes_dat_filtered,
            file = 'results/S8.txt',
            row.names = F, sep = '\t', quote = F)

com_go_kegg <- plot_grid(com_BP, 
                         com_MF,
                         com_CC,
                         com_kegg,
                         ncol=2, 
                         labels = LETTERS[1:4],
                         align = 'hv')

com_go_kegg
library(ggplot2)
ggsave(plot = com_go_kegg,
       filename = 'PDFs/com_go_kegg.pdf',
       width = 18.5,height = 10)

save.image('OP_001.Rdata')

# STRING PPI 分析 ###############
dir.create('files/STRING')
write.table(com_genes,
            file = 'files/STRING/com_genes.txt',
            col.names = F, row.names = F, quote = F)

Betweenness <- read.csv('files/STRING/Betweenness_top15.csv')
Betweenness_genes <- Betweenness$name

Closeness <- read.csv('files/STRING/Closeness_top15.csv')
Closeness_genes <- Closeness$name

Degree <- read.csv('files/STRING/Degree_top15.csv')
Degree_genes <- Degree$name

MCC <- read.csv('files/STRING/MCC_top15.csv')
MCC_genes <- MCC$name


library(VennDiagram)
hub.venn.plot <- venn.diagram(x = list(Betweenness = Betweenness_genes,
                                       Closeness = Closeness_genes,
                                       Degree = Degree_genes,
                                       MCC = MCC_genes),
                              filename = NULL,
                              fill = color8[1:4], # 可以自定义颜色
                              alpha = 0.5,
                              euler.d = TRUE,
                              cex = 1,
                              cat.cex = 1,
                              cat.pos = 0,
                              fontfamily = "serif",
                              cat.fontfamily = "serif")
grid.draw(hub.venn.plot)
dev.off()

library(venn)
library(ggplot2)
pdf('PDFs/hubgenes_venn.pdf', width = 5, height = 5)
venn(list(Betweenness = Betweenness_genes,
          Closeness = Closeness_genes,
          Degree = Degree_genes,
          MCC = MCC_genes),
     ilabels = TRUE,
     # ellipse = TRUE,
     zcolor = color9[1:4],
     box = FALSE)
dev.off()

intersect(intersect(Betweenness_genes, Closeness_genes),
          intersect(Degree_genes, MCC_genes))

hub_genes <- intersect(intersect(Betweenness_genes, Closeness_genes),
                       intersect(Degree_genes, MCC_genes))



# 6.1、数据集 GSE35959 的 ROC 曲线 ##########
rownames(GSE35959_cli1) <- GSE35959_cli1$Samples
GSE35959_roc_dat <- cbind(GSE35959_cli1, 
                          t(GSE35959_exp[hub_genes, GSE35959_cli1$Samples]))

library(factoextra)
library(FactoMineR)
library(ggbiplot)
GSE35959.pca <- prcomp(GSE35959_roc_dat[, 5:10],  scale = TRUE)
fviz_pca_ind(GSE35959.pca,  axes = c(1, 2), label="none", 
             addEllipses = T, ellipse.type="norm", ellipse.level=0.95,
             habillage = GSE35959_roc_dat$Type, palette = "jco",
             mean.point=F, legend.title = "Type")

ggbiplot(GSE35959.pca, scale=1, groups = GSE35959_roc_dat$Type,
         ellipse = TRUE,ellipse.prob=0.5, circle = F,var.axes=F) +
  scale_color_discrete(name = 'Type') +
  # scale_fill_manual(values = c("#00468BFF", "#ED0000FF", "#42B540FF")) +
  xlim(-3, 3) + ylim(-3, 3) + theme_bw() +
  theme(legend.direction = 'horizontal', legend.position = 'top')


ggplotViolin_group(GSE35959_roc_dat[, hub_genes],
                   GSE35959_roc_dat$Type,
                   cols = c("#E64B35", "#4DBBD5"),
                   title = 'Cluster',
                   xlab = '',
                   ylab = 'Score',
                   angle = 60)


# GSE35959_roc_dat$Type <- as.numeric(factor(GSE35959_roc_dat$Type))
plot_ROC <- function(dat, genes, outdir = './', prefix = '') {
  dir.create(outdir, recursive = T)
  library(pROC)
  dat$Type1 <- as.numeric(factor(dat$Type))
  Combine <- glm(as.formula(paste0('Type1 ~ ', paste0(genes, collapse = " + "))), 
                 data = dat)
  Combine <- predict(Combine)
  dat$Combine <- as.numeric(Combine)
  
  for (ge in c(genes, 'Combine')) {
    p_roc <- roc(as.formula(paste0('Type1 ~ ', paste0(ge, collapse = " + "))), 
                 data = dat, auc=TRUE, ci=TRUE)
    p_roc_plot <- ggroc(p_roc, size = 1.2, legacy.axes = TRUE)
    p_roc_plot <- p_roc_plot + 
      theme_bw() + labs(title = ge) + 
      theme(plot.title = element_text(hjust = 0.5)) +
      annotate("text", x = .75, y = .25, 
               label = paste("AUC = ", round(p_roc$auc, 3), '(', 
                             round(p_roc$ci[1], 3), '-',  
                             round(p_roc$ci[3], 3), ')'),
               size = 5)
    ggsave(plot = p_roc_plot,
           filename = paste0(outdir, prefix, '_', ge, '_ROC.pdf'),
           width = 5, height = 5)
    ggsave(plot = p_roc_plot,
           filename = paste0(outdir, prefix, '_', ge, '_ROC.png'),
           width = 5, height = 5)
  }
  return(dat)
}
GSE35959_roc_dat1 <- plot_ROC(dat = GSE35959_roc_dat,
                              genes = hub_genes,
                              outdir = 'PDFs/GSE35959/',
                              prefix = 'GSE35959')
testdat <- GSE35959_roc_dat 
testdat$Type1 <- as.numeric(factor(GSE35959_roc_dat$Type))
Combine <- glm(as.formula(paste0('Type1 ~ ', paste0(hub_genes, collapse = " + "))), 
               data = testdat)
Combine <- predict(Combine)

library(e1071)
library(affy)
set.seed(12345)
GSE35959.tObj <- tune.svm(GSE35959_roc_dat[, hub_genes],
                          factor(GSE35959_roc_dat$Type)
                          ,type="C-classification"
                          ,kernel="radial"
                          ,probability = TRUE
                          , cost=c(0.001,0.01,0.1,1,5,10,100,1000),scale=T)
GSE35959.BestSvm<-GSE35959.tObj$best.model
GSE35959.pre_label <- predict(GSE35959.BestSvm, 
                              GSE35959_roc_dat[, hub_genes], 
                              probability = F,
                              decision.values=T)

# dev.off()
GSE35959_svm_res <- cal_auc(GSE35959.pre_label,factor(GSE35959_roc_dat$Type))
table(GSE35959.pre_label,factor(GSE35959_roc_dat$Type))

GSE35959.pre_dat <- read.csv('files/svm/GSE35959.svm.csv')
GSE35959.pre_dat <- merge(GSE35959_roc_dat, GSE35959.pre_dat,
                          by = 'Samples')
GSE35959.pre_dat$Type1 <- as.numeric(factor(GSE35959.pre_dat$Type))

GSE35959_roc <- roc(Type1 ~ Values, 
                    data = GSE35959.pre_dat, auc=TRUE, ci=TRUE)
GSE35959_roc_plot <- ggroc(p_roc, size = 1.2, legacy.axes = TRUE)
GSE35959_roc_plot <- GSE35959_roc_plot + 
  theme_bw() + labs(title = 'GSE35959') + 
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate("text", x = .75, y = .25, 
           label = paste("AUC = ", round(p_roc$auc, 3), '(', 
                         round(p_roc$ci[1], 3), '-',  
                         round(p_roc$ci[3], 3), ')'),
           size = 5)
GSE35959_roc_plot
ggsave(plot = GSE35959_roc_plot,
       filename = 'PDFs/GSE35959_roc_plot.pdf',
       width = 5, height = 5)


# 6.2、数据集 GSE7429 的 ROC 曲线 ##########
library(stringr)
GSE7429_cli <- read.delim('raw_dat/osteoporosis/GSE7429_family/SampleInfo_contrib1-GPL96.txt')
GSE7429_cli1 <- GSE7429_cli[, c("Accession", "Title")]
colnames(GSE7429_cli1) <- c('Samples', 'Type')
GSE7429_cli1$Type <- str_split_fixed(GSE7429_cli1$Type, '-', 3)[, 2]
table(GSE7429_cli1$Type)
# GSE7429_cli1 <- GSE7429_cli1[GSE7429_cli1$Type != '', ]
GSE7429_cli1$Type <- ifelse(GSE7429_cli1$Type == 'lowBMD', 'OP', 'Control')
rownames(GSE7429_cli1) <- GSE7429_cli1$Samples

GSE7429_exp <- read.delim('raw_dat/osteoporosis/GSE7429_family/MergeExpro_contrib1-GPL96.txt',
                          row.names = 1)
GSE7429_exp <- GSE7429_exp[which(str_split_fixed(rownames(GSE7429_exp), ' /// ', 2)[, 2] == ''), ]
# GSE7429_exp <- log2(GSE7429_exp + 1)
boxplot(GSE7429_exp, las = 2)

GSE7429_genes <- intersect(hub_genes, rownames(GSE7429_exp))

rownames(GSE7429_cli1) <- GSE7429_cli1$Samples
GSE7429_roc_dat <- cbind(GSE7429_cli1, 
                         t(GSE7429_exp[GSE7429_genes, 
                                       GSE7429_cli1$Samples]))

GSE7429.pca <- prcomp(GSE7429_roc_dat[, 3:6],  scale = TRUE)
fviz_pca_ind(GSE7429.pca,  axes = c(1, 2), label="none", 
             addEllipses = T, ellipse.type="norm", ellipse.level=0.95,
             habillage = GSE7429_roc_dat$Type, palette = "jco",
             mean.point=F, legend.title = "Type")

ggbiplot(GSE7429.pca, scale=1, groups = GSE7429_roc_dat$Type,
         ellipse = TRUE,ellipse.prob=0.5, circle = F,var.axes=F) +
  scale_color_discrete(name = 'Type') +
  # scale_fill_manual(values = c("#00468BFF", "#ED0000FF", "#42B540FF")) +
  xlim(-3, 3) + ylim(-3, 3) + theme_bw() +
  theme(legend.direction = 'horizontal', legend.position = 'top')

ggplotViolin_group(GSE7429_roc_dat[, GSE7429_genes],
                   GSE7429_roc_dat$Type,
                   cols = c("#E64B35", "#4DBBD5"),
                   title = 'Cluster',
                   xlab = '',
                   ylab = 'Score',
                   angle = 60)

cal_auc=function(pre_label1,group){
  stat_res=table(pre_label1,group)
  acu=(stat_res[1, 1] + stat_res[2, 2])/length(pre_label1)
  sens=stat_res[1, 1]/(stat_res[1, 1] + stat_res[2, 1])
  spec=stat_res[2, 2]/(stat_res[1, 2] + stat_res[2, 2])
  plot(c(0,0,spec),c(0,sens,1),type='l',xlim = c(0,1),col='blue',xlab='1-Specificity',ylab='Sensitivity')
  # polygon(x=c(1,spec,0,0,1),y=c(0,sens,1,0,0),col="#66CCFF",border=NA)
  auc=spec+(1-spec)*sens/2-spec*(1-sens)/2
  print(stat_res)
  return(c(acu,sens,spec,auc))
}

set.seed(12345)
GSE7429.tObj <- tune.svm(GSE7429_roc_dat[, GSE7429_genes],
                         factor(GSE7429_roc_dat$Type)
                         ,type="C-classification"
                         ,kernel="radial"
                         ,probability = TRUE
                         , cost=c(0.001,0.01,0.1,1,5,10,100,1000),scale=T)
GSE7429.BestSvm<-GSE7429.tObj$best.model
GSE7429.pre_label <- predict(GSE7429.BestSvm, 
                             GSE7429_roc_dat[, GSE7429_genes], 
                             probability = F,
                             decision.values=T)
GSE7429.pre_label
cal_auc(GSE7429.pre_label,factor(GSE7429_roc_dat$Type))
table(GSE7429.pre_label,factor(GSE7429_roc_dat$Type))

GSE7429.pre_dat <- read.csv('files/svm/GSE7429.svm.csv')
GSE7429.pre_dat <- merge(GSE7429_roc_dat, GSE7429.pre_dat,
                         by = 'Samples')
GSE7429.pre_dat$Type1 <- as.numeric(factor(GSE7429.pre_dat$Type))

GSE7429_roc <- roc(Type1 ~ Values, 
                   data = GSE7429.pre_dat, auc=TRUE, ci=TRUE)
GSE7429_roc_plot <- ggroc(p_roc, size = 1.2, legacy.axes = TRUE)
GSE7429_roc_plot <- GSE7429_roc_plot + 
  theme_bw() + labs(title = 'GSE7429') + 
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate("text", x = .75, y = .25, 
           label = paste("AUC = ", round(p_roc$auc, 3), '(', 
                         round(p_roc$ci[1], 3), '-',  
                         round(p_roc$ci[3], 3), ')'),
           size = 5)
GSE7429_roc_plot
ggsave(plot = GSE7429_roc_plot,
       filename = 'PDFs/GSE7429_roc_plot.pdf',
       width = 5, height = 5)

# GSE7429_roc_dat$Type <- as.numeric(factor(GSE7429_roc_dat$Type))
plot_ROC(dat = GSE7429_roc_dat,
         genes = GSE7429_genes,
         outdir = 'PDFs/GSE7429/',
         prefix = 'GSE7429')




# 6.3、GSE56116 的 ROC 曲线 ###############
GSE56116_cli <- read.delim('raw_dat/osteoporosis/GSE56116_family/SampleInfo_contrib1-GPL4133.txt')
GSE56116_cli1 <- GSE56116_cli[, c("Accession", "X1.CLINICAL.PHENOTYPE")]
colnames(GSE56116_cli1) <- c('Samples', 'Type')
GSE56116_cli1$Type <- str_split_fixed(GSE56116_cli1$Type, ' ', 2)[, 2]
table(GSE56116_cli1$Type)
GSE56116_cli1$Type <- ifelse(GSE56116_cli1$Type == 'osteoporosis', 'OP', 'Control')
rownames(GSE56116_cli1) <- GSE56116_cli1$Samples

GSE56116_exp <- read.delim('raw_dat/osteoporosis/GSE56116_family/MergeExpro_contrib1-GPL4133.txt',
                           row.names = 1)
GSE56116_exp <- log2(GSE56116_exp + 1)
boxplot(GSE56116_exp, las = 2)


GSE56116_genes <- intersect(hub_genes, rownames(GSE56116_exp))

rownames(GSE56116_cli1) <- GSE56116_cli1$Samples
GSE56116_roc_dat <- cbind(GSE56116_cli1, 
                          t(GSE56116_exp[GSE56116_genes, 
                                         GSE56116_cli1$Samples]))

GSE56116.pca <- prcomp(GSE56116_roc_dat[, 3:8],  scale = TRUE)
fviz_pca_ind(GSE56116.pca,  axes = c(1, 2), label="none", 
             addEllipses = F, ellipse.type="norm", ellipse.level=0.95,
             habillage = GSE56116_roc_dat$Type, palette = "jco",
             mean.point=F, legend.title = "Type")

ggbiplot(GSE56116.pca, scale=1, groups = GSE56116_roc_dat$Type,
         ellipse = TRUE,ellipse.prob=0.5, circle = F,var.axes=F) +
  scale_color_discrete(name = 'Type') +
  # scale_fill_manual(values = c("#00468BFF", "#ED0000FF", "#42B540FF")) +
  xlim(-3, 3) + ylim(-3, 3) + theme_bw() +
  theme(legend.direction = 'horizontal', legend.position = 'top')

ggplotViolin_group(GSE56116_roc_dat[, GSE56116_genes],
                   GSE56116_roc_dat$Type,
                   cols = c("#E64B35", "#4DBBD5"),
                   title = 'Cluster',
                   xlab = '',
                   ylab = 'Score',
                   angle = 60)

set.seed(12345)
GSE56116.tObj <- tune.svm(GSE56116_roc_dat[, GSE56116_genes],
                          factor(GSE56116_roc_dat$Type)
                          ,type="C-classification"
                          ,kernel="radial"
                          ,probability = TRUE
                          , cost=c(0.001,0.01,0.1,1,5,10,100,1000),scale=T)
GSE56116.BestSvm<-GSE56116.tObj$best.model
GSE56116.pre_label <- predict(GSE56116.BestSvm, 
                              GSE56116_roc_dat[, GSE56116_genes], 
                              probability = F,
                              decision.values=T)
cal_auc(GSE56116.pre_label,factor(GSE56116_roc_dat$Type))
table(GSE56116.pre_label,factor(GSE56116_roc_dat$Type))

# GSE56116.pre_dat <- read.csv('files/svm/GSE56116.svm.csv')
# GSE56116.pre_dat <- merge(GSE56116_roc_dat, GSE56116.pre_dat,
#                           by = 'Samples')
# GSE56116.pre_dat$Type1 <- as.numeric(factor(GSE56116.pre_dat$Type))
# 
# GSE56116_roc <- roc(Type1 ~ Values, 
#                     data = GSE56116.pre_dat, auc=TRUE, ci=TRUE)
# GSE56116_roc_plot <- ggroc(p_roc, size = 1.2, legacy.axes = TRUE)
# GSE56116_roc_plot <- GSE56116_roc_plot + 
#   theme_bw() + labs(title = 'GSE56116') + 
#   theme(plot.title = element_text(hjust = 0.5)) +
#   annotate("text", x = .75, y = .25, 
#            label = paste("AUC = ", round(p_roc$auc, 3), '(', 
#                          round(p_roc$ci[1], 3), '-',  
#                          round(p_roc$ci[3], 3), ')'),
#            size = 5)
# GSE56116_roc_plot
# ggsave(plot = GSE56116_roc_plot,
#        filename = 'PDFs/GSE56116_roc_plot.pdf',
#        width = 5, height = 5)

GSE56116_roc_dat$Type <- as.numeric(factor(GSE56116_roc_dat$Type))
plot_ROC(dat = GSE56116_roc_dat,
         genes = GSE56116_genes,
         outdir = 'PDFs/GSE56116/',
         prefix = 'GSE56116')


# 6.4、GSE13850 的 ROC 曲线 ###############
GSE13850_cli <- read.delim('raw_dat/osteoporosis/GSE13850_family/SampleInfo_contrib1-GPL96.txt')
GSE13850_cli1 <- GSE13850_cli[, c("Accession", "Title")]
colnames(GSE13850_cli1) <- c('Samples', 'Type')
GSE13850_cli1$Type <- str_split_fixed(GSE13850_cli1$Type, '-', 3)[, 2]
table(GSE13850_cli1$Type)
GSE13850_cli1 <- GSE13850_cli1[GSE13850_cli1$Type != '', ]
GSE13850_cli1$Type <- ifelse(GSE13850_cli1$Type == 'lowBMD', 'OP', 'Control')
rownames(GSE13850_cli1) <- GSE13850_cli1$Samples

GSE13850_exp <- read.delim('raw_dat/osteoporosis/GSE13850_family/MergeExpro_contrib1-GPL96.txt',
                           row.names = 1)
GSE13850_exp <- GSE13850_exp[which(str_split_fixed(rownames(GSE13850_exp), ' /// ', 2)[, 2] == ''), ]
# GSE13850_exp <- log2(GSE13850_exp + 1)
boxplot(GSE13850_exp, las = 2)

GSE13850_genes <- intersect(hub_genes, rownames(GSE13850_exp))

rownames(GSE13850_cli1) <- GSE13850_cli1$Samples
GSE13850_roc_dat <- cbind(GSE13850_cli1, 
                          t(GSE13850_exp[GSE13850_genes, 
                                         GSE13850_cli1$Samples]))

GSE13850.pca <- prcomp(GSE13850_roc_dat[, 3:6],  scale = TRUE)
fviz_pca_ind(GSE13850.pca,  axes = c(1, 2), label="none", 
             addEllipses = F, ellipse.type="norm", ellipse.level=0.95,
             habillage = GSE13850_roc_dat$Type, palette = "jco",
             mean.point=F, legend.title = "Type")

ggbiplot(GSE13850.pca, scale=1, groups = GSE13850_roc_dat$Type,
         ellipse = TRUE,ellipse.prob=0.5, circle = F,var.axes=F) +
  scale_color_discrete(name = 'Type') +
  # scale_fill_manual(values = c("#00468BFF", "#ED0000FF", "#42B540FF")) +
  xlim(-3, 3) + ylim(-3, 3) + theme_bw() +
  theme(legend.direction = 'horizontal', legend.position = 'top')

ggplotViolin_group(GSE13850_roc_dat[, GSE13850_genes],
                   GSE13850_roc_dat$Type,
                   cols = c("#E64B35", "#4DBBD5"),
                   title = 'Cluster',
                   xlab = '',
                   ylab = 'Score',
                   angle = 60)


set.seed(12345)
GSE13850.tObj <- tune.svm(GSE13850_roc_dat[, GSE13850_genes],
                          factor(GSE13850_roc_dat$Type)
                          ,type="C-classification"
                          ,kernel="radial"
                          ,probability = TRUE
                          , cost=c(0.001,0.01,0.1,1,5,10,100,1000),scale=T)
GSE13850.BestSvm<-GSE13850.tObj$best.model
GSE13850.pre_label <- predict(GSE13850.BestSvm, 
                              GSE13850_roc_dat[, GSE13850_genes], 
                              probability = F,
                              decision.values=T)
cal_auc(GSE13850.pre_label,factor(GSE13850_roc_dat$Type))
table(GSE13850.pre_label,factor(GSE13850_roc_dat$Type))

GSE13850.pre_dat <- read.csv('files/svm/GSE13850.svm.csv')
GSE13850.pre_dat <- merge(GSE13850_roc_dat, GSE13850.pre_dat,
                          by = 'Samples')
GSE13850.pre_dat$Type1 <- as.numeric(factor(GSE13850.pre_dat$Type))

GSE13850_roc <- roc(Type1 ~ Values, 
                    data = GSE13850.pre_dat, auc=TRUE, ci=TRUE)
GSE13850_roc_plot <- ggroc(p_roc, size = 1.2, legacy.axes = TRUE)
GSE13850_roc_plot <- GSE13850_roc_plot + 
  theme_bw() + labs(title = 'GSE13850') + 
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate("text", x = .75, y = .25, 
           label = paste("AUC = ", round(p_roc$auc, 3), '(', 
                         round(p_roc$ci[1], 3), '-',  
                         round(p_roc$ci[3], 3), ')'),
           size = 5)
GSE13850_roc_plot
ggsave(plot = GSE13850_roc_plot,
       filename = 'PDFs/GSE13850_roc_plot.pdf',
       width = 5, height = 5)

GSE13850_roc_dat$Type <- as.numeric(factor(GSE13850_roc_dat$Type))
plot_ROC(dat = GSE13850_roc_dat,
         genes = GSE13850_genes,
         outdir = 'PDFs/GSE13850/',
         prefix = 'GSE13850')

# 6.5、GSE62402 的 ROC 曲线 ###############
GSE62402_cli <- read.delim('raw_dat/osteoporosis/GSE62402_family/SampleInfo_contrib1-GPL5175.txt')
View(GSE62402_cli)
GSE62402_cli1 <- GSE62402_cli[, c("Accession", "Title")]
colnames(GSE62402_cli1) <- c('Samples', 'Type')
GSE62402_cli1$Type <- str_split_fixed(GSE62402_cli1$Type, ' ', 4)[, 3]
table(GSE62402_cli1$Type)
GSE62402_cli1 <- GSE62402_cli1[GSE62402_cli1$Type != '', ]
GSE62402_cli1$Type <- ifelse(GSE62402_cli1$Type == 'low', 'OP', 'Control')
rownames(GSE62402_cli1) <- GSE62402_cli1$Samples

GSE62402_exp <- read.delim('raw_dat/osteoporosis/GSE62402_family/MergeExpro_contrib1-GPL5175.txt',
                           row.names = 1)
GSE62402_exp$genes <- str_split_fixed(rownames(GSE62402_exp), ' // ', 3)[, 2]
GSE62402_exp <- aggregate(.~genes, data = GSE62402_exp, FUN = median)

rownames(GSE62402_exp) <- GSE62402_exp$genes
GSE62402_exp <- GSE62402_exp[, -1]
# GSE62402_exp <- log2(GSE62402_exp + 1)
boxplot(GSE62402_exp, las = 2)

GSE62402_genes <- intersect(hub_genes, rownames(GSE62402_exp))

rownames(GSE62402_cli1) <- GSE62402_cli1$Samples
GSE62402_roc_dat <- cbind(GSE62402_cli1,
                          t(GSE62402_exp[GSE62402_genes,
                                         GSE62402_cli1$Samples]))

GSE62402.pca <- prcomp(GSE62402_roc_dat[, 3:8],  scale = TRUE)
fviz_pca_ind(GSE62402.pca,  axes = c(1, 2), label="none",
             addEllipses = F, ellipse.type="norm", ellipse.level=0.95,
             habillage = GSE62402_roc_dat$Type, palette = "jco",
             mean.point=F, legend.title = "Type")

ggbiplot(GSE62402.pca, scale=1, groups = GSE62402_roc_dat$Type,
         ellipse = TRUE,ellipse.prob=0.5, circle = F,var.axes=F) +
  scale_color_discrete(name = 'Type') +
  # scale_fill_manual(values = c("#00468BFF", "#ED0000FF", "#42B540FF")) +
  xlim(-3, 3) + ylim(-3, 3) + theme_bw() +
  theme(legend.direction = 'horizontal', legend.position = 'top')


ggplotViolin_group(GSE62402_roc_dat[, GSE62402_genes],
                   GSE62402_roc_dat$Type,
                   cols = c("#E64B35", "#4DBBD5"),
                   title = 'Cluster',
                   xlab = '',
                   ylab = 'Score',
                   angle = 60)

set.seed(12345)
GSE62402.tObj <- tune.svm(GSE62402_roc_dat[, GSE62402_genes],
                          factor(GSE62402_roc_dat$Type)
                          ,type="C-classification"
                          ,kernel="radial"
                          ,probability = TRUE
                          , cost=c(0.001,0.01,0.1,1,5,10,100,1000),scale=T)
GSE62402.BestSvm<-GSE62402.tObj$best.model
GSE62402.pre_label <- predict(GSE62402.BestSvm, 
                              GSE62402_roc_dat[, GSE62402_genes], 
                              probability = F,
                              decision.values=T)
cal_auc(GSE62402.pre_label,factor(GSE62402_roc_dat$Type))
table(GSE62402.pre_label,factor(GSE62402_roc_dat$Type))

GSE62402.pre_dat <- read.csv('files/svm/GSE62402.svm.csv')
GSE62402.pre_dat <- merge(GSE62402_roc_dat, GSE62402.pre_dat,
                          by = 'Samples')
GSE62402.pre_dat$Type1 <- as.numeric(factor(GSE62402.pre_dat$Type))

GSE62402_roc <- roc(Type1 ~ Values, 
                    data = GSE62402.pre_dat, auc=TRUE, ci=TRUE)
GSE62402_roc_plot <- ggroc(p_roc, size = 1.2, legacy.axes = TRUE)
GSE62402_roc_plot <- GSE62402_roc_plot + 
  theme_bw() + labs(title = 'GSE62402') + 
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate("text", x = .75, y = .25, 
           label = paste("AUC = ", round(p_roc$auc, 3), '(', 
                         round(p_roc$ci[1], 3), '-',  
                         round(p_roc$ci[3], 3), ')'),
           size = 5)
GSE62402_roc_plot
ggsave(plot = GSE62402_roc_plot,
       filename = 'PDFs/GSE62402_roc_plot.pdf',
       width = 5, height = 5)

GSE62402_roc_dat$Type <- as.numeric(factor(GSE62402_roc_dat$Type))
plot_ROC(dat = GSE62402_roc_dat,
         genes = GSE62402_genes,
         outdir = 'PDFs/GSE62402/',
         prefix = 'GSE62402')

# 6.6、GSE7158 的 ROC 曲线 ###############
library(stringr)
GSE7158_cli <- read.delim('raw_dat/osteoporosis/GSE7158_family/SampleInfo_contrib1-GPL570.txt')
View(GSE7158_cli)
GSE7158_cli1 <- GSE7158_cli[, c("Accession", "Title")]
colnames(GSE7158_cli1) <- c('Samples', 'Type')
GSE7158_cli1$Type <- str_split_fixed(GSE7158_cli1$Type, ' ', 3)[, 2]
table(GSE7158_cli1$Type)
GSE7158_cli1 <- GSE7158_cli1[GSE7158_cli1$Type != '', ]
GSE7158_cli1$Type <- ifelse(GSE7158_cli1$Type == 'Low', 'OP', 'Control')
rownames(GSE7158_cli1) <- GSE7158_cli1$Samples

GSE7158_exp <- read.delim('raw_dat/osteoporosis/GSE7158_family/MergeExpro_contrib1-GPL570.txt',
                          row.names = 1)
GSE7158_exp <- GSE7158_exp[which(str_split_fixed(rownames(GSE7158_exp), ' /// ', 2)[, 2] == ''), ]
GSE7158_exp <- log2(GSE7158_exp + 1)
boxplot(GSE7158_exp, las = 2)

GSE7158_genes <- intersect(hub_genes, rownames(GSE7158_exp))

rownames(GSE7158_cli1) <- GSE7158_cli1$Samples
GSE7158_roc_dat <- cbind(GSE7158_cli1, 
                         t(GSE7158_exp[GSE7158_genes, 
                                       GSE7158_cli1$Samples]))

GSE7158.pca <- prcomp(GSE7158_roc_dat[, 3:8],  scale = TRUE)
fviz_pca_ind(GSE7158.pca,  axes = c(1, 2), label="none", 
             addEllipses = F, ellipse.type="norm", ellipse.level=0.95,
             habillage = GSE7158_roc_dat$Type, palette = "jco",
             mean.point=F, legend.title = "Type")

ggbiplot(GSE7158.pca, scale=1, groups = GSE7158_roc_dat$Type,
         ellipse = TRUE,ellipse.prob=0.5, circle = F,var.axes=F) +
  scale_color_discrete(name = 'Type') +
  # scale_fill_manual(values = c("#00468BFF", "#ED0000FF", "#42B540FF")) +
  xlim(-3, 3) + ylim(-3, 3) + theme_bw() +
  theme(legend.direction = 'horizontal', legend.position = 'top')

ggplotViolin_group(GSE7158_roc_dat[, GSE7158_genes],
                   GSE7158_roc_dat$Type,
                   cols = c("#E64B35", "#4DBBD5"),
                   title = 'Cluster',
                   xlab = '',
                   ylab = 'Score',
                   angle = 60)

set.seed(12345)
GSE7158.tObj <- tune.svm(GSE7158_roc_dat[, GSE7158_genes],
                         factor(GSE7158_roc_dat$Type)
                         ,type="C-classification"
                         ,kernel="radial"
                         ,probability = TRUE
                         , cost=c(0.001,0.01,0.1,1,5,10,100,1000),scale=T)
GSE7158.BestSvm<-GSE7158.tObj$best.model
GSE7158.pre_label <- predict(GSE7158.BestSvm, 
                             GSE7158_roc_dat[, GSE7158_genes], 
                             probability = F,
                             decision.values=T)
cal_auc(GSE7158.pre_label,factor(GSE7158_roc_dat$Type))
table(GSE7158.pre_label,factor(GSE7158_roc_dat$Type))

GSE7158.pre_dat <- read.csv('files/svm/GSE7158.svm.csv')
GSE7158.pre_dat <- merge(GSE7158_roc_dat, GSE7158.pre_dat,
                         by = 'Samples')
GSE7158.pre_dat$Type1 <- as.numeric(factor(GSE7158.pre_dat$Type))

GSE7158_roc <- roc(Type1 ~ Values, 
                   data = GSE7158.pre_dat, auc=TRUE, ci=TRUE)
GSE7158_roc_plot <- ggroc(p_roc, size = 1.2, legacy.axes = TRUE)
GSE7158_roc_plot <- GSE7158_roc_plot + 
  theme_bw() + labs(title = 'GSE7158') + 
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate("text", x = .75, y = .25, 
           label = paste("AUC = ", round(p_roc$auc, 3), '(', 
                         round(p_roc$ci[1], 3), '-',  
                         round(p_roc$ci[3], 3), ')'),
           size = 5)
GSE7158_roc_plot
ggsave(plot = GSE7158_roc_plot,
       filename = 'PDFs/GSE7158_roc_plot.pdf',
       width = 5, height = 5)

GSE7158_roc_dat$Type <- as.numeric(factor(GSE7158_roc_dat$Type))
plot_ROC(dat = GSE7158_roc_dat,
         genes = GSE7158_genes,
         outdir = 'PDFs/GSE7158/',
         prefix = 'GSE7158')


GSE_ROC_dat <- cowplot::plot_grid(GSE35959_roc_plot,
                                  GSE7429_roc_plot,
                                  GSE13850_roc_plot,
                                  GSE7158_roc_plot,
                                  GSE62402_roc_plot,
                                  ncol = 3, nrow = 2,
                                  align = 'hv', 
                                  labels = LETTERS[1:5])
GSE_ROC_dat
ggsave(plot = GSE_ROC_dat,
       filename = 'PDFs/GSE_ROC_dat.pdf',
       width = 12, height = 8)
