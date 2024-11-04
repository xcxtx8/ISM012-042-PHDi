library(ggplot2)
library(DESeq2)
library(dplyr)
library(ggfortify)
library(edgeR)
library(scales)
library(fgsea)
library("RColorBrewer")
library(pheatmap)


################################################
# read in data
df.count = read.csv('exp_mtx/count.csv',
                    as.is = T,
                    row.names = 1)


# generate metadata
meta.df = data.frame(name = colnames(df.count))

# decomposite metadata
meta.df$days = as.character(mapply(function(x) strsplit(x,'_')[[1]][1], meta.df$name))
meta.df$treatment = as.character(mapply(function(x) strsplit(x,'_')[[1]][2], meta.df$name))


# convert to int count matrix
df.int = data.frame(lapply(df.count, as.integer))
rownames(df.int) = rownames(df.count)
colnames(df.int) = colnames(df.count)

# build DESeq object
meta.df$group = paste0(meta.df$days, meta.df$treatment)
dds = DESeqDataSetFromMatrix(countData = data.matrix(df.int),
                             colData = meta.df,
                             design = ~ group)


####################################
# DEG pairwise comparsion
system('mkdir DEG_Deseq2_with_sample_Exp')

f.Deseq_compare = function(dds, group2, group1){
  dds$group = relevel(dds$group, group2, group1)
  dds = DESeq(dds)
  out = results(dds, contrast = c('group', group1, group2))
  countMtx = counts(dds, normalized = T)
  group1Sample = colnames(dds)[dds$group == group1]
  group2Sample = colnames(dds)[dds$group == group2]
  countMtx = countMtx[,c(group1Sample, group2Sample)]
  if(all(rownames(countMtx) == rownames(out))){
    out = cbind(out, countMtx)
  }
  
  
  outName = paste0(group1,'_vs_',group2,'.csv')
  write.csv(out,paste0('DEG_Deseq2_with_sample_Exp/',outName))
  return(out)
}

DEGList[['day2ISM012_vs_day2Vehicle']] = f.Deseq_compare(dds, 'day2Vehicle', 'day2ISM012')
DEGList[['day2Vehicle_vs_day2sham']] = f.Deseq_compare(dds, 'day2sham', 'day2Vehicle')

#########################
# GSEA
system('mkdir GSEA_pathway_analysis')
pathwayDB_kegg = gmtPathways('pathway_db/kegg_mouse_symbol.gmt')
gseaList_kegg = list()

for(comparison in c('day2ISM012_vs_day2Vehicle', 'day2Vehicle_vs_day2sham')){
  
  geneList = DEGList[[comparison]]$log2FoldChange
  names(geneList) = rownames(DEGList[[comparison]])
  geneList[is.na(geneList)] = 0
  geneList = geneList[-which(geneList == 0)]
  geneList = sort(geneList, decreasing = TRUE)
  
  fgRes = fgsea(pathways = pathwayDB_kegg, 
                stats = geneList,
                minSize=15,
                maxSize=600) %>% 
    as.data.frame()
  fgRes$leadingEdge = mapply(function(x) paste0(fgRes$leadingEdge[x]), 1:nrow(fgRes))
  write.csv(fgRes,paste0('GSEA_pathway_analysis/',comparison,'.csv'))
  
  gseaList_kegg[[comparison]] = fgRes
  
}


######################### get top 20 pathway
day2_topPathway_list = list()
for(comparison in names(gseaList_kegg)){
  
  p.df = gseaList_kegg[[comparison]]
  p.df = p.df %>% filter(padj < 0.05) %>% arrange(pval)
  p.df = p.df[1:20,]
  p.df = p.df[order(p.df$NES, decreasing = T),]
  p.df$pathway = factor(p.df$pathway, levels = rev(p.df$pathway))

  if(grepl('day2',comparison)){day2_topPathway_list[[comparison]] = p.df$pathway}
}


commonPath = day2_topPathway_list$day2Vehicle_vs_day2sham
p.df1 = gseaList_kegg$day2Vehicle_vs_day2sham
p.df1$group = 'Veh_vs_sham'
p.df2 = gseaList_kegg$day2ISM012_vs_day2Vehicle
p.df2$group = 'ISM042_vs_Veh'
p.df = rbind(p.df1, p.df2)
p.df$col = ifelse(p.df$NES > 0, 'pos' ,'neg')
ord = p.df %>% filter(group == 'Veh_vs_sham') %>% arrange(-padj) %>% pull(pathway)
p.df$pathway = factor(p.df$pathway, levels = ord)
p.df$group = factor(p.df$group, levels = c('Veh_vs_sham','ISM042_vs_Veh'))
p.df = p.df %>% filter(pathway %in% commonPath)
p.df$col_con = mapply(function(x) ifelse(p.df$NES[x] > 0, -log10(p.df$padj[x]), log10(p.df$padj[x]) * 2), 1:nrow(p.df)) %>% as.numeric()
p = ggplot(p.df, aes(x = group, y = pathway)) +
  geom_point(aes(col = col_con, size = abs(NES)), shape = 16) +
  labs(size = 'ABS (NES)') +
  theme_classic() +
  labs(col = 'Enrichment') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = 'bold', size = 14),
        axis.text.y = element_text(face = 'bold', size = 14, color = 'black'),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  ggtitle('Day 2') +
  scale_color_gradient2(low = 'blue',mid = 'white',high = 'red', name = '-log10 (p.adj)',
                        breaks = c(10,5,0,-5,-10), labels = c(10,5,0,2.5,5)) +
  ggnewscale::new_scale_color() +
  geom_point(data = data.frame(t = c('Pos Enrich','Neg Enrich'), group = c(1,2), pathway = c(1,2)), aes(col = t), size = 5) +
  scale_color_manual(values = c('Pos Enrich' = 'red', 'Neg Enrich' = 'blue'), name = '') +
  theme(
    legend.title = element_text(size=14),
    plot.title = element_text(size = 18)
  )
p
