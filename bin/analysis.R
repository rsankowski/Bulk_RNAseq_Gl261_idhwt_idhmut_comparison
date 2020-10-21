library(tidyverse)
library(biomaRt)
require(clusterProfiler)
require(org.Mm.eg.db)
keytypes(org.Mm.eg.db)
library(pheatmap)
library(viridis)
library(RColorBrewer)

source("bin/functions.R")

counts <- read_delim("data/Galaxy28-[Column_Join_on_data_26_and_data_24].tabular", 
                     delim="\t",
                     col_types = c("cdddd"))[,c(1:2, 4:5)]
colnames(counts) <- c("geneid", "IDHwt", "IDHmut", "geneLength")

#convert gene ids
gene <- bitr(counts$geneid, fromType = "ENTREZID",
                   toType = c("SYMBOL", "ENTREZID"),
                   OrgDb = org.Mm.eg.db)

counts <- counts %>%
  left_join(gene, by=c("geneid" = "ENTREZID")) %>%
  na.omit() %>%
  mutate(IDHwt_rpkm = IDHwt / ( geneLength/1000 * sum(IDHwt)/1000000 ),
         IDHmut_rpkm = IDHmut / ( geneLength/1000 * sum(IDHmut)/1000000)) %>%
  as.data.frame()

rownames(counts) <- counts$SYMBOL
counts2 <- counts[,c("IDHwt_rpkm", "IDHmut_rpkm")]
colnames(counts2) <- gsub("_.*", "", colnames(counts2))
counts2["Cd44",]

#define gene signatures
genes <- c( "Cd44", "Tnfrsf1a", "Mgmt", "Col4a2", "Timp1", "Timp4", "Myc", "Fam20c", "Serpine1", "Angptl4",
            "Abcc3", "Dll3", "Srrm2",  "Ncam1", "Scg3", "Sox2","Sox8", "Nes", "Olig2", "Tagln","Lif", "Fosl2" 
            )


heatm <- pheatmap(counts2[genes,],
         color = colorRampPalette(rev(brewer.pal(n =3, name = "RdYlBu")))(10),
         cluster_rows = F, 
         cluster_cols = F, 
         scale = "row", 
         cellwidth = 19, 
         cellheight = 19,
         lwd=0.25,
         border_color = 'black',
         angle_col=45)

pdf("plots/proneur_mesench_heatmap_idhwt_vs_mut.pdf")
heatm
dev.off()

counts2[genes,] %>%
  rownames_to_column(var="Gene") %>%
  pivot_longer(IDHwt_rpkm :IDHmut_rpkm, names_to = "Genotype", values_to="RPKM") %>%
  mutate(Genotype = gsub("_.*", "", .$Genotype),
         Gene = factor(.$Gene, levels = genes)) %>%
  ggplot(aes(x=Genotype, y=Gene, fill=log2(RPKM))) +
  geom_tile() 
