### The network of SARS-CoV-2—cancer molecular interactions and pathways
### Pau Erola (2021)
### Plots Venn diagrams comparing the number of cancer genes and pathways potentially perturbed in each cancer type

# install.packages("VennDiagram")
# install.packages("venn")
# install.packages("readxl")
library("VennDiagram")
library("venn")
library("readxl")   


palette = c("#377EB8","#FFFF33","#E41A1C","#4DAF4A","#984EA3","#FF7F00")
#c("blue","yellow","red","green","purple","orange")

fname2=  c("lung", "testis", "colon", "intestinal", "kidney", "stomach")
fname2=  c("lung", "testis", "colon", "small\n intestine", "kidney", "stomach")


###
### Venn for genes
###
x = list()
x1 = list()
x2 = list()
for (i in 1:length(fname)) {

	tmpf = read.table(paste0("out/net_",fname[i],".txt"), header=TRUE, sep="\t")
	x[[fname2[i]]] = unique(c(as.character(tmpf$dest), as.character(tmpf$orig)))
	x1[[fname2[i]]] = unique(as.character(tmpf$dest))
	x2[[fname2[i]]] = unique(as.character(tmpf$orig))
}


tiff(paste0("out/venn_genes.tif"), 1000, 1000)
venn::venn( x1,
		cexil = 2,
		cexsn = 2,
		ilcs = 2.3,
		sncs = 2.5,
		zcolor=palette[1:6]
)
dev.off()

capture.output(get.venn.partitions(x1), file = "out/venn_genes.txt")




###
### Venn for pathways
###
excelf = list()
excelf[["lung"]] = as.character(read_excel("cytoscape/gsea_cancer.xlsx", sheet="lung")$`GeneSet`)
excelf[["testis"]] = as.character(read_excel("cytoscape/gsea_cancer.xlsx", sheet="testis")$`GeneSet`)
excelf[["colon"]] = as.character(read_excel("cytoscape/gsea_cancer.xlsx", sheet="colon")$`GeneSet`)
excelf[["small\n intestine"]] = as.character(read_excel("cytoscape/gsea_cancer.xlsx", sheet="intestinal")$`GeneSet`)
excelf[["kidney"]] = as.character(read_excel("cytoscape/gsea_cancer.xlsx", sheet="kidney")$`GeneSet`)
excelf[["stomach"]] = as.character(read_excel("cytoscape/gsea_cancer.xlsx", sheet="stomach")$`GeneSet`)


tiff(paste0("out/venn_pathway.tif"), 1000, 1000)
venn::venn( excelf,
		cexil = 2,
		cexsn = 2,
		ilcs = 2.3,
		sncs = 2.5,
		zcolor=palette[1:6]
)
dev.off()

capture.output(get.venn.partitions(excelf), file = "out/venn_pathway.txt")




