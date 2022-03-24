### The network of SARS-CoV-2—cancer molecular interactions and pathways
### Pau Erola (2021)
### We use binary matrix format to show the shared cancer genes and pathways potentially perturbed in each cancer type

if (!require("stringr")) install.packages("stringr")
if (!require("ggcorrplot")) install.packages("ggcorrplot")
if (!require("readxl")) install.packages("readxl")
library("ggcorrplot")
library("readxl")   
library("stringr")


###
### Binary matrix for pathways
###

cancer = read.table("datasets/cancermine/cancermine_collated.tsv", sep="\t", quote = "", header=T)
# note we are removing low confidence entries with less than 2 citations
sel = which(cancer$citation_count>2)
cancer = cancer[sel,]

tissues=c("lung", "testic", "colon", "intestinal cancer", "kidney|renal c|renal p|nephr", "stomach|gastric")
fname=  c("lung", "testis", "colon", "intestinal", "kidney", "stomach")
fname2=  c("lung", "testis", "colon", "small\n intestine", "kidney", "stomach")
x = list()
x2 = list()
for (t in 1:length(tissues)) {

	sel = grep(tissues[t], cancer$cancer_normalized)
	cancer_subset = cancer[sel,]
	gcancer = unique(cancer_subset$gene_normalized)

  
	tmpf = read.table(paste0("out/net_",fname[t],".txt"), header=TRUE, sep="\t")
	x[[fname[t]]]  = intersect( gcancer,unique(c(as.character(tmpf$dest), as.character(tmpf$orig))) )
	x2[[fname[t]]] = unique(c(as.character(tmpf$dest) ) )
}

# create a comma-separated string of genes for each tissue 
# we use comma as a delimiter to make sure str_count doesn't count for partial names		
glist = list()
glist[["lung"]] = paste0(",",gsub(" ", "", toString(exceln[["lung"]]$Nodes,sep=",") ),",") 
glist[["testis"]] = paste0(",",gsub(" ", "", toString(exceln[["testis"]]$Nodes,sep=",") ),",") 
glist[["colon"]] = paste0(",",gsub(" ", "", toString(exceln[["colon"]]$Nodes,sep=",") ),",") 
glist[["intestinal"]] = paste0(",",gsub(" ", "", toString(exceln[["intestinal"]]$Nodes,sep=",") ),",") 
glist[["kidney"]] = paste0(",",gsub(" ", "", toString(exceln[["kidney"]]$Nodes,sep=",") ),",") 
glist[["stomach"]] = paste0(",",gsub(" ", "", toString(exceln[["stomach"]]$Nodes,sep=",") ),",") 


genes = sort(unique(as.character(unlist(x))))
m = matrix(0, 6, length(genes))
m2 = matrix(0, 6, length(genes))
colnames(m) = genes
rownames(m) = fname
colnames(m2) = genes
rownames(m2) = fname
for (i in 1:6) {
	for (p in genes) {
		if (p %in% x[[i]]) {
			m[i,p]=1
			m2[i,p] = str_count( glist[[i]], paste0(",",p,",") )
		}
	}
}
  
m=t(m)
m2=t(m2)
sel1=which(rowSums(m)>1)
sel2=which(rowSums(m2)>1)
#sel=intersect(sel1,sel2)
sel=sel1
m2[m2==0] <- NA

colnames(m) = c("lung", "testis", "colon", "small intestine", "kidney", "stomach")
colnames(m2) = c("lung", "testis", "colon", "small intestine", "kidney", "stomach")

write.table(m,"out/binmat_genes.txt", sep="\t", quote=F, col.names=T, row.name=T)


tiff("out/binmat_gene.tif", heigh=500, width=500+nrow(m)*13)
ggcorrplot(data.frame((m[sel,])),
		tl.cex = 26,
		tl.col = "black",
		tl.srt = 45,
		digits = 2,
		show.legend = FALSE,
		colors = c("yellow", "white", "black"),
		outline.color = "white") +
geom_vline(xintercept=1:nrow(m[sel,])-0.5, colour="white", size=4) +
geom_hline(yintercept=1:ncol(m[sel,])-0.5, colour="white", size=4) +
aes(label=as.character(m2[sel,])) +
geom_text(size=7, colour="yellow")
dev.off()




###
### Binary matrix for pathways
###

excelf = list()
excelf[["lung"]] = as.character(read_excel("cytoscape/gsea_cancer.xlsx", sheet="lung")$`GeneSet`)
excelf[["testis"]] = as.character(read_excel("cytoscape/gsea_cancer.xlsx", sheet="testis")$`GeneSet`)
excelf[["colon"]] = as.character(read_excel("cytoscape/gsea_cancer.xlsx", sheet="colon")$`GeneSet`)
excelf[["intestinal"]] = as.character(read_excel("cytoscape/gsea_cancer.xlsx", sheet="intestinal")$`GeneSet`)
excelf[["kidney"]] = as.character(read_excel("cytoscape/gsea_cancer.xlsx", sheet="kidney")$`GeneSet`)
excelf[["stomach"]] = as.character(read_excel("cytoscape/gsea_cancer.xlsx", sheet="stomach")$`GeneSet`)

exceln = list()
exceln[["lung"]] = data.frame(read_excel("cytoscape/gsea_cancer.xlsx", sheet="lung"))
exceln[["testis"]] = data.frame(read_excel("cytoscape/gsea_cancer.xlsx", sheet="testis"))
exceln[["colon"]] = data.frame(read_excel("cytoscape/gsea_cancer.xlsx", sheet="colon"))
exceln[["intestinal"]] = data.frame(read_excel("cytoscape/gsea_cancer.xlsx", sheet="intestinal"))
exceln[["kidney"]] = data.frame(read_excel("cytoscape/gsea_cancer.xlsx", sheet="kidney"))
exceln[["stomach"]] = data.frame(read_excel("cytoscape/gsea_cancer.xlsx", sheet="stomach"))


pathways = unique(as.character(unlist(excelf)))
m = matrix(0, 6, length(pathways))
m2 = matrix(0, 6, length(pathways))
colnames(m) = pathways
rownames(m) = fname
colnames(m2) = pathways
rownames(m2) = fname
for (i in 1:6) {
	
	for (p in pathways) {
		if (p %in% excelf[[i]]) {
			m[i,p]=1
			m2[i,p]= max(exceln[[i]][which(exceln[[i]]$GeneSet==p),"ProteinsFromModule"])
		}
	}
}

# manually shortening the names for better display
colnames(m) = c("cytokine-mediated sig. pathway", "reg. of DNA-binding tr. fact. activity", "cellular response to cadmium ion", "neg. reg. of autophagy", "telomere maint. via semi-conservative repl.", "DNA replication, synthesis of RNA primer", "protein ubiquitination", "pos. reg. of kinase activity", "pos. reg. of proteasomal ubiquitin-dependent protein catabolic pr.", "neg. reg. of canonical Wnt sig. pathway", "mRNA splicing, via spliceosome", "pos. reg. of NF-kappaB tr. fact. activity", "hippo sig.", "reg. of cell diff. inv. in embryonic placenta dev.", "primitive hemopoiesis", "endocardium dev.", "cell diff. inv. in embryonic placenta dev.", "neural tube formation", "neg. reg. of organ growth", "hepatocyte apoptotic pr.", "pos. reg. of extrinsic apoptotic sig. pathway via death…", "stress-activated protein kinase sig. cascade", "U4 snRNA 3'-end pr.ing", "nuclear-transcribed mRNA catabolic pr., exonucleolytic, 3'-5'", "exonucleolytic catabolism of deadenylated mRNA", "nuclear polyadenylation-dependent tRNA catabolic pr.", "polyadenylation-dependent snoRNA 3'-end pr.ing", "exonucleolytic trimming to generate mature 3'-end of 5.8S rRNA…", "reg. of mRNA stability", "nuclear mRNA surveillance", "rRNA pr.ing", "rRNA catabolic pr.", "reg. of signal transduction by p53 class mediator", "G1/S transition of mitotic cell cycle", "DNA replication initiation", "pos. reg. of fat cell diff.", "Wnt sig. pathway", "beta-catenin-TCF complex assembly", "reg. of cell motility", "transmembrane receptor protein tyrosine kinase…", "peptidyl-tyrosine phosphorylation", "Golgi organization", "endoplasmic reticulum to Golgi vesicle-mediated transport", "activation of protein kinase A activity", "cellular response to glucagon stimulus", "renal water homeostasis", "post-translational protein modification", "reg. of transcription from RNA polymerase II promoter in hypoxia", "protein polyubiquitination", "reg. of ERK1 and ERK2 cascade", "pos. reg. of snRNA transcription by RNA polymerase II", "neg. reg. of type I interferon production")

m=t(m)
m2=t(m2)
m2[m2==0] <- NA

colnames(m) = c("lung", "testis", "colon", "small intestine", "kidney", "stomach")

write.table(m,"out/binmat_pathways.txt", sep="\t", quote=F, col.names=T, row.name=T)


tiff("out/binmat_pathway.tif", heigh=800, width=700+nrow(m)*25)
ggcorrplot(data.frame((m[rowSums(m)>1,])), # only shared pathways are shown
		tl.cex = 26,
		tl.col = "black",
		tl.srt = 45,
		digits = 2,
		show.legend = FALSE,
		colors = c("yellow", "white", "black"),
		outline.color = "white") +
geom_vline(xintercept=1:nrow(m[rowSums(m)>1,])-0.5, colour="white", size=4) +
geom_hline(yintercept=1:ncol(m[rowSums(m)>1,])-0.5, colour="white", size=4) +
aes(label=as.character(m2[rowSums(m)>1,])) +
geom_text(size=7, colour="yellow")
dev.off()







