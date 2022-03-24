### The network of SARS-CoV-2—cancer molecular interactions and pathways
### Pau Erola (2021)
### This code retrives PPI between SARS-CoV-2—targeted genes, cancer genes and genetic risk factors


if (!requireNamespace("epigraphdb", quietly = TRUE)) install.packages("epigraphdb")
library(epigraphdb)

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("biomaRt", quietly = TRUE)) BiocManager::install("biomaRt")
library(biomaRt)


# load biomart
mart_grch37 <- useEnsembl(biomart="ensembl",GRCh=37)
mart_grch37 <- useDataset("hsapiens_gene_ensembl", mart_grch37)


###
### SARS-COV-2
###
sarscov2 = read.table("datasets/gordon/SARS-CoV-2_HighConfidence.txt", sep="\t", quote="", header=T)
prey = unique(sarscov2$PreyGene)

# convert gene names to uniprot using biomart
attributes <- c("hgnc_symbol","uniprotswissprot")
filters <- c("hgnc_symbol")
values <- list(hgnc_symbol=prey)
q_covid <- getBM(attributes=attributes, filters=filters, values=values, mart=mart_grch37)
q_covid <- q_covid[which(q_covid$uniprotswissprot!=""),]



###
### Genetic risk factors
### source https://doi.org/10.1038/s41586-020-03065-y
###
risk = matrix( c( "rs73064425", 3, 45901089, "LZTFL1",
		"rs9380142", 6, 29798794, "HLA-G",
		"rs143334143", 6, 31121426, "CCHCR1",
		"rs3131294", 6, 32180146, "NOTCH4",
		"rs10735079", 12, 113380008, "OAS1–OAS3",
		"rs2109069", 19, 4719443, "DPP9",
		"rs74956615", 19, 10427721, "TYK2",
		"rs2236757", 21, 34624917, "IFNAR2",
		"rs71325088", 3, 45862952, "LZTFL1",
		"rs6489867", 12, 113363550, "OAS1–OAS3",
		"rs11085727", 19, 10466123, "TYK2",
		"rs13050728", 21, 34615210, "IFNAR2" ),
		ncol=4, byrow=TRUE )
colnames(risk) = c("SNP","chr","pos","gene_symbol")
risk = as.data.frame(risk)

# convert gene names to uniprot using biomart
attributes <- c("hgnc_symbol","uniprotswissprot")
filters <- c("hgnc_symbol")
values <- list(hgnc_symbol=risk$gene_symbol)
q_risk <- getBM(attributes=attributes, filters=filters, values=values, mart=mart_grch37)
q_risk <- q_risk[which(q_risk$uniprotswissprot!=""),]

q_covid = rbind(q_covid, q_risk)



###
### Cancer genes
### source http://bionlp.bcgsc.ca/cancermine/
###
cancer = read.table("datasets/cancermine/cancermine_collated.tsv", sep="\t", quote = "", header=T)
# remove low confidence
sel = which(cancer$citation_count>1) 
cancer = cancer[sel,]

# convert gene names to uniprot using biomart
attributes <- c("hgnc_symbol","uniprotswissprot")
filters <- c("hgnc_symbol")
values <- list(hgnc_symbol=cancer$gene_normalized)
q_cancer <- getBM(attributes=attributes, filters=filters, values=values, mart=mart_grch37)
q_cancer <- q_cancer[which(q_cancer$uniprotswissprot!=""),]




###
### Find PPI using EpiGraphDB
###
ppi <- cypher(paste0("MATCH (p1:Protein)-[]-(p2:Protein)",
			" WHERE p1.uniprot_id IN ['", paste(q_covid$uniprotswissprot, collapse="','"),
			"'] AND p2.uniprot_id IN ['", paste(q_cancer$uniprotswissprot, collapse="','"),
			"'] RETURN p1.uniprot_id, p2.uniprot_id") )

ppi <- data.frame(ppi)
ppi = unique(as.matrix(ppi, ncol=2))

q_allgenes = rbind(q_covid, q_cancer)
ppidf = data.frame( gname = rep(NA,nrow(ppi)), uniprot = rep(NA,nrow(ppi)) )

# convert uniprot to gene names
for (i in 1:nrow(ppi)) {

	sel = which(q_allgenes$uniprotswissprot == ppi[i,1])[1]
	ppidf[i,1] = q_allgenes[sel,"hgnc_symbol"]
	
	sel = which(q_allgenes$uniprotswissprot == ppi[i,2])[1]
	ppidf[i,2] = q_allgenes[sel,"hgnc_symbol"]
	
}


write.table(ppidf, "out/ppi_covidcancer.txt", row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\t')


