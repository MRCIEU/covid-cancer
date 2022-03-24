### The network of SARS-CoV-2—cancer molecular interactions and pathways
### Pau Erola (2021)
### Reconstructs tissue-specific networks of SARS-CoV-2—targeted genes, cancer genes and genetic risk factors



###
### PPI
###
human_ppi = read.table("out/ppi_covidcancer.txt", sep="\t", header=F)
colnames(human_ppi) = c("gene1","gene2")


###
### SARS-COV-2
###
sarscov2 = read.table("datasets/gordon/SARS-CoV-2_HighConfidence.txt", sep="\t", quote="", header=T)


###
### Genetic risk factors
### source https://doi.org/10.1038/s41586-020-03065-y
###
# ! rs74956615 mapped to RAVER1
riskfactor = matrix( c( "rs73064425", 3, 45901089, "LZTFL1",
		"rs9380142", 6, 29798794, "HLA-G",
		"rs143334143", 6, 31121426, "CCHCR1",
		"rs3131294", 6, 32180146, "NOTCH4",
		"rs10735079", 12, 113380008, "OAS1–OAS3",
		"rs2109069", 19, 4719443, "DPP9",
		"rs74956615", 19, 10427721, "RAVER1",
		"rs2236757", 21, 34624917, "IFNAR2",
		"rs71325088", 3, 45862952, "LZTFL1",
		"rs6489867", 12, 113363550, "OAS1–OAS3",
		"rs11085727", 19, 10466123, "TYK2",
		"rs13050728", 21, 34615210, "IFNAR2" ),
		ncol=4, byrow=TRUE )
colnames(riskfactor) = c("SNP","chr","pos","gene_symbol")
riskfactor = as.data.frame(riskfactor)


###
### Cancer genes
### source http://bionlp.bcgsc.ca/cancermine/
###
cancer = read.table("datasets/cancermine/cancermine_collated.tsv", sep="\t", quote = "", header=T)
# remove low confidence
sel = which(cancer$citation_count>1) 
cancer = cancer[sel,]



###
### Tissue definitions (one tissue might have different cancer types)
###
tissues=c("lung", "testic", "colon", "intestinal cancer", "kidney|renal c|renal p|nephr", "stomach|gastric")
fname=  c("lung", "testis", "colon", "intestinal", "kidney", "stomach")
#> unique(cancer$cancer_normalized[grep("lung", cancer$cancer_normalized)])
#[1] "lung non-small cell carcinoma" "lung cancer"                  
#[3] "lung adenocarcinoma"           "lung small cell carcinoma"    
#[5] "lung squamous cell carcinoma"  "lung carcinoma"               
#[7] "lung giant cell carcinoma"     "lung large cell carcinoma"    
#> unique(cancer$cancer_normalized[grep("testic", cancer$cancer_normalized)])
#[1] "testicular germ cell cancer" "testicular cancer"          
#> unique(cancer$cancer_normalized[grep("stomach|gastric", cancer$cancer_normalized)])
#[1] "stomach cancer"                    "stomach carcinoma"                
#[3] "hereditary diffuse gastric cancer" "diffuse gastric cancer"           
#[5] "gastric adenocarcinoma"           
#> unique(cancer$cancer_normalized[grep("colon", cancer$cancer_normalized)])
#[1] "colon cancer"         "colon carcinoma"      "colon adenocarcinoma"
#> unique(cancer$cancer_normalized[grep("intestinal cancer", cancer$cancer_normalized)])
#[1] "intestinal cancer"
#> unique(cancer$cancer_normalized[grep("kidney|renal c|renal p|nephr", cancer$cancer_normalized)])
#[1] "nephroblastoma"                     "clear cell renal cell carcinoma"   
#[3] "renal cell carcinoma"               "kidney cancer"                     
#[5] "familial renal papillary carcinoma" "childhood kidney cancer"           
#[7] "papillary renal cell carcinoma"     "renal carcinoma"                   
#[9] "adrenal cortex cancer"             



###
### Reconstruction of tissue-specific networks
###
for (t in 1:length(tissues)) {
	network = c("orig","orig_type","dest","dest_type", "int")

	sel = grep(tissues[t], cancer$cancer_normalized)
	cancer_subset = cancer[sel,]

	gvirus = unique(sarscov2$PreyGene)
	gcancer = unique(cancer_subset$gene_normalized)
	grisk = unique(riskfactor$gene_symbol)

	# ppi between target gene and cancer gene
	for (i in 1:nrow(human_ppi)) {

		if ((human_ppi$gene1[i] %in% gvirus) & (human_ppi$gene2[i] %in% gcancer)) { 
			network = rbind(network, c(human_ppi$gene1[i], "gvirus", human_ppi$gene2[i], "gcancer", "cov2") )
			
		} else if ((human_ppi$gene2[i] %in% gvirus) & (human_ppi$gene1[i] %in% gcancer)) {
			network = rbind(network, c(human_ppi$gene2[i], "gvirus", human_ppi$gene1[i], "gcancer", "cov2") )
		}	
	}

	network = network[-1,]
	colnames(network) = c("orig","orig_type","dest","dest_type", "int")



	# ppi between cancer gene and genetic risk factor
	for (i in 1:nrow(human_ppi)) {

		if ((human_ppi$gene1[i] %in% grisk) & (human_ppi$gene2[i] %in% gcancer)) { 
			network = rbind(network, c(human_ppi$gene1[i], "grisk", human_ppi$gene2[i], "gcancer", "risk") )
			
		} else if ((human_ppi$gene2[i] %in% grisk) & (human_ppi$gene1[i] %in% gcancer)) {
			network = rbind(network, c(human_ppi$gene2[i], "grisk", human_ppi$gene1[i], "gcancer", "risk") )
		}	
	}


	# this tables can be imported in Cytoscape
	write.table(network,paste0("out/net_",fname[t],".txt"), sep="\t", quote=F, col.names=T, row.name=F)

}




