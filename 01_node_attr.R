### The network of SARS-CoV-2—cancer molecular interactions and pathways
### Pau Erola (2021)
### This code extracts tissue-specific cancer genes and their attributes (oncogene, tumor suppressor gene, driver, citations)


# data source http://bionlp.bcgsc.ca/cancermine/
cancer = read.table("datasets/cancermine/cancermine_collated.tsv", sep="\t", quote = "", header=T)
sel = which(cancer$citation_count>1)
cancer = cancer[sel,]

# We are using grep with different expressions to get 6 tissue-specific cancers
tissues=c("lung", "testic", "colon", "intestinal cancer", "kidney|renal c|renal p|nephr", "stomach|gastric")
fname=  c("lung", "testis", "colon", "intestinal", "kidney", "stomach")


for (t in 1:length(tissues)) {

	g = c("gene", "oncg", "tsg", "driv", "count")

	# select tissue
	sel = grep(tissues[t], cancer$cancer_normalized)
	cancer_subset = cancer[sel,]


	# for each oncogenic gene in that tissue
	gcancer = unique(cancer_subset$gene_normalized)
	m = matrix(0, 4, length(gcancer))
	colnames(m) = gcancer
	rownames(m) = c("oncg", "tsg", "driv", "count")
	m=t(m)
	for (i in 1:nrow(cancer_subset)) {
	
		#"Driver" "Tumor_Suppressor" "Oncogene"
		if ( cancer_subset$role[i] == "Oncogene") {
			m[cancer_subset$gene_normalized[i],1] = 1
			m[cancer_subset$gene_normalized[i],4] = m[cancer_subset$gene_normalized[i],4] + cancer_subset$citation_count[i]
		
		} else if ( cancer_subset$role[i] == "Tumor_Suppressor") {
			m[cancer_subset$gene_normalized[i],2] = 1
			m[cancer_subset$gene_normalized[i],4] = m[cancer_subset$gene_normalized[i],4] + cancer_subset$citation_count[i]

		} else { #"Driver"
			m[cancer_subset$gene_normalized[i],3] = 1
			m[cancer_subset$gene_normalized[i],4] = m[cancer_subset$gene_normalized[i],4] + cancer_subset$citation_count[i]

		}
	}

	# remove oncogenes with less than 3 citations - low confidence
	sel = which(m[,4]>2)
	m = m[sel,]

  	m = cbind( gene=rownames(m), m)
  	write.table(m, paste0("out/node_attr_",fname[t],".txt"), row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')

}





