##merge Xtr Xla gene expression

gene_trans_table_file <- "gene_id_trans_table_three_species.txt"
Xtr_gene_exp_table_file <- "Hybridize.featureCount.Xtr.matrix"
Xla_gene_exp_table_file <- "Hybridize.featureCount.Xla.matrix"
gene_trans_table <- read.table(file = gene_trans_table_file, sep = "\t", quote = "", header = T, colClasses = "character")
Xtr_gene_exp_table <- read.table(file = Xtr_gene_exp_table_file, header = F, sep = "\t", colClasses = "character")
Xla_gene_exp_table <- read.table(file = Xla_gene_exp_table_file, header = F, sep = "\t", colClasses = "character")
names(Xtr_gene_exp_table) <- Xtr_gene_exp_table[1, ]
Xtr_gene_exp_table <- Xtr_gene_exp_table[-1, ]
names(Xla_gene_exp_table) <- Xla_gene_exp_table[1, ]
Xla_gene_exp_table <- Xla_gene_exp_table[-1, ]
################### Xla gene merge
Xla_gene_exp_table$gene <- gene_trans_table[match(Xla_gene_exp_table$gene, gene_trans_table$XB_GENE_NAME), 4]
Xla_gene_exp_table <- na.omit(Xla_gene_exp_table)
sample_names <- names(Xla_gene_exp_table)
Xla_gene_exp_table[, -1] <- as.data.frame(lapply(Xla_gene_exp_table[, -1], as.numeric))
Xla_gene_exp_table <- aggregate(as.matrix(Xla_gene_exp_table[, -1]), by = list("gene" = Xla_gene_exp_table$gene), FUN = sum)

################# Xtr gene convert
Xtr_gene_exp_table$gene <- gene_trans_table[match(Xtr_gene_exp_table$gene, gene_trans_table$XB_GENE_NAME), 4]
Xtr_gene_exp_table <- na.omit(Xtr_gene_exp_table)
Xtr_gene_exp_table[, -1] <- as.data.frame(lapply(Xtr_gene_exp_table[, -1], as.numeric))

################## merge Xtr and Xla gene
merged_gene_exp_table <- rbind(Xtr_gene_exp_table, Xla_gene_exp_table[na.omit(match(Xtr_gene_exp_table$gene, Xla_gene_exp_table$gene)), ])
merged_gene_exp_table[, -1] <- as.data.frame(lapply(merged_gene_exp_table[, -1], as.numeric))
merged_gene_exp_table <- aggregate(merged_gene_exp_table[, -1], by = list("gene" = merged_gene_exp_table$gene), FUN = sum)


write.table(merged_gene_exp_table, "Hybridize.count.merge.matrix", quote = F, sep = "\t", row.names = F)