# #Xenopus tropicalis genome V10.0
# #Xenopus laevis genome V10.1


###############remove sequence not in chromosome
sed -n '1,14493212p' XENTR_10.0_genome.fasta > XENTR_10.0_genome.clean.fasta
sed -n '1,45693362p' XENLA_10.1_genome.fa > XENLA_10.1_genome.clean.fa

###############create fusion genome
cat XENTR_10.0_genome.clean.fasta XENLA_10.1_genome.clean.fa > XenTr_v10.0_XenLa_v10.1_clean.fasta

###############bwa index
bwa index -p XenTr_v10.0_XenLa_v10.1_clean.bwa_index XenTr_v10.0_XenLa_v10.1_clean.fasta

###############hisat2 index
hisat2-build XenTr_v10.0_XenLa_v10.1_clean.fasta XenTr_v10.0_XenLa_v10.1_clean.hisat2_index