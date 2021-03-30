GFF=GWHAZHV00000000.gff
FASTA=GWHAZHV00000000.genome.fasta.gz

##GFF=../PPYR_OGS/PPYR_OGS1.1.gff3
##FASTA=../PPYR_OGS/Genome_release.fa

##GFF=GCF_008802855.1_Ppyr1.3_genomic.gff.gz
##FASTA=GCF_008802855.1_Ppyr1.3_genomic.fna.gz

nextflow run kozak_calculate.nf --gff ${GFF} --fasta ${FASTA} -resume
