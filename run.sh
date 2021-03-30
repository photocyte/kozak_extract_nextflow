##GFF=~/scratch/Projects/20200925_UTEX2797_assembly_v0.2/UTEX2797_genomic.v0.2.gff
##FASTA=~/scratch/Projects/20200925_UTEX2797_assembly_v0.2/UTEX2797_genome.v0.2.fasta

##GFF=../PPYR_OGS/PPYR_OGS1.1.gff3
##FASTA=../PPYR_OGS/Genome_release.fa

##GFF=GCF_008802855.1_Ppyr1.3_genomic.gff.gz
##FASTA=GCF_008802855.1_Ppyr1.3_genomic.fna.gz

GFF=GWHAZHV00000000.gff
FASTA=GWHAZHV00000000.genome.fasta.gz

nextflow run kozak_calculate.nf --gff ${GFF} --fasta ${FASTA} -resume
