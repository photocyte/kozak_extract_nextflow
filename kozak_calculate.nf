nextflow.enable.dsl=2

process makeHighlyExpressedList {

input:
  path kallistoIsoformTMM from kallistoIsoformTMM_ch
output:
  path 'highly_expressed_transcripts_regex.txt' into highlyExpressed_ch

script:
"""
cat ${kallistoIsoformTMM} | sort -k2nr,2 | head -n 1500 | cut -f 1 | sed 's/\$/[_.]/g'> highly_expressed_transcripts_regex.txt
""" 
}

process kozakGFF3Extract {
publishDir 'results', mode: 'link', overwrite: true
conda "biopython gffutils"
input:
  path transdecoderGFF

output:
  path 'kozaks.gff.gz'

script:
"""
set -o pipefail
kozak_gff3_extract.py -g ${transdecoderGFF} --in_memory --i 20 | gzip > kozaks.gff.gz
"""
}

process gff_standardize {
conda "genometools-genometools"
input:
  path transdecoderGFF

output:
  path "gt.$transdecoderGFF"

script:
"""
gt gff3 -tidy -sort -retainids $transdecoderGFF | gzip > gt.$transdecoderGFF
"""
}


process kozakFastaExtract {
publishDir 'results', mode: 'link', overwrite: true
conda "seqkit"
input:
  tuple path(inputGff),path(inputFasta)

output:
  path 'all_kozaks.fa'

script:
"""
set -o pipefail
seqkit subseq --gtf ${inputGff} ${inputFasta} | seqkit rmdup -s >  all_kozaks.fa
"""
}

process kozakFastaExtract_2 {
publishDir 'results', mode: 'link', overwrite: true
conda "genometools-genometools"
input:
  tuple path(inputGff),path(inputFasta)

output:
  path 'all_kozaks.fa'

script:
"""
set -o pipefail
gt extractfeat -seqid -matchdescstart -retainids -coords -type kozakseq -seqfile ${inputFasta} ${inputGff} > all_kozaks.fa
"""
}



process kozakFastaExtractHighExp {
input:
  path inputFasta from kozaksFasta_ch
  path highlyExpressedRegexFile from highlyExpressed_ch
output:
  path 'high_only_kozaks.fa' into kozaksFastaHigh_ch

script:
"""
set -o pipefail
seqkit grep -n -r -f ${highlyExpressedRegexFile} ${inputFasta} > high_only_kozaks.fa
"""
}

process plot_kpLogo {
publishDir 'results', mode: 'link', overwrite: true
tag {'kozak_kpLogo.output.pdf'}
input:
  path inputFasta
output:
  path 'kozak_kpLogo.output.pdf'
script:
"""
kpLogo ${inputFasta} -o kozak_kpLogo.output -seq 1 -weight 2 -alphabet dna -max_k 4 -max_shift 1 -startPos 21 -gapped -minCount 0.01 -pseudo 1 -region 1,0 -plot p -pc 0.01 -stack_order 1 -fix 0.75
"""
}

process plot_weblogo {
publishDir 'results', mode: 'link', overwrite: true
conda "weblogo"
input:
  path inputFasta
output:
  path "kozak_weblogo.output.pdf"
script:
"""
weblogo -F pdf < ${inputFasta} > kozak_weblogo.output.pdf
"""
}

workflow fromTransxpress {
kallistoIsoformTMM_ch = Channel.fromPath(params.transXpressResults+'kallisto.isoform.TPM.not_cross_norm')
transdecoderGFF_ch =    Channel.fromPath(params.transXpressResults+'*.transdecoder.gff3')
trinityFasta_ch1 =      Channel.fromPath(params.transXpressResults+'*Trinity.fasta')

println kallistoIsoformTMM_ch
println transdecoderGFF_ch
println trinityFasta_ch1


//Note, if the channel is empty, it won't execute the downstream processes
//and nextflow gives no error
}

workflow {
theGff = Channel.fromPath(params.gff)
theFasta = Channel.fromPath(params.fasta) 
kozakGFF3Extract(theGff) | gff_standardize
gff_standardize.out.combine(theFasta) | kozakFastaExtract_2 | (plot_kpLogo & plot_weblogo)

}
