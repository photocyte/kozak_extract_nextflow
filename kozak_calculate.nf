
kallistoIsoformTMM_ch = file(params.transXpressResults+'kallisto.isoform.TPM.not_cross_norm')
transdecoderGFF_ch =    file(params.transXpressResults+'*.transdecoder.gff3')
trinityFasta_ch1 =      file(params.transXpressResults+'*Trinity.fasta')

//Note, if the channel is empty, it won't execute the downstream processes
//and nextflow gives no error

println kallistoIsoformTMM_ch
println transdecoderGFF_ch
println trinityFasta_ch1

process makeHighlyExpressedList {

input:
  file kallistoIsoformTMM from kallistoIsoformTMM_ch
output:
  file 'highly_expressed_transcripts_regex.txt' into highlyExpressed_ch

script:
"""
cat ${kallistoIsoformTMM} | sort -k2nr,2 | head -n 1500 | cut -f 1 | sed 's/\$/[_.]/g'> highly_expressed_transcripts_regex.txt
""" 
}

process kozakGFF3Extract {
scratch "/lab/solexa_weng/tmp/"
input:
  file transdecoderGFF from transdecoderGFF_ch

output:
  file 'kozaks.gff.gz' into kozaksGff_ch

script:
"""
kozak_gff3_extract.py -g ${transdecoderGFF} --in_memory | gzip > kozaks.gff.gz
"""
}

process kozakFastaExtract {
scratch "/lab/solexa_weng/tmp/"
input:
  file inputGff from kozaksGff_ch
  file inputFasta from trinityFasta_ch1

output:
  file 'all_kozaks.fa.gz' into kozaksFasta_ch

script:
"""
seqkit subseq --gtf ${inputGff} ${inputFasta} | seqkit rmdup -s | gzip >  all_kozaks.fa.gz
"""
}

process kozakFastaExtractHighExp {
scratch "/lab/solexa_weng/tmp/"
input:
  file inputFasta from kozaksFasta_ch
  file highlyExpressedRegexFile from highlyExpressed_ch
output:
  file 'high_only_kozaks.fa' into kozaksFastaHigh_ch

script:
"""
seqkit grep -n -r -f ${highlyExpressedRegexFile} ${inputFasta} > high_only_kozaks.fa
"""
}

process plotKozak {
publishDir params.transXpressResults, mode: 'copy', overwrite: true
tag {params.transXpressResults+'kozak_kpLogo.output.pdf'}
input:
  file inputFasta from kozaksFastaHigh_ch
output:
  file 'kozak_kpLogo.output.pdf' into kozakOutputs
script:
"""
/lab/solexa_weng/testtube/kpLogo-1.1/bin/kpLogo ${inputFasta} -o kozak_kpLogo.output -seq 1 -weight 2 -alphabet dna -max_k 4 -max_shift 1 -startPos 21 -gapped -minCount 0.01 -pseudo 1 -region 1,0 -plot p -pc 0.01 -stack_order 1 -fix 0.75
"""
}

