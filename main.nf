#!/usr/bin/env nextflow

params.enzyme = 'GATC'

reference_file = file(params.reference)
r1_reads = file(params.r1_reads)
r2_reads = file(params.r2_reads)

reads_channel = Channel.from(r1_reads, r2_reads)

process bwa_index {
    publishDir 'bwa_index'

    input:
    file reference from reference_file

    output:
    file "${reference}.*" into genome_index

    """
    bwa index ${reference}
    """
}

process fasta_index {
    input:
    file reference from reference_file

    output:
    file "${reference}.fai" into faidx

    "samtools faidx ${reference}"
}

process bwa_mem {
    cpus 16

    input:
    file ref from reference_file
    file index from genome_index
    file reads from reads_channel

    output:
    file 'out.bam' into mapped

    """
    bwa mem -t 16 ${ref} ${reads} | samtools view -bh - | \
        filter_chimeras.py - > out.bam
    """
}

filtered_pairs = mapped.buffer(size: 2)

process combine {
    publishDir 'alignments'

    input:
    file 'r?.bam' from filtered_pairs

    output:
    file 'combined.bam' into combinedbam

    """
    combine_ends.py r1.bam r2.bam | samtools fixmate -m - - \
        | samtools sort - | samtools markdup -r - combined.bam
    """
}

process bam2bed {
    publishDir 'alignments'

    input:
    file 'combined.bam' from combinedbam

    output:
    file 'combined.bed' into combinedbed

    """
    bedtools bamtobed -i combined.bam | sort -k 4 -T . > combined.bed
    """
}

process salsa {
    publishDir 'salsa_out', mode: 'move'

    input:
    file 'ref.fa' from reference_file
    file 'ref.fa.fai' from faidx
    file 'combined.bed' from combinedbed

    output:
    file 'salsa/*' into salsa_out

    """
    python \$SALSA_DIR/run_pipeline.py -a ref.fa -l ref.fa.fai \
        -b combined.bed -e ${params.enzyme} -o salsa -m yes -i 10 -p yes
    """
}
