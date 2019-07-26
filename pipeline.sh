#!/bin/bash

# environment variables to set: $reference, $reads_r1, $reads_$2,
#   $program_path, $workdir, $enzyme

set -e

module load bwa/bwa-0.7.17 samtools/samtools-1.9 bedtools/bedtools-2.26.0

# generate index if it does not exist
if [ ! -f "${reference}.amb" ]; then
    bwa index $reference
fi

# align reads and postprocess alignments
let "half_of_cpus = SLURM_CPUS_PER_TASK / 2"
bwa mem -t $half_of_cpus $reference $reads_r1 | \
    samtools view -@ $half_of_cpus -bh - | \
    $program_path/filter-chimeras.py - > $workdir/r1.f.bam &
bwa mem -t $half_of_cpus $reference $reads_r2 | \
    samtools view -@ $half_of_cpus -bh - | \
    $program_path/filter-chimeras.py - > $workdir/r2.f.bam &
wait

$program_path/combine_ends.py $workdir/r1.f.bam $workdir/r2.f.bam \
    | samtools fixmate -m - - \
    | samtools sort -@ $SLURM_CPUS_PER_TASK - \
    | samtools markdup -r - $workdir/combined.bam

# make links files
bamToBed -i $workdir/combined.bam | sort -k 4 > $workdir/combined.bed

# run SALSA
python2.7 $salsa/run_pipeline.py -a $reference -l ${reference}.fai \
    -b $workdir/combined.bed -e $enzyme -o $workdir -m yes
