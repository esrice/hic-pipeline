# hic-pipeline
A pipeline for scaffolding a genome assembly using Hi-C and SALSA2,
generalizable to whatever job scheduler or local system you are using. Loosely
based on the [Arima Hi-C mapping pipeline](https://github.com/ArimaGenomics/mapping_pipeline).

*N.B.* (30 September 2022): I have developed a new pipeline for scaffolding with
Hi-C reads. It produces much better results thanks to improved algorithms for
mapping (chromap instead of bwa + postprocessing) and scaffolding (YAHS instead
of SALSA2). I strongly recommend using that one now:

<https://github.com/WarrenLab/hic-scaffolding-nf>

## Options for running
I use [nextflow](https://www.nextflow.io/) to run this pipeline, because it
automatically deals with submitting cluster jobs and containerization. You can run
it manually too, though, or write your own scripts. There is a section below for
both of these options.

## Requirements
### Recommended: use provided docker image
If you can run this pipeline on a system/cluster that supports docker or
singularity (e.g., Lewis), you don't need to install any requirements. The
nextflow pipeline will take care of loading the image and running everything
inside a container.

### Not recommended: manual installation
If you're not using docker or singularity, you'll need the following software
installed:
* [bwa](https://github.com/lh3/bwa) (tested with 0.7)
* [samtools](http://www.htslib.org/) (tested with 1.9)
* [bedtools](https://bedtools.readthedocs.io/en/latest/) (tested with 2.26)
* python3 and [pysam](https://pysam.readthedocs.io/en/latest/api.html)
* [SALSA2](https://github.com/marbl/SALSA) (requires python2.7)

Unfortunately, SALSA2 is written in python2.7, but other parts of the pipeline
are in python3, so you'll need to figure out a workaround for this problem.
Take a look at the Dockerfile to see how I set up the container.

No installation is necessary for this pipeline; just download it and put the
resulting folder in your path:
```bash
git clone https://github.com/esrice/slurm-hic.git
export PATH=$PATH:$(PWD)/slurm-hic
```

## Running using nextflow (recommended)
You can run this pipeline with a single command using nextflow. If you have
a way to run commands inside containers, such as singularity or docker, you
don't even need to worry about installing requirements! Nextflow is easy to
install:
```
wget -qO- https://get.nextflow.io | bash
```
This command creates the `nextflow` binary in the current directory. You can
put it somewhere in your `$PATH` if you want.

### Configuring nextflow
If you're running this on Lewis, the configuration included with the pipeline
should work out of the box, so you can probably skip this section.

Otherwise, copy `nextflow.config` to the folder in which you are going to run
the pipeline and edit accordingly for your system.  You can read the nextflow
[documention](https://www.nextflow.io/docs/latest/config.html) on configuration
for more information. Nextflow supports SLURM, SGE, LSF, PBS, AWS, and other
batch and cloud systems.

### Running the pipeline
Just run the following command:
```
nextflow run esrice/hic-pipeline --reference contigs.fa \
    --r1_reads HiC-reads_R1.fastq.gz --r2_reads HiC-reads_R2.fastq.gz \
    --enzyme GATC
```
This will download the pipeline from github, distribute jobs on your cluster
as specified in `nextflow.config`, and put the output in `salsa_out/`. If the
pipeline gets cut off before finishing, you can start from where it left off
by adding the `-resume` flag.


## Running manually
### Mapping reads
First, index your contigs if you haven't already:
```bash
bwa index contigs.fa
samtools faidx contigs.fa
```

Now, align the R1 and R2 reads to your contigs separately, filtering the output
through the `filter_chimeras.py` script to remove experimental artifacts from
the alignments:
```bash
bwa mem contigs.fa HiC-reads_R1.fastq.gz | samtools view -bh - \
    | filter-chimeras.py - > r1.bam
bwa mem contigs.fa HiC-reads_R2.fastq.gz | samtools view -bh - \
    | filter-chimeras.py - > r2.bam
```

### Combining paired reads
Combine the r1 and r2 files using `combine_ends.py`, and then fix mates, sort,
and remove PCR duplicates:
```bash
combine_ends.py r1.bam r2.bam | samtools fixmate -m - - | samtools sort - \
    | samtools markdup -r - combined.bam
```

### Run SALSA2
SALSA2 takes alignments in bed rather than bam format, sorted by read name.
Here's the command to do that:
```bash
bamToBed -i combined.bam | sort -k 4 > combined.bed
```

Now, just run SALSA. This is the command I use, but see
[the manual](https://github.com/marbl/SALSA) for more options:
```bash
python2.7 run_pipeline.py -a contigs.fa -l contigs.fa.fai -b combined.bed \
    -e GATC -m yes -o salsa -p yes
```

