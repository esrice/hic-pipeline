# slurm-hic
A pipeline for scaffolding a genome assembly using Hi-C and SALSA2 on SLURM,
generalizable to whatever job scheduler or local system you are using. Loosely
based on the [Arima Hi-C mapping pipeline](https://github.com/ArimaGenomics/mapping_pipeline).

## Options for running
I use [nextflow](https://www.nextflow.io/) to run this pipeline, because it
automatically deals with submitting SLURM jobs and containerization. You can run
it manually too, though, or write your own scripts. There is a section below for
both of these options.

## Requirements
If you're not using docker or singularity, you'll need the following software
installed:
* [bwa](https://github.com/lh3/bwa) (tested with 0.7)
* [samtools](http://www.htslib.org/) (tested with 1.9)
* [bedtools](https://bedtools.readthedocs.io/en/latest/) (tested with 2.26)
* python3
* [SALSA2](https://github.com/marbl/SALSA) (requires python2.7)

No installation is necessary for this pipeline; just download it and put the
resulting folder in your path:
```bash
git clone https://github.com/esrice/slurm-hic.git
export PATH=$PATH:$(PWD)/slurm-hic
```

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
```
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

## Running using nextflow
You can also run this pipeline with a single command using nextflow. If you have
a way to run commands inside containers, such as singularity or docker, you
don't even need to worry about installing requirements! Nextflow is easy to
install:
```
wget -qO- https://get.nextflow.io | bash
```
This command creates the `nextflow` binary in the current directory. You can
put it somewhere in your `$PATH` if you want.

### Configuring nextflow
Nextflow looks for the file `nextflow.config` in the current directory for
instructions on how to run the pipeline, e.g., whether to run it locally or
using a cluster scheduler. You can read the nextflow
[documention](https://www.nextflow.io/docs/latest/config.html) on configuration
for more information, but here's what I use to run the pipeline on SLURM:
```
process {
    executor = 'slurm'
    queue = 'BioCompute'
    time = '2d'
    memory = '48 GB'
}
```
Nextflow also supports SGE, LSF, PBS, AWS, and several other batch systems.

#### Configuring to use locally installed software
Make sure the software specified in the Requirements section is installed and
in your path. If you're in an HPC/cluster environment that uses the `module`
system for software, make sure to tell nextflow to load the correct modules,
e.g., by putting the following line in `nextflow.config`:
```
process.module = 'bwa/0.7:samtools/1.9:bedtools/2.26:salsa/2.2'
```

#### Configuring to use containers
If you can use docker or singularity, you don't need to worry about installing
any requirements. Just add the following lines to `nextflow.config`:
```
process.container = 'esrice/hic-pipeline:latest'
singularity.enabled = true
```
substituting `docker` for `singularity` if using docker. Depending on how you
have singularity set up, you may need to add `singularity.autoMount = true` or
`singularity.runOptions = "--bind /storage"` where `/storage` is the root of
your data filesystem. See nextflow documentation on
[docker](https://www.nextflow.io/docs/latest/docker.html) and
[singularity](https://www.nextflow.io/docs/latest/singularity.html) for more
information.

### Running the pipeline
Just run the following command:
```
nextflow run hic-pipeline.nf --reference contigs.fa \
    --r1_reads HiC-reads_R1.fastq.gz --r2_reads HiC-reads_R2.fastq.gz \
    --enzyme GATC
```
The output will be in `salsa_out/`. If the pipeline gets cut off before
finishing, you can start from where it left off by adding the `-resume` flag.

