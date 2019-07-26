#!/usr/bin/env python3
"""
Sets up and submits a SLURM script to run the whole Hi-C scaffolding
process, from 
"""

import argparse
import os
import subprocess

class FileNotReadableError(Exception):
    def __init__(self, path):
        self.message = "File at {} not readable or doesn't exist".format(path)


def check_path_readable(path):
    if not os.access(path, os.R_OK):
        raise FileNotReadableError(path)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-w', '--work-dir', type=os.path.abspath,
                        default=os.path.abspath('.'),
                        help='Directory for intermediate and output files')
    parser.add_argument('-m', '--mem', default='48G',
                        help='Amount of memory to give SLURM task')
    parser.add_argument('-p', '--cpus', default='16',
                        help='Number of CPUs to give SLURM task')
    parser.add_argument('-t', '--time', default='2-00:00:00',
                        help='Amount of time to give SLURM task')
    parser.add_argument('-e', '--enzyme', default='GATC',
    parser.add_argument('contigs', type=os.path.abspath,
                        help='fasta containing contigs to scaffold')
    parser.add_argument('r1_reads', type=os.path.abspath,
                        help='R1 (forward) reads in fastq(.gz) format')
    parser.add_argument('r2_reads', type=os.path.abspath,
                        help='R2 (reverse) reads in fastq(.gz) format')
    args = parser.parse_args()

    # we don't actually want to open the files, because that's for the
    # programs we call to do, but we do want to make sure they exist
    # and are readable/writable before going on
    for f in [args.contigs, args.r1_reads, args.r2_reads]:
        check_path_readable(f)

    return args


def main():
    args = parse_args()

    logdir = os.path.join(args.work_dir, 'logs')
    os.makedirs(logdir, exist_ok = True)

    cmd = ['sbatch', '--mem', args.mem,
           '--cpus-per-task', args.cpus,
           '--time', args.time,
           '-o', os.path.join(logdir, 'hic.out'),
           '-e', os.path.join(logdir, 'hic.err')]
    export_dict = {
            'reference': args.contigs,
            'reads_r1': args.r1_reads,
            'reads_r2': args.r2_reads,
            'program_path': os.path.abspath(os.path.dirname(__file__)),
            'workdir': args.work_dir,
            'enzyme': args.enzyme,
    }
    exports = ','.join(['{}={}'.format(k, v) for k, v in export_dict.items()])
    cmd += ['--export', exports]
    if args.partition:
        cmd += ['-p', args.partition]
    subprocess.run(cmd, check=True, stdout=sys.stdout, stderr=sys.stderr)


if __name__ == '__main__':
    main()


