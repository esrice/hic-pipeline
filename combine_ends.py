#!/usr/bin/env python3
"""
Given an R1 and an R2 bam file that have been filtered to remove the non-5' ends
of chimeric alignments (with remove_chimeras.py), combine them into a single
bam file.

Loosely based on
https://github.com/ArimaGenomics/mapping_pipeline/blob/master/two_read_bam_combiner.pl
"""

import argparse
import sys

import pysam

class MismatchedReadsError(Exception):
    def __init__(self, line_number, r1_name, r2_name):
        self.message = """Mismatched R1 and R2 read names on the {}th read of
            R1 ({}) and R2 ({})""".format(line_number, r1_name, r2_name)
        super().__init__(self.message)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    # cannot use type=lambda because template is needed, which comes
    # from the infile
    parser.add_argument('-o', '--outfile', default='-', nargs='?',
                        help='where to put filtered bam')
    parser.add_argument('-q', '--mapq', type=int, default=20,
                        help='minimum MAPQ for both ends to keep a pair [20]')
    parser.add_argument('r1_bam', type=lambda f: pysam.AlignmentFile(f, 'rb'),
                        help='bam file containing R1 alignments')
    parser.add_argument('r2_bam', type=lambda f: pysam.AlignmentFile(f, 'rb'),
                        help='bam file containing R2 alignments')
    return parser.parse_args()


def pair_reads(r1, r2):
    """
    Given bam entries for two ends of the same read (R1/R2) that have
    been aligned as if they were single-end reads, set the fields in the
    entries so that they are now properly paired.

    Args:
        r1, r2 (pysam.AlignedSegment): two ends of the same read, both
            aligned and having the same read name

    Returns:
        r1, r2 (pysam.AlignedSegment): the same two reads that were
            provided as arguments, but with the FLAG, RNEXT, PNEXT, and
            TLEN fields modified to make them a proper pair
    """
    # if the two ends map to the same reference, we need to set RNEXT
    # to '=' for both reads and also calculate the TLEN
    if r1.reference_name == r2.reference_name:
        r1.next_reference_name = '='
        r2.next_reference_name = '='
        if r1.reference_start < r2.reference_start:
            tlen = (r2.reference_start
                    + max(r2.get_reference_positions())
                    - r1.reference_start)
            r1.template_length = tlen
            r2.template_length = -1 * tlen
        else:
            tlen = (r1.reference_start
                    + max(r1.get_reference_positions())
                    - r2.reference_start)
            r1.template_length = -1 * tlen
            r2.template_length = tlen
    else: # ends map to different references, so just set RNEXT
        r1.next_reference_name = r2.reference_name
        r2.next_reference_name = r1.reference_name

    # set PNEXT
    r1.next_reference_start = r2.reference_start
    r2.next_reference_start = r1.reference_start

    # set some bits in the FLAG
    r1.is_paired, r2.is_paired = True, True
    r1.is_proper_pair, r2.is_proper_pair = True, True
    r1.mate_is_unmapped, r2.mate_is_unmapped = False, False
    r1.mate_is_reverse = r2.is_reverse
    r2.mate_is_reverse = r1.is_reverse
    r1.is_read1, r2.is_read1 = True, False
    r1.is_read2, r2.is_read2 = False, True
    r1.is_secondary, r2.is_secondary = False, False
    r1.is_supplementary, r2.is_supplementary = False, False

    return r1, r2


def main():
    args = parse_args()
    outfile = pysam.AlignmentFile(args.outfile, 'wb', template=args.r1_bam)

    for i, (r1, r2) in enumerate(zip(args.r1_bam, args.r2_bam)):
        if r1.query_name != r2.query_name:
            raise MismatchedReadsError(i, r1.query_name, r2.query_name)

        if (not r1.is_unmapped and not r2.is_unmapped
                and r1.mapping_quality >= args.mapq
                and r2.mapping_quality >= args.mapq):

            r1, r2 = pair_reads(r1, r2)

            outfile.write(r1)
            outfile.write(r2)


if __name__ == '__main__':
    main()


