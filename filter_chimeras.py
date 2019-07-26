#!/usr/bin/env python3
"""
Given a bam file containing alignments of either R1 or R2 of a Hi-C
library to a reference genome, select the single aligned segment for
each read that is most useful for scaffolding, if any.

Loosely based on
https://github.com/ArimaGenomics/mapping_pipeline/blob/master/filter_five_end.pl
"""

import argparse
import sys
from enum import Enum

import pysam

# this is the int code used by pysam to indicate an 'M' in a CIGAR
# string, but does not seem to be exposed as a constant anywhere
BAM_CMATCH = 0

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    # cannot use type=lambda because template is needed, which comes
    # from the infile
    parser.add_argument('-o', '--outfile', default='-', nargs='?',
                        help='where to put filtered bam')
    parser.add_argument('in_bam', type=lambda f: pysam.AlignmentFile(f, 'rb'),
                        nargs='?',
                        help='bam to remove chimeric alignments from')
    return parser.parse_args()


def is_useful_segment(segment):
    """
    Returns True if this segment is useful for scaffolding -- i.e.,
    it is mapped and not clipped or gapped at the 5' end.

    Args:
        segment (pysam.AlignedSegment): the segment to examine

    Returns: True if the segment is useful for scaffolding,
        false otherwise.
    """
    if segment.is_unmapped:
        return False
    if ((not segment.is_reverse
            and segment.cigartuples[0][0] == BAM_CMATCH)
            or (segment.is_reverse
            and segment.cigartuples[-1][0] == BAM_CMATCH)):
        return True
    else:
        return False


def choose_segment_to_write(segments):
    """
    Given a list of alignment segments from the same read, determines
    which should be written based on usefulness for scaffolding.

    To choose a segment, this function looks for segments with unclipped
    5' ends. If there is exactly one segment like this, it is returned;
    otherwise, the first segment is returned with the 'unmapped' flag
    set to True.

    Args:
        segments (list): a list of pysam.AlignedSegment items coming
            from the same read

    Returns:
        best_segment (pysam.AlignedSegment): the segment that should
            be written for this read.
    """
    # make a list of segments in this read that have unclipped 5' ends
    useful_segments = list(filter(is_useful_segment, segments))
    if len(useful_segments) == 1:
        return useful_segments[0]
    else:
        segments[0].is_unmapped = True
        return segments[0]


def main():
    args = parse_args()
    outfile = pysam.AlignmentFile(args.outfile, 'wb', template=args.in_bam)

    segments_in_read = []
    for segment in args.in_bam:
        # if we are looking at a new read now, we need to output the
        # good part of the previous read before moving on
        if (len(segments_in_read) != 0 and
                segments_in_read[-1].query_name != segment.query_name):
            outfile.write(choose_segment_to_write(segments_in_read))
            segments_in_read = []

        segments_in_read.append(segment)

    outfile.write(choose_segment_to_write(segments_in_read))


if __name__ == '__main__':
    main()


