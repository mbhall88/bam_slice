import pysam
import os
import sys
import bisect
import math
import resource
import gzip
import click
import io
import pandas as pd
from typing import List, Tuple

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
FILE_LIMIT = resource.getrlimit(resource.RLIMIT_NOFILE)[0] - 5
OUT = io.StringIO()


def merge_overlap_intervals(intervals: List[List[int]]) -> List[List[int]]:
    """Checks consecutive intervals and if they overlap it merges them into a
    single interval.

    Args:
        intervals: A list of intervals where each interval is a List with two
        elements corresponding to the start and end of the interval
        respectively.

    Returns:
        A new intervals list where any intervals that overlapped have been
        merged into a single interval.

    Example:
        >>> intervals = [[1, 4], [3, 7], [10, 14]]
        >>> merge_overlap_intervals(intervals)
        [[1, 7], [10, 14]]

    """
    merged_intervals = []
    cached_interval = None

    for interval in intervals:
        if cached_interval is None:
            cached_interval = interval
            continue

        if outside_interval(cached_interval, interval):
            merged_intervals.append(cached_interval)
            cached_interval = interval
        else:
            cached_interval = extend_interval(cached_interval, interval)

    if cached_interval is not None:
        merged_intervals.append(cached_interval)

    return merged_intervals


def outside_interval(first_interval: List[int],
                     second_interval: List[int]) -> bool:
    """Determines whether two intervals overlap.

    Args:
        first_interval: The interval with the lower start index.
        second_interval: The interval with the higher start index.

    Returns:
        Whether the start of the second interval is less than the end of the
        first interval. i.e do they overlap?

    Notes:
        If the end index of the first interval is equal to the start of the
        second interval than they are deemed to NOT be overlapping.

    Example:
        >>> first_interval = [0, 4]
        >>> second_interval = [3, 7]
        >>> outside_interval(first_interval, second_interval)
        False
    """
    return second_interval[0] >= first_interval[1]


def extend_interval(interval_to_extend: List[int],
                    interval: List[int]) -> List[int]:
    """Extends an interval to encompass another.

    Args:
        interval_to_extend: The interval to extend.
        interval: The interval to extend by.

    Returns:
        A new interval with the same start as interval_to_extend and the same
        end as interval.

    """
    interval_to_extend[1] = interval[1]
    return interval_to_extend


def handle_nones(aligned_pairs: List[Tuple[int, int]]) -> \
        Tuple[List[int], List[int]]:
    """Replaces all instances of None with the preceeding element.

    Args:
        aligned_pairs: A list of tuples of two integers.

    Returns:
        Two lists (in a tuple) that are the "unzipping" of the list of tuples
        passed in. All instances of None are replaced with the previous element
        in that position of the tuple (or -1 if at the start).

    Example:
        >>> aligned_pairs = [(1, 1), (None, 2), (2, None)]
        >>> handle_nones(aligned_pairs)
        ([1, 1, 2], [1, 2, 2])

    """
    handled_ref_positions = []
    handled_read_positions = []
    previous_ref_pos = -1
    previous_read_pos = -1
    for read_pos, ref_pos in aligned_pairs:
        if ref_pos is None:
            handled_ref_positions.append(previous_ref_pos)
        else:
            handled_ref_positions.append(ref_pos)
            previous_ref_pos = ref_pos

        if read_pos is None:
            handled_read_positions.append(previous_read_pos)
        else:
            handled_read_positions.append(read_pos)
            previous_read_pos = read_pos
    return handled_read_positions, handled_ref_positions


def ref_interval_indexes(interval: List[int], ref: List[int]) -> \
        Tuple[int, int]:
    """Gets interval's indexes within ref.

    Args:
        interval: A list of two integers.
        ref: A list of integers.

    Returns:
        A tuple of two integers where the first element is the left-most index
        witin ref of the first element in interval and the second element is
        the right-most index within ref of the second element in interval.

    Example:
        >>> interval = [2, 5]
        >>> ref = [1, 1, 2, 2, 2, 3, 4, 5, 5]
        >>> ref_interval_indexes(interval, ref)
        (2, 8)
        >>> ref = [2, 3, 4]
        >>> ref_interval_indexes(interval, ref)
        (0, 3)
        >>> ref = [3, 4, 4, 4, 4]
        >>> ref_interval_indexes(interval, ref)
        (0, 5)
        >>> ref = [10, 11, 12]
        >>> ref_interval_indexes(interval, ref)
        (0, 0)
    """
    start, end = interval
    return bisect.bisect_left(ref, start), bisect.bisect_right(ref, end)


def get_subsequence_indexes(interval: List[int],
                            clean_aligned_pairs: Tuple[List[int], List[int]]) \
        -> Tuple[int, int]:
    """Given an interval (assumed to be from the second list in
    clean_aligned_pairs) will return the corresponding interval in the first
    list in clean_aligned_pairs.

    Args:
        interval: An interval of two integers assumed to be from the second
        list in clean_aligned_pairs (reference positions).
        clean_aligned_pairs: A tuple of two lists that each contain list of two
        integers. There are assumed to be no Nones in either list.

    Returns:
        A list of two integers (an interval) that is the interval corresponding
        to the passed interval. Returns None if interval is outside alignment.

    Notes:
        It is effectively assumed that clean_aligned_pairs is a variable
        returned from the handle_nones method.

    Example:
        >>> interval = [2, 5]
        >>> clean_aligned_pairs = ([30, 31, 32, 33, 34, 35],
                                   [1, 2, 3, 4, 5, 6])
        >>> get_subsequence_indexes(interval, clean_aligned_pairs)
        [31, 34]

    """
    read_positions, ref_positions = clean_aligned_pairs

    # binary search for intervals start and end indexes within ref positions
    start_index, end_index = ref_interval_indexes(interval, ref_positions)

    # get the corresponding read indexes for given ref indexes
    read_indexes = read_positions[start_index: end_index]

    # if no match was found return nothing
    if not read_indexes:
        return None

    read_start_index = read_indexes[0]
    read_end_index = read_indexes[-1]

    # if there is a deletion in start of read and it crosses the start of the
    # interval, round it up to one
    if read_start_index < 0:
        read_start_index = 0

    return read_start_index, read_end_index


def dump_interval_fastq(interval: List[int],
                        clean_aligned_pairs: Tuple[List[int], List[int]],
                        read: pysam.AlignedSegment,
                        fastq_file: gzip.GzipFile) -> None:
    """For a given interval, and aligned pairs for a read to a reference, pull
    out the read subsequence within that interval and the sub-quality sequence
    and write an entry to fastq. If interval does not occur within read, then
    nothing is written.

    Args:
        interval: A list of two integers.
        clean_aligned_pairs: Two lists of integers in a tuple.
        read: An alignment read within a samfile from pysam.
        fastq_file: File handle to the file for the interval.

    """
    sequence = read.query_sequence
    quality = read.query_qualities
    subsequence_indexes = get_subsequence_indexes(interval, clean_aligned_pairs)
    if subsequence_indexes is not None:  # interval is in read
        start, end = subsequence_indexes
        read_subsequence = sequence[start: end + 1]
        quality_subsequence = quality[start: end + 1]
        append_fastq(read_subsequence, quality_subsequence,
                     read.query_name, fastq_file)


def append_fastq(subsequence: str, quality: List[int],
                 header: str, fastq_file: gzip.GzipFile) -> None:
    """Writes an entry to a fastq file.

    Args:
        subsequence: The DNA sequence to write.
        quality: A list of phred quality scores as integers.
        header: the read name from bam/sam file.
        fastq_file: A file handle for the fastq file to write.

    """
    quality_string = "".join([chr(x + 33) for x in quality])
    fastq_file.write("@{0}\n{1}\n+\n{2}\n".format(header, subsequence,
                                                  quality_string))


def get_filename_prefix(fname: str) -> str:
    """Pulls out the filename without the extension(s). Also removed .gz
    extension before parsing.

    Args:
        fname: A filename or path.

    Returns:
        A string of the filename without it's extensions (or prefixed path).

    Example:
        >>> fname = 'path/to/my.fastq.gz'
        >>> get_filename_prefix(fname)
        'my'

    """
    if fname.endswith('.gz'):
        fname_prefix = os.path.splitext(os.path.basename(fname[:-3]))[0]
    else:
        fname_prefix = os.path.splitext(os.path.basename(fname))[0]

    return fname_prefix


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('sam',
                type=click.Path(exists=True, dir_okay=False, resolve_path=True))
@click.argument('positions_file')
@click.option('--output', '-o',
              default='.',
              type=click.Path(dir_okay=True, resolve_path=True,
                              writable=True),
              help="Directory to save the fastq files to. If not specified,"
                   " will use the current directory. [.]")
@click.option('--column', '-c',
              default=0,
              type=int,
              help="Column within positions file containing the positions. [0]")
@click.option('--delimiter', '-d',
              default='\t',
              type=str,
              help="Delimiter for the positions file. [\\t]")
@click.option('--padding', '-p',
              default=100,
              type=int,
              help="Take this number of bases either side of position in the "
                   "slice.")
def main(sam, positions_file, output, column, delimiter, padding):
    """A script to slice regions out of a [BS]AM file.

    Positional Arguments (required):\n
        SAM: A BAM or SAM file to slice positions from. Must be indexed.\n
        POSITIONS_FILE: A file containing positions to cut out from SAM. The
        delimiter of the file and the column containing the positions can be
        specified by the --delimiter and --column options respectively. Will
        read from STDIN if passed the - character instead of a path. If using
        STDIN it is assumed you are piping only positions and no other data.
        The file should also not contain a header.

    """
    fname_prefix = os.path.join(output, get_filename_prefix(sam))

    # read in snps file and get reference positions for all variants
    if positions_file == '-':
        positions = sorted(sys.stdin.readlines())
    else:
        data = pd.read_csv(positions_file, sep=delimiter, header=None)
        positions = sorted(data[column].drop_duplicates())

    # add padding to either side of positions to create intervals and then
    # merge overlapping intervals.
    intervals = []
    for pos in positions:
        # if taking padding off position goes below 0, make it 0 instead
        start = pos - padding if pos - padding >= 0 else 0
        intervals.append([start, pos + padding])
    merged_intervals = merge_overlap_intervals(intervals)

    # batch intervals into lots of FILE_LIMIT as you cannot have more than that
    # number of files open at once
    batches = math.ceil(len(merged_intervals) / FILE_LIMIT)

    # loop through each batch and write fastq files
    for i in range(int(batches)):
        batch_fastq_files = {}  # holds all open files for batch
        batch_merged_intervals = merged_intervals[i * FILE_LIMIT:
                                                  (i + 1) * FILE_LIMIT]

        with pysam.AlignmentFile(sam, 'rb') as samfile:
            # loop through each read in the BAM/SAM file
            for read in samfile:
                # skip non-primary alignments
                if (read.is_unmapped or
                        read.is_secondary or
                        read.is_supplementary):
                    continue
                # get the index correspondences between read and reference
                aligned_pairs = read.get_aligned_pairs()
                # strip out None values and set them to previous non-None value
                clean_aligned_pairs = handle_nones(aligned_pairs)

                # for each interval within this read, write entry to fastq
                for interval in batch_merged_intervals:
                    filename = "{0}_{1}-{2}.fastq.gz".format(fname_prefix,
                                                             interval[0],
                                                             interval[1])
                    if filename not in batch_fastq_files:
                        file_ = gzip.GzipFile(filename=filename, mode='w',
                                              fileobj=OUT)
                        batch_fastq_files[filename] = file_
                    fastq_file = batch_fastq_files[filename]

                    # check if interval is in read and write fastq entry for it
                    dump_interval_fastq(interval, clean_aligned_pairs,
                                        read, fastq_file)

            # close all files in batch
            for open_file in batch_fastq_files.values():
                open_file.close()


if __name__ == "__main__":
    main()
