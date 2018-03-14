# bam_slice
Cut slices out of a BAM file given a list of positions.  
 
Briefly, this script takes a position or list of positions and generates an interval of 100bp either side of this position 
(padding can be specified with `--padding` argument). If any intervals overlap, they are merged.  
Next, the list of intervals to slice out are batched into lots dependent on the limit of open files on the machine the 
script is running on (As a fastq file needs to be opened for each interval).  
It then iterates over every read in a given BAM file. For each read in the BAM file, if the read is not a primary alignment 
it is skipped. Otherwise, the subsequence of sequence and quality is sliced out of the read for each interval is falls 
within and written as a record in a fastq file for the interval.  

The help menu for the script looks like 

```
Usage: bam_slice.py [OPTIONS] SAM POSITIONS_FILE

  A script to slice regions out of a [BS]AM file.

  Positional Arguments (required):

      SAM: A BAM or SAM file to slice positions from. Must be indexed.

      POSITIONS_FILE: A file containing positions to cut out from SAM. The
      delimiter of the file and the column containing the positions can be
      specified by the --delimiter and --column options respectively. Will
      read from STDIN if passed the - character instead of a path. If using
      STDIN it is assumed you are piping only positions and no other data.
      The file should also not contain a header.

Options:
  -o, --output PATH      Directory to save the fastq files to. If not
                         specified, will use the current directory. [.]
  -c, --column INTEGER   Column within positions file containing the
                         positions. [0]
  -d, --delimiter TEXT   Delimiter for the positions file. [\t]
  -p, --padding INTEGER  Take this number of bases either side of position in
                         the slice.
  -h, --help             Show this message and exit.
```

An example of how to run the script is

```sh
POSITIONS_FILE=positions.snps
SAM=reads.bam
OUT=slice_fastq_files/ 
python bam_slice/bam_slice.py -o "$OUT" "$BAM" "$POSITIONS_FILE"
```

After running this, the output directory contains a fastq file fo each interval/position. The interval is built into the 
file name like this `reads_262-462.fastq`  

As the help file mentions you could also pipe your positions into the script if you don't have them stored in a nice neat 
form. A simplistic version of this would be a file `positions.txt` that looks like

```
some_text	12
some_text	222
some_text	353
some_text	45435
```

Two ways you could run this file.  
Pipe into the script

```sh
POSITIONS_FILE=positions.txt
SAM=reads.bam
OUT=slice_fastq_files/ 
cut -f2 "$POSITIONS_FILE" | python bam_slice/bam_slice.py -o "$OUT" "$BAM" -
```

Or via the CLI options

```sh
POSITIONS_FILE=positions.txt
SAM=reads.bam
OUT=slice_fastq_files/ 
python bam_slice/bam_slice.py --delimiter '\t' --column 1 -o "$OUT" "$BAM" "$POSITIONS_FILE"
```
