#!/usr/bin/env python

import argparse
import re
import sys


def get_dinucleotide(fasta_filename, dinucleotide):
    """Parses FASTA for dinucleotide loci to print as BED format.
    """
    # Make sure the file exists, otherwise throw an error
    try:
        f = open(fasta_filename, 'r')
        f.close()
    except FileNotFoundError as e:
        sys.stderr.write("FileNotFoundError: %s\n" % str(e))
    current_contig = None
    loci_tracker = 0
    with open(fasta_filename, 'r') as handle:
        for line in handle:
            # Checks for first line / first contig
            # In case there's unnecessary header before the contigs
            # (Should not happen, but just in case.)
            if current_contig == "None" and line[0] == ">":
                # Strinp the new line, skip the '>' and split 
                # the line by white space
                # Contig name is first words before white space
                current_contig = line.strip('\n')[1:].split(' ')[0]
                # Reset tracking your position in chromosome
                loci_tracker = 0
                continue # continue to next loop
            elif line[0] == ">": # Checks if new contig
                current_contig = line.strip('\n')[1:].split(' ')[0]
                loci_tracker = 0
                continue # continue to next loop
            seq = line.strip('\n') # remove the new line character
            # Find location of all dinucleotides in the single line
            match_indices = [m.start() for m in re.finditer(dinucleotide, seq)]
            # Loop through all the found dinucleotides
            for match in match_indices:
                # Create the columns of the BED file
                out_cols = [
                    current_contig,
                    str(loci_tracker + match),
                    str(loci_tracker + match + len(dinucleotide))
                ]
                # Print out to columns tab-delimited to the terminal
                sys.stdout.write('\t'.join(out_cols) + '\n')
            # Before moving to the next loop, update loci_tracker
            loci_tracker += len(seq)



def main():
    # Argument parsing
    parser = argparse.ArgumentParser(
        "Finds dinucleotide loci in FASTA and outputs loci to BED"
    )
    parser.add_argument(
        "--dinucleotide", metavar="diNT", default="TA", 
        help="Dinucleotide to search"
    )
    parser.add_argument(
        "fasta", metavar="FASTA",
        help="FASTA file to search"
    )
    args = parser.parse_args()
    # Check that the dinucleotide is only two nts long
    if len(args.dinucleotide) != 2:
        sys.stderr.write("Input dinucleotide is of the wrong size.\n")
        sys.exit()
    # Start the actual finding
    get_dinucleotide(args.fasta, args.dinucleotide)

if __name__ == "__main__":
    main()
