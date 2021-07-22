#!/usr/bin/env python

import argparse
import re
import sys


def find_indices(prev, current, nuc):
    """Finds matching indices from previous seq plus one less than nuc
    seq from next line. (Finds overlap between lines.)
    """
    if current != None:
        buf = current[0: len(nuc)-1]
    else:
        buf = ""
    return [
        m.start() 
        for m in re.finditer(
            nuc, 
            prev + buf
        )
    ]


def print_bed(indices, contig, loci, stepper):
    for match in indices:
        out_cols = [
            contig,
            str(loci + match),
            str(loci + match + stepper)
        ]
        sys.stdout.write('\t'.join(out_cols) + '\n')


def get_dinucleotide(fasta_filename, nuc_seq):
    """Parses FASTA for dinucleotide loci to print as BED format.

    Note: we delay printing until the next line so that we can
    add a buffer to the line to account for spanning nucleotide sequence.
    """
    nuc_comp_map = {
        "A" : "T",
        "T" : "A",
        "C" : "G",
        "G" : "C"
    }
    comp_nuc_seq = ''.join([nuc_comp_map[c] for c in nuc_seq])
    # Make sure the file exists, otherwise throw an error
    try:
        f = open(fasta_filename, 'r')
        f.close()
    except FileNotFoundError as e:
        sys.stderr.write("FileNotFoundError: %s\n" % str(e))
    current_contig = None
    loci_tracker = 0
    first_line = True
    with open(fasta_filename, 'r') as handle:
        for line in handle:
            # Checks for first line / first contig
            # In case there's unnecessary header before the contigs
            # (Should not happen, but just in case.)
            if current_contig == None and line[0] == ">":
                # Strip the new line, skip the '>' and split 
                # the line by white space
                # Contig name is first words before white space
                current_contig = line.strip('\n')[1:].split(' ')[0]
                # Reset tracking your position in chromosome
                loci_tracker = 0
                continue # continue to next loop
            elif line[0] == ">": # Checks if new contig
                # print last sequence line from previous contig
                indices = sorted(
                    find_indices(
                        prev_seq, None, nuc_seq
                    ) + find_indices(
                        prev_seq, None, comp_nuc_seq
                    )
                )
                print_bed(
                    indices,
                    current_contig,
                    loci_tracker,
                    len(nuc_seq)
                )
                # then reset everything for the new contig
                current_contig = line.strip('\n')[1:].split(' ')[0]
                loci_tracker = 0
                first_line = True
                continue # continue to next loop
            elif first_line == True:
                # if first line of sequence, store and move to next line
                prev_seq = line.strip('\n')
                first_line = False
                continue
            seq = line.strip('\n') # remove the new line character
            # Find location of all dinucleotides in the previous line
            indices = sorted(
                find_indices(
                    prev_seq, seq, nuc_seq
                ) + find_indices(
                    prev_seq, seq, comp_nuc_seq
                )
            )
            print_bed(
                indices,
                current_contig,
                loci_tracker,
                len(nuc_seq)
            )
            # Before moving to the next loop, update loci_tracker
            loci_tracker += len(prev_seq)
            prev_seq = seq
    # Print final line of sequence b/c we delay the seq addition
    indices = sorted(
        find_indices(
            prev_seq, None, nuc_seq
        ) + find_indices(
            prev_seq, None, comp_nuc_seq
        )
    )
    print_bed(
        indices,
        current_contig,
        loci_tracker,
        len(nuc_seq)
    )



def main():
    # Argument parsing
    parser = argparse.ArgumentParser(
        "Finds dinucleotide loci in FASTA and outputs loci to BED"
    )
    parser.add_argument(
        "--nucleotides", metavar="NT", default="TA", 
        help="Nucleotide sequence to search"
    )
    parser.add_argument(
        "fasta", metavar="FASTA",
        help="FASTA file to search"
    )
    args = parser.parse_args()
    # Check that the dinucleotide is only two nts long
    #if len(args.dinucleotide) != 2:
    #    sys.stderr.write("Input dinucleotide is of the wrong size.\n")
    #    sys.exit()
    # Start the actual finding
    get_dinucleotide(args.fasta, args.nucleotides)

if __name__ == "__main__":
    main()
