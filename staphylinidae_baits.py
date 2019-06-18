# coding: utf8
"""
Ortholog Based Bait Design Script for creating Staphylinidae ortholog based baits suitable submission to myBaits
Compares t_coffee AA alignment scores with nucleotide tranalignments to find conserved blocks
Author Jackson Eyres jackson.eyres@canada.ca
License: MIT
Copywright: Government of Canada
"""

import glob
import os
from Bio import AlignIO, SeqIO
import time
import argparse


def main():
    """
    Main Function to run Staphylinidae Bait Designer
    :return:
    """
    parser = argparse.ArgumentParser(description='Processes T_Coffee AA alignments to generate a ortholog bait set')
    parser.add_argument('-o', type=str, required=True,
                        help='Output Directory')
    parser.add_argument('-i', type=str, required=True,
                        help='T_Coffee Directory containing aa based .score_ascii files')
    parser.add_argument('-n', type=str, required=True,
                        help='Directory containing tranalign nucleotide alignments')
    parser.add_argument('-p', type=str, required=True,
                        help='Priorities File for Staphylinidae')

    args = parser.parse_args()
    print("Starting Staphylinidae Ortholog Bait Design".format(args.o))

    print(args.o, args.i, args.n)

    dict_of_max_sums = longest_exon_length(args.i)
    sum_file = write_sums(args.o, dict_of_max_sums)

    blocks_dir = extract_conserved_blocks(sum_file, args.n, args.o)

    window_ranges = [500, 400, 300]
    for window in window_ranges:
        filtered_blocks_dir = filter_blocks(blocks_dir, args.p, args.o, window)
        processed_blocks_dir = filtered_blocks_dir

        # Original was going to stagger tile the baits, but bait manufacturer inherently does this
        # tiled_blocks_dir = tile_blocks(filtered_blocks_dir, args.o, window)
        # processed_blocks_dir =  tiled_blocks_dir

        merge_baits(processed_blocks_dir, args.o, "Staphylinidae", window)


def extract_conserved_blocks(sum_file, alignment_directory, results_directory):
    """
    Takes an AA T_coffee alignment score_ascii file, the corresponding nt fasta tranalign file, and the sum file to
    Extract out a conserved block
    :param sum_file:
    :param alignment_directory:
    :param results_directory:
    :return: Output Directory of conserved blocks
    """
    output_directory = os.path.join(results_directory, "conserved_blocks")
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    with open(sum_file) as f:
        lines = f.readlines()
        lines.pop(0)

        for line in lines:
            list_of_seqs = []
            split = line.rstrip().split(",")

            name = split[0].replace(".aa.summarized.score_ascii", "_nt_renamed.fasta")
            window_range = int(split[2])*3
            index = int(split[3])*3

            file_path = os.path.join(alignment_directory, name)
            if os.path.isfile(file_path):
                with open(file_path) as g:
                    alignments = AlignIO.read(g, "fasta")
                    for alignment in alignments:
                        list_of_seqs.append(alignment[index:index + window_range])

            orthogroup = split[0].split(".")[0]
            file_name = "{}_block.fasta".format(orthogroup)
            file_path = os.path.join(output_directory, file_name)

            with open(file_path, "w") as h:
                for seq in list_of_seqs:
                    h.write(seq.format("fasta"))

    return output_directory


def longest_exon_length(directory):
    """
    Scans t_coffee alignments in score_ascii format for a region of between 75-2000 positions in length that is
    highly conserved, and sorts by the degree of conservation into an output file
    :param directory: Directory of T_coffee results (containing score_ascii and aln files)
    :return: Dictionary of Orthogroups with a 300bp region TCS scores above 2400
    """
    increments = [75, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350, 375, 400, 450, 500, 550, 600,
                  650, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000]
    increments_rev = increments[::-1]

    dict_of_max_sums = {}

    files = glob.glob(os.path.join(directory, "*.score_ascii"))

    for file in files:
        # Scans an alignment and converts the cons string of numbers into a continous list of numbers
        number_string = ""
        with open(file) as f:
            number_of_specimens = f.read().count(":") - 4
            f.seek(0)
            if number_of_specimens < 5:
                print("Skipping {} Due to Low Specimen Count".format(file))
                continue

            for line in f:
                if line.startswith("cons") and ":" not in line:
                    number = line.rstrip().split(" ")[-1]
                    number_string += number

            number_list = [int(i) for i in number_string]

            # Scans number list for sequence containing the highest window range of conserved bases within 95% of max
            # TCS score for said window range aka 9*Window Range

            # Sort the list so the highest score block within the window range is first. If the window range
            # has 95% quality or higher, add it to dictionary and move on to next file, otherwise decrease
            # window range and try again
            for window_range in increments_rev:

                list_of_sums = []
                if len(number_list) > window_range:
                    for i in range(0, len(number_list) - window_range):
                        the_sum = sum(number_list[i:i + window_range])
                        list_of_sums.append((the_sum, window_range, i))

                    sorted_list = sorted(list_of_sums, reverse=True, key=lambda element: (element[0]))
                    if float(sorted_list[0][0]) >= float(9 * window_range * .95):
                        if os.path.basename(file) not in dict_of_max_sums:
                            dict_of_max_sums[os.path.basename(file)] = sorted_list[0]
                            break

    return dict_of_max_sums


def write_sums(directory, dict_of_max_sums):
    """
    Writes the dictionary of all ortholog T_coffee scores/sums to csv file
    :param directory:
    :param dict_of_max_sums:
    :return:
    """
    if not os.path.exists(directory):
        os.makedirs(directory)

    timestr = time.strftime("%Y%m%d-%H%M%S")
    file_name = "Conserved_Exons_Sums_{}.csv".format(timestr)
    file_path = os.path.join(directory, file_name)

    # Sorts dictionary into a list by score sum and then window length
    sorted_x = sorted(dict_of_max_sums.items(), reverse=True, key=lambda x: (x[1][0], x[1][1]))

    print("Writing T_Coffee score analysis to {}".format(file_path))
    with open(file_path, "w") as f:
        f.write("Orthogroup,Sum,Window,Index\n")
        for entry in sorted_x:
            f.write("{},{},{},{}\n".format(entry[0], entry[1][0], entry[1][1], entry[1][2]))

    return file_path


def filter_blocks(directory, priorities_file, results_dir, window):
    """
    Filters blocks generated by longest exon length and write sum functions based on various criteria
    :param directory: Directory of fasta blocks to filter
    :param priorities_file: Species priorities to focus on
    :param results_dir: Parent Result Folder
    :param window: Minimum length of a conserved block in basepairs
    :return: Output Directory of filtered blocks
    """
    fastas = glob.glob(os.path.join(directory, "*.fasta"))

    output_dir = os.path.join(results_dir, "filtered_blocks_{}".format(window))
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    # Priority list of species to be concerned about in the filtering of sequences
    priorities = {}
    with open(priorities_file) as f:
        lines = f.readlines()
        for line in lines:
            split = line.rstrip().split(",")
            key = split[2]
            primary_priority = int(split[0])
            secondary_priority = int(split[1])
            priorities[key] = (primary_priority, secondary_priority)

    total_seq_length = 0
    total_after_gap_removal = 0
    total_sequences = 0
    gene_count = 0

    # For each block/file extract out sequences that meet the following critiera:
    # Part of Priority List = 1
    # Minimum Length of Window size in basepairs
    # Gaps represent less than 20% of sequence
    # Block contains atleast 5 sequences from priority list = 1

    for fasta in fastas:
        seqs = []
        with open(fasta) as f:
            file_name = os.path.basename(fasta).replace(".fasta", "_filtered.fasta")

            for seq in SeqIO.parse(f, 'fasta'):
                species = " ".join(seq.id.split("_")[1:3])
                if species in priorities:
                    if priorities[species][0] == 1:
                        gaps = seq.seq.count("-")
                        gap_percent = float(gaps / len(seq.seq))

                        if gap_percent > 0.20:
                            pass
                        else:
                            if len(seq.seq) >= window:
                                seqs.append(seq)
        if len(seqs) < 5:
            pass
        else:
            gene_count += 1

            # Chooses sequences based on secondary priority
            # only required if need to curate much more strictly due to budget considerations

            # final_seqs = []
            # for seq in seqs:
            #     species = " ".join(seq.id.split("_")[1:3])
            #     if len(final_seqs) < 5:
            #         if priorities[species][1] == 1:
            #             final_seqs.append(seq)
            # for seq in seqs:
            #     species = " ".join(seq.id.split("_")[1:3])
            #     if len(final_seqs) < 5:
            #         if priorities[species][1] == 2:
            #             final_seqs.append(seq)
            # for seq in seqs:
            #     species = " ".join(seq.id.split("_")[1:3])
            #     if len(final_seqs) < 5:
            #         if priorities[species][1] == 3:
            #             final_seqs.append(seq)
            #
            # if len(final_seqs) == 5:
            #     seqs = final_seqs

            total_sequences += len(seqs)
            for seq in seqs:
                total_seq_length += len(seq.seq)
                seq.seq = seq.seq.ungap(gap="-")
                total_after_gap_removal += len(seq.seq)
            new_file = os.path.join(output_dir, file_name)

            with open(new_file, "w") as g:
                SeqIO.write(seqs, g, "fasta")

    print("Total Genes: {}, "
          "Total Sequences: {}, "
          "Total Length in bp: {}, "
          "After Gap Removal: {}".format(gene_count, total_sequences, total_seq_length, total_after_gap_removal))

    return output_dir


def tile_blocks(directory, results_dir, window):
    """
    Takes a prefiltered block generated by the filtered_blocks function and tiles each bait
    The first 0, 40 or 80 basepairs of each sequence are removed so the baits tile amongst each other
    :param directory:
    :param results_dir:
    :param window:
    :return:
    """
    fastas = glob.glob(os.path.join(directory, "*.fasta"))

    output_dir = os.path.join(results_dir, "tiled_blocks_{}".format(window))
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    for fasta in fastas:
        seqs = []
        with open(fasta) as f:
            count = 0
            for seq in SeqIO.parse(f, 'fasta'):
                seq.description = ""

                # Remove the first 0, 40 or 80 basepairs of the sequence every 3rd time
                count += 1
                if count == 1:
                    pass
                if count == 2:
                    seq.seq = seq.seq[40:]
                if count == 3:
                    seq.seq = seq.seq[80:]
                    count = 0
                seqs.append(seq)

        file_name = os.path.basename(fasta).replace("_block_filtered", "_block_tiled")
        new_file = os.path.join(output_dir, file_name)
        with open(new_file, "w") as g:
            SeqIO.write(seqs, g, "fasta")

    return output_dir


def merge_baits(directory, results_dir, prefix, window):
    """
    Merges multifastas in the input directory into a single multi fasta file. Can be accomplished with bash cat, but
    using biopython ensures each fasta entry is formatted correctly
    :param directory: Input directory of fastas
    :param results_dir: Output Parent directory
    :param prefix: Name of the output file
    :param window:
    :return:
    """
    output_dir = os.path.join(results_dir, "final_baits")
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    fastas = glob.glob(os.path.join(directory, "*.fasta"))
    seqs = []
    for fasta in fastas:
        with open(fasta) as f:
            for seq in SeqIO.parse(f, 'fasta'):
                seq.description = ""
                seqs.append(seq)

    file_name = "{}-{}-final-baits.fasta".format(prefix, window)
    new_file = os.path.join(output_dir, file_name)
    with open(new_file, "w") as g:
        SeqIO.write(seqs, g, "fasta")

    return output_dir


if __name__ == "__main__":
    main()
