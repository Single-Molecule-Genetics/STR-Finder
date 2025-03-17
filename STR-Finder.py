import re
import argparse
import pickle
from pyfasta import Fasta
import os


# Function to parse command-line arguments
def parse_arguments():
    parser = argparse.ArgumentParser(description="Find STRs in a given genome sequence.")
    parser.add_argument("-g", "--genome", required=True, help="Path to the genome FASTA file")
    parser.add_argument("-o", "--output", required=True, help="Output path for the .pickle file")
    parser.add_argument("--min_mono", type=int, default=5, help="Minimum length for mononucleotide STRs (default: 5)")
    parser.add_argument("--min_di", type=int, default=6, help="Minimum length for dinucleotide STRs (default: 6)")
    parser.add_argument("--min_tri", type=int, default=4, help="Minimum length for trinucleotide STRs (default: 4)")
    return parser.parse_args()


# Regex patterns for STRs
def generate_str_regex(min_mono, min_di, min_tri):
    mono_regex = rf"(?:A{{{min_mono},}}|T{{{min_mono},}}|C{{{min_mono},}}|G{{{min_mono},}}|a{{{min_mono},}}|t{{{min_mono},}}|c{{{min_mono},}}|g{{{min_mono},}})"
    di_regex = rf"((?:([ACGTacgt])(?!\2)[ACGTacgt]))(?:\1){{{min_di - 1},}}"
    tri_regex = rf"([ACGT]{{3}})\1{{{min_tri - 1},}}"
    return mono_regex, di_regex, tri_regex


# Find repeats in the sequence
def find_repeats(query_seq, chromosome, repeat_type, repeat_integer, min_mono, min_di, min_tri, ref_str,
                 str_region_annotation, chromosome_STRdict):
    mono_regex, di_regex, tri_regex = generate_str_regex(min_mono, min_di, min_tri)

    if repeat_type == 'mono':
        matches = [repeat for repeat in re.finditer(mono_regex, query_seq.upper()) if len(repeat.group()) >= min_mono]
    elif repeat_type == 'di':
        matches = [repeat for repeat in re.finditer(di_regex, query_seq.upper())]
    elif repeat_type == 'tri':
        matches = [repeat for repeat in re.finditer(tri_regex, query_seq.upper())]

    for match in matches:
        repeat_start = match.start()
        repeat_end = match.end()
        repeat = match.group()[:repeat_integer].upper()

        if len(set(repeat)) == 1 and repeat_type != 'mono':
            continue

        repeat_length = len(match.group()) // repeat_integer
        chr_key = range(repeat_start, repeat_end + 1)
        ref_str[chromosome].append(chr_key)
        annotation = f"{repeat_start}-{repeat_end}-{repeat_length}{repeat}"
        str_region_annotation[chromosome][chr_key] = annotation
        annotation_strcalling = f"{chromosome}-{repeat_start}-{repeat_end - 1}-{repeat_length}{repeat}"
        chromosome_STRdict[chromosome].update({num: annotation_strcalling for num in chr_key})


# Main function
def main():
    args = parse_arguments()

    # Load genome
    if not os.path.exists(args.genome):
        raise FileNotFoundError(f"Genome file not found: {args.genome}")

    genome = Fasta(args.genome)
    ref_info = {i[0].split(' ')[0]: {'seq': str(i[1])} for i in genome.items()}

    # Data structures
    ref_str = {i: [] for i in ref_info}
    str_region_annotation = {i: {} for i in ref_info}
    chromosome_STRdict = {i: {} for i in ref_info}

    # Find STRs
    for repeat_type in ['mono', 'di', 'tri']:
        repeat_integer = 1 if repeat_type == 'mono' else 2 if repeat_type == 'di' else 3
        for chromosome in ref_str.keys():
            find_repeats(ref_info[chromosome]['seq'], chromosome, repeat_type, repeat_integer,
                         args.min_mono, args.min_di, args.min_tri, ref_str, str_region_annotation, chromosome_STRdict)

    # Save results
    output_filename = os.path.join(args.output, f"CuteCV_STRs_{args.min_mono}m_{args.min_di}d_{args.min_tri}t.pickle")
    data_to_save = {
        'ref_str': ref_str,
        'ref_info': ref_info,
        'chromosome_STRdict': chromosome_STRdict,
        'str_region_annotation': str_region_annotation
    }

    with open(output_filename, "wb") as f:
        pickle.dump(data_to_save, f)

    print(f"Results saved to {output_filename}")


if __name__ == "__main__":
    main()