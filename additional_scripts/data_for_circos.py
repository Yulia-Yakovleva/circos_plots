import re
import argparse
import numpy as np
import pandas as pd
from Bio import SeqIO


def make_karyotype(rec, out_f):
    with open(out_f, 'w') as out:
        print(f"Writing...\t{out_f}")
        line = f"chr\t-\t{rec.id}\t{rec.id}\t0\t{len(rec.seq)}\tlgrey"
        out.write(line)


def make_rRNA(in_f, out_c, out_l):
    with open(in_f, 'r') as input, open(out_c, 'w') as output_c, open(out_l, 'w') as output_l:
        print(f"Writing...\t{out_c}, {out_l}")
        for in_line in input:
            if in_line.startswith('#'):
                pass
            else:
                in_line = in_line.strip()
                cols = in_line.split('\t')
                if cols[2] == 'rRNA':
                    outc_line = f"{cols[0]}\t{cols[3]}\t{cols[4]}\tcolor=red\n"
                    outl_line = f"{cols[0]}\t{cols[3]}\t{cols[4]}\trRNAs\n"
                    output_c.write(outc_line)
                    output_l.write(outl_line)


def GC_skew_window(s):
    g = s.count('G') + s.count('g')
    c = s.count('C') + s.count('c')
    try:
        skew = (g - c) / float(g + c)
    except ZeroDivisionError:
        skew = 0
    return round(skew, 4)


def make_gc_skew(rec, out_f, window=10000, step=1000):
    with open(out_f, 'w') as out:
        print(f"Writing...\t{out_f}")
        for i in range(0, len(rec), step):
            subseq = rec.seq[i:i + window]
            gc_skew = GC_skew_window(subseq)
            start = (i + 1 if (i + 1 <= len(rec)) else i)
            end = (i + step if (i + step <= len(rec)) else len(rec))
            if gc_skew > 0:
                line = f"{rec.id}\t{start}\t{end}\t{gc_skew}\tfill_color=orange\n"
                out.write(line)
            else:
                line = f"{rec.id}\t{start}\t{end}\t{gc_skew}\tfill_color=blue\n"
                out.write(line)


def get_is_intervals(in_f):
    df = pd.read_csv(in_f)
    intervals = []
    for orf_L, orf_R in zip(df.orf_L, df.orf_R):
        if orf_L < orf_R:
            start = orf_R
            end = orf_R
        else:
            start = orf_R
            end = orf_L
        intervals.append((start, end))
    return tuple(intervals)


def split_genome_on_intervals(rec, interval):
    intervals = []
    genome_size = len(rec.seq)
    for i in range(0, genome_size, interval):
        intervals.append((i + 1, i + interval))
    return tuple(intervals)


def count_intervals(genome_interval, feat_intervals):
    passed_intervals = 0
    failed_intervals = 0
    g_range = range(genome_interval[0], genome_interval[1])
    for interval in feat_intervals:
        start = interval[0]
        end = interval[1]
        if start in g_range or end in g_range:
            passed_intervals += 1
        elif start not in g_range and end not in g_range:
            failed_intervals += 1
    return passed_intervals, failed_intervals


def intersect_intervals(g_intervals, f_intervals, rec_id):
    strings = []
    check_results = {}
    counts = 0
    for g_interval in g_intervals:
        counts += 1
        passed, failed = count_intervals(g_interval, f_intervals)
        check_results[counts] = passed
    for g_int, count in zip(g_intervals, check_results):
        strings.append(f"{rec_id}\t{g_int[0]}\t{g_int[1]}\t{check_results[count]}\n")
    return strings


def write_is_counts(rec, df, out_prefix, interval):
    outf_name = f"{out_prefix}_{interval // 1000}kb.txt"
    with open(df, 'r') as is_seqs, open(outf_name, 'w') as out_f:
        print(f"Writing...\t{outf_name}")
        genome_intervals = split_genome_on_intervals(rec, interval)
        insertion_seqs_intervals = get_is_intervals(is_seqs)
        out_f.writelines(intersect_intervals(genome_intervals, insertion_seqs_intervals, rec.id))


# here starts the part of parse_to_df.py script
p = "([A-Za-z0-9_\(\)\/\s]+)\.([A-Za-z0-9_\.]+):(\d+)-(\d+) ([+|-]).*"
pattern = re.compile(p)
columns = ["block", "species", "chr", "chr_beg", "chr_end", "orientation"]


def find_indices(lst, condition):
    return [i for i, elem in enumerate(lst) if condition(elem)]


def parse_to_df(file_name):
    with open(file_name) as f:
        lines = f.readlines()
    last_line = len(lines) - 1
    while lines[last_line] == '\n': last_line -= 1
    n_at_end = len(lines) - 1 - last_line
    for _ in range(1 - n_at_end): lines.append('\n')
    bs = np.split(lines, find_indices(lines, lambda x: x[0] == ">"))
    temp = []
    for i, b in enumerate(bs):
        if len(b) == 0: continue
        b_i = int(b[0][1:])
        for oc in b[1:-1]:
            m = pattern.match(oc)
            temp.append([b_i, m.group(1), m.group(2), int(m.group(3)), int(m.group(4)), m.group(5)])
    return pd.DataFrame(temp, columns=columns)


# here ends the part of parse_to_df.py script


def filter_common_blocks(df):
    allowed_blocks = set()
    all_sp = len(df['species'].unique())
    for block, df_block in df.groupby('block'):
        if len(df_block) == len(df_block['species'].unique()) == all_sp:
            allowed_blocks.add(block)
    return df[df.apply(lambda x: x['block'] in allowed_blocks, axis=1)]


def get_blocks_intervals(df):
    intervals = []
    for start, end in zip(df.chr_beg, df.chr_end):
        intervals.append((start, end))
    return tuple(intervals)


def write_blocks_counts(rec, in_f, out_postfix, btype, acc_id, window=5000):
    outf_name = f"{btype}_{out_postfix}"
    with open(outf_name, 'w') as out_f:
        print(f"Writing...\t{outf_name}")
        df = parse_to_df(in_f)
        if btype == 'all_founded':
            df = df[df.species == acc_id].sort_values(by=['chr_beg'])
        elif btype == 'common':
            df = filter_common_blocks(df)[filter_common_blocks(df).species == acc_id].sort_values(by=['chr_beg'])
        else:
            raise ValueError('Key "type" could be only "common" or "all_founded"')
        genome_intervals = split_genome_on_intervals(rec, window)
        blocks_intervals = get_blocks_intervals(df)
        out_f.writelines(intersect_intervals(genome_intervals, blocks_intervals, rec.id))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Create files for Circos plot')
    parser.add_argument('-chr', '--chr', type=str,
                        help='Path to genome sequence in .fasta format. First sequence will be recognized as chromosome')
    parser.add_argument('-gw', '--genome_window', type=int, help='Window size (in bp) to split genome on intervals')
    parser.add_argument('-gff', '--gff', type=str, help='Path to genome annotation in .gff format')
    parser.add_argument('-is', '--is_data', type=str, help='Path to IS table in .csv format')
    parser.add_argument('-b', '--blocks', type=str, help='Path to blocks_coords.infercars file')
    parser.add_argument('-id', '--assembly_id', type=str,
                        help='GenBankID of assembly for filter blocks data. Example: GCA_001578125')
    parser.add_argument('-k', '--karyotype', type=str, default='karyotype.txt',
                        help='Path to output karyotype.txt file')
    parser.add_argument('-rc', '--rrna_coords', type=str, default='rRNA_coords.txt',
                        help='Path to output rRNA_coords.txt file')
    parser.add_argument('-rl', '--rrna_labels', type=str, default='rRNAs_labels.txt',
                        help='Path to output rRNAs_labels.txt file')
    parser.add_argument('-gc', '--gc_skew', type=str, default='gc_skew.txt', help='Path to output gc_skew.txt file')
    parser.add_argument('-isc', '--is_counts', type=str, default='insertion_sequences_counts',
                        help='Prefix for insertion_sequences_counts_*kb.txt output file')
    parser.add_argument('-bc', '--blocks_counts', default='blocks_counts.txt',
                        help='Postfix for all founded and common blocks counts file')
    args = parser.parse_args()

    chromosome = [rec for rec in SeqIO.parse(args.chr, 'fasta')][0]

    make_karyotype(chromosome, args.karyotype)
    make_rRNA(args.gff, args.rrna_coords, args.rrna_labels)
    make_gc_skew(chromosome, args.gc_skew)
    write_is_counts(chromosome, args.is_data, args.is_counts, interval=1000)
    write_is_counts(chromosome, args.is_data, args.is_counts, interval=10000)
    write_blocks_counts(chromosome, args.blocks, args.blocks_counts, 'all_founded', args.assembly_id, args.genome_window)
    write_blocks_counts(chromosome, args.blocks, args.blocks_counts, 'common', args.assembly_id, args.genome_window)
