import pandas as pd
import glob
import sys
import subprocess as sp
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import os
import numpy as np

three_to_one = {
    'ALA': 'A',
    'ARG': 'R',
    'ASN': 'N',
    'ASP': 'D',
    'CYS': 'C',
    'GLU': 'E',
    'GLN': 'Q',
    'GLY': 'G',
    'HIS': 'H',
    'ILE': 'I',
    'LEU': 'L',
    'LYS': 'K',
    'MET': 'M',
    'PHE': 'F',
    'PRO': 'P',
    'SER': 'S',
    'THR': 'T',
    'TRP': 'W',
    'TYR': 'Y',
    'VAL': 'V'
}

def run_NACCESS(pdb):
    pdbid_chain = pdb.split('.')[0]
    sp.call(f'naccess {pdb}', shell = True)
    return f'{pdbid_chain}.rsa'


def find_hits(pdb_chain):
    df = pd.read_csv(f'./{pdb_chain}_BLAST_HITS.csv',header=None)
    filtered_df = df[(df[1].str.len() <= 6) & (df[1] != pdb)]
    num_rows = min(len(filtered_df), 10)

    formatted_vals = []
    values = [filtered_df.iloc[i, 1] for i in range(num_rows)]
    for val in values:
        split = val.split('_')
        formatted_vals.append(f'{split[0]} {split[1]}')

    return formatted_vals


def run_clean(hit):
    split = hit.split()
    pdbid_chain = f'{split[0]}_{split[1]}'
    sp.call(f'python2 ~/rosetta/tools/protein_tools/scripts/clean_pdb.py {hit}', shell = True)
    return f'{pdbid_chain}.pdb', f'{pdbid_chain}.fasta'

def run_alignment(pdb,target_fasta, hit_fasta):
    seq1 = ''
    with open(target_fasta,'r') as infile:
        for line in infile.readlines()[1:]:
            seq1 = seq1.join(line.strip())

    seq2 = ''

    with open(hit_fasta,'r') as infile:
        for line in infile.readlines()[1:]:
            seq2 = seq2.join(line.strip())


    alignments = pairwise2.align.globalxx(seq1, seq2)

    alignment = alignments[0]
    print(format_alignment(*alignment))
    target_aligned = alignment.seqA
    query_aligned = alignment.seqB
    alignment_dict = {}

    target_pos = 1
    query_pos = 1

    for t_res, q_res in zip(target_aligned, query_aligned):
        if t_res != '-' and q_res != '-':
            if t_res == q_res:
                alignment_dict[target_pos] = query_pos
            else:
                print('Mismatch. This should not happen.')
        if t_res != '-':
            target_pos += 1
        if q_res != '-':
            query_pos += 1
    return alignment_dict

def assign_se(relSASA, polar_frac):
    if relSASA >= 36.0:
        return 'E'
    if 9.0 <= relSASA < 36.0:
        if polar_frac < 0.67:
            return 'P1'
        if polar_frac >= 0.67:
            return 'P2'
    if relSASA < 9.0:
        if polar_frac < 0.45:
            return 'B1'
        if 0.45 <= polar_frac < 0.58:
            return 'B2'
        if polar_frac >= 0.58:
            return 'B3'


def remap_data(res_sanity_check, alignment, hit_SASA, id_chain):
    mapped = []
    for resno, abbrv in res_sanity_check.items():
        try:
            mapped_no = alignment[resno]
        except KeyError:
            mapped.append('X')
            continue
        with open(hit_SASA,'r') as infile:
            for line in infile.readlines()[4:]:
                spline = line.split()
                if spline[3] == str(mapped_no):
                    if abbrv == spline[1]:
                        ts_abs = float(spline[6])
                        ap_abs = float(spline[-2])
                        try:
                            polar_frac = ap_abs / ts_abs
                        except ZeroDivisionError:
                            polar_frac = 0.0
                        rel = float(spline[7])
                        mapped.append(assign_se(rel, polar_frac))
                        break
                    else:
                        mapped.append('ERR')
                        break
    return mapped

def fetch_seq(fasta):
    seq = ''
    with open(fasta,'r') as infile:
        lines = infile.readlines()
        for line in lines[1:]:
            seq += line.strip()
    return seq

def main():
    pdbid, chain = sys.argv[1],sys.argv[2]

    pdb_fasta = f'./{pdbid}_{chain}.fasta'
    seq = fetch_seq(pdb_fasta)

    res_dict = {}
    res = []

    for i,c in enumerate(seq):
        for k, v in three_to_one.items():
            if v == c:
                res_dict[i + 1] = k
        res.append(i+1)

    df = pd.DataFrame(columns=res)
    hits = find_hits(f'{pdbid}_{chain}')
    for hit in hits:
        hit_id, hit_fasta = run_clean(hit)
        hit_SASA = run_NACCESS(hit_id)
        alignment = run_alignment(pdb_chain,pdb_fasta, hit_fasta)
        mapped = remap_data(res_dict, alignment, hit_SASA, pdb_chain)
        df.loc[len(df)] = mapped
    df.to_excel(f'{pdbid}_{chain}_mapped.xlsx')
main()

