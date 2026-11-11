#!/usr/bin/env python3

from argparse import ArgumentParser
from Bio import Entrez
from tqdm import tqdm

def parse_args():
    parser = ArgumentParser()
    parser.add_argument('--infile', required=True)
    parser.add_argument('--outfile', required=True)
    parser.add_argument('--email', required=True)
    return parser.parse_args()

def get_accessions(infile):
    return [line.strip().split(' ')[0][1:] for line in open(infile) if '>' in line]

def gen_chunk_list(accessions, n):
    return [accessions[i:i + n] for i in range(0, len(accessions), n)]

def get_sources(accessions):
    acc_to_source = {}
    try:
        handle = Entrez.efetch(db='nucleotide', rettype='gb', id=','.join(accessions))
        current = None
        for line in handle:
            if line.startswith('LOCUS'):
                current = line.split()[1].split('.')[0]   # accession, version stripped
                acc_to_source[current] = 'NA'             # placeholder until we see its SOURCE
            elif current and 'SOURCE' in line:
                acc_to_source[current] = line.strip().split('     ')[-1][1:]
        handle.close()
    except Exception:
        pass
    return [acc_to_source.get(a.split('.')[0], 'NA') for a in accessions]

def process_chunks(accession_chunk_list):
    results = []
    [results.extend(el) for el in (get_sources(chunk) for chunk in tqdm(accession_chunk_list))]
    return results

def write_out(accessions, sources, outfile):
    with open(outfile, 'w') as o:
        [o.write(f'{a}\t{s}\n') for a,s in zip(*[accessions, sources])]

def main():
    args = parse_args()
    Entrez.email = args.email
    accessions = get_accessions(args.infile)
    accession_chunk_list = gen_chunk_list(accessions, 1000)
    sources = process_chunks(accession_chunk_list)
    write_out(accessions, sources, args.outfile)

if __name__ == '__main__':
    main()
