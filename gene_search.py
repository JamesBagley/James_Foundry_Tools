# -*- coding: utf-8 -*-
"""
Created on Wed Sep 30 08:33:49 2020

@author: jbag2
"""
import pandas as pd
import plotly.express as px
import numpy as np
from os import listdir
from pathlib import Path
from Bio.Seq import Seq
from collections import namedtuple
from Bio.SeqIO import parse, write
from Bio.SeqRecord import SeqRecord


ROOT = Path(r'C:\Users\jbag2\OneDrive - Concordia University - Canada\Genomes\Genomes and GFFs')
fastas, gffs = [], []
for file in [ROOT.joinpath(x) for x in listdir(ROOT)]:
    #print(file)
    if file.suffix == '.fasta':
        fastas.append(file)
    elif file.suffix == '.gff':
        gffs.append (file)
        
fastas = {x.parts[-1].split(';')[0].strip(): str(x) for x in fastas}
gffs = {x.parts[-1].split('.')[0].split(' ')[-1].strip(): str(x) for x in gffs}

g_and_g_map = pd.DataFrame(fastas.values(),
                           index=fastas.keys(),
                           columns=['fasta']).join(
                               pd.DataFrame(gffs.values(),
                                            index=gffs.keys(),
                                            columns=['gff']))

def find_genes(gene_names):
    cds_names, strains = [], []
    results = dict()
    for gene in gene_names:
        cds_name, strain = correct_name(gene)
        cds_names.append(cds_name)
        strains.append(strain)
    
    if all([strains[0] == strain for strain in strains]):
        
        strain = strains[0]
        gff, genome = load_genome(strain)
        
        sequences = {gene:find_gene(gene, genome, gff) for gene in cds_names}
        results = sequences
        
    else:
        for working_strain in set(strains):
            if working_strain == '04Q':
                continue
            working_set = [] 
            for gene, strain in zip(cds_names, strains):
                if strain != working_strain:
                    continue
                else:
                    working_set.append(gene)
            gff, genome = load_genome(working_strain)
            sequences = {gene:find_gene(gene, genome, gff) for gene in working_set}
            results[working_strain] = sequences
    
    return(results)
        
                               
                               
def find_gene(gene_name, genome, gff):
    
    exons = find_CDS_coords(gene_name, gff)
    exon_sequences = []
    for key, exon in exons.items():
        coords, strand, frame, contig = exon
        fasta_coords = [int(coord) + genome[contig].start_char for coord in coords]
        direction = {'+':+1, '-':-1}[strand]
        sequence = genome.fasta.replace('\n', '')[fasta_coords[0]-1:fasta_coords[1]]
        if direction == -1:
            sequence = str(Seq(sequence[::-1]).complement())
        #exon_sequences.append(sequence[frame:])    
        exon_sequences.append(sequence)    
    sequence = ''.join(exon_sequences)
    return sequence
    

def correct_name(name):
    breakpoints = ['|','c','g']
    locs = [name.find(x) for x in breakpoints]
    subs = [name[:locs[0]], name[locs[1]+1:locs[2]], name[locs[2]+1:]]
    strain = subs[0]
    contig = subs[1]
    gene = subs[2]
    
    return f'{contig.zfill(6)}F|arrow.g{gene}', strain.strip('JB_')

def find_CDS_coords(name, gff):
    found = 0
    exons = {}
    target_attribute = f'ID={name}.t1.cds'
    with open(gff) as file:
        for line in file:
            if line[0]=='#':
                continue
            entries = line.split('\t')
            attributes = entries[8].split(';')
            if target_attribute in attributes:
                if not found:
                    found += 1
                    coords, strand, frame, contig = (entries[3], entries[4]), entries[6], int(entries[7]), entries[0]
                    exons[found] = coords, strand, frame, contig
                else:
                    #multiple exons
                    found += 1
                    #print(f'Found {found} instances of {name} in {gff}')
                    coords, strand, frame, contig = (entries[3], entries[4]), entries[6], int(entries[7]), entries[0]
                    exons[found] = coords, strand, frame, contig
                    
    if found:
        return exons
    else:
        print(f'No instances of "{target_attribute}" found in {gff}')


def load_fasta(path):
    
    Contig_details = namedtuple('Contig', ['name', 'start_row', 'end_row',
                                           'start_char', 'end_char'])
    class Genome():
        # Genome class stores contig headers as keys and returns full contig
        # when indexed by a key
        
        def __init__(self, path):
            contig_names = []
            contig_starts = []
            contig_char_starts = []
            contig_ends = []
            contig_char_ends = []
            char_count = 0
            with open(path) as file :
                for line_num, line in enumerate(file):
                    if line[0] != '>':
                        contig_ends[-1] = line_num
                        char_count += len(line.strip('\n'))
                        contig_char_ends[-1] = char_count
                        continue
                    else:
                        contig_names.append(line[1:].split(' ')[0].strip('\n'))
                        contig_starts.append(line_num+1)
                        contig_ends.append(line_num)
                        char_count += len(line.strip('\n'))
                        contig_char_starts.append(char_count)
                        contig_char_ends.append(char_count)

            self.contigs = \
                {contig_names[i]:
                 Contig_details(contig_names[i], contig_starts[i],
                                contig_ends[i], contig_char_starts[i],
                                contig_char_ends[i])
                 for i in range(len(contig_names))}
            self.path = path
        
        def __getitem__(self, item):
            return self.contigs[item]
        
        @property
        def fasta(self):
            return(open(self.path, 'r').read())
        
    return(Genome(path))

def load_genome(assembly_name):
    locs = g_and_g_map.loc[assembly_name].values
    gff = locs[1]
    genome = load_fasta(locs[0])
    
    return (gff, genome)

def results_to_fasta(results, name):
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.SeqIO import write
    
    reformatted = []
    for strain in results.keys():
        for gene in results[strain].keys():
            reformatted.append(['|'.join([strain, gene]), results[strain][gene]])
    out = [SeqRecord(Seq(sequence), id=full_name) for full_name, sequence in reformatted]
    
    write(out, f'{name}.fasta', 'fasta')
    
def from_aybrah(FOG, amino_acid=True):
    aybrah = pd.read_excel(r'C:\Users\jbag2\OneDrive - Concordia University - Canada\Aybrah annotation dataset  (2).xlsx')
    aybrah = aybrah.set_index('FOG')
    SOI = 'JB_02W,JB_04R,JB_02M,pku_NG7,JB_02G,pfe_madison,JB_02L,JB_02Q,JB_019,JB_01V,JB_01X,pnor_UCD,sce,opm,ppa,dbx,yli'.split(',')
    GOI = aybrah.loc[FOG].loc[SOI]
    
    JB = 'JB_02W,JB_04R,JB_02M,JB_02G,JB_02L,JB_02Q,JB_019,JB_01V,JB_01X'.split(',')
    AYB = 'sce,opm,ppa,dbx,yli'.split(',')
    unimplemented = set(GOI.index).difference(set(AYB)).difference(set(JB))
    print(f'strains {", ".join(unimplemented)} NOT YET IMPLEMENTED')
    
    JB_genes = find_genes(GOI.loc[JB].str.split(';', expand=True).melt().dropna().value.values)
    ayb_fasta = parse(r'C:\Users\jbag2\OmicsBoxWorkspace\yeast_proteome_aybrah.faa', 'fasta')
    ayb_goi = GOI.loc[AYB].str.split(';', expand=True).melt().dropna().value.values
    ayb_out = []
    for seq in ayb_fasta:
        if seq.id in ayb_goi:
            ayb_out.append([seq.id, seq.seq])
    JB_out = []
    for strain in JB_genes.keys():
            try:
                for gene in JB_genes[strain].keys():
                    if amino_acid:
                        JB_out.append(['|'.join([strain, gene]), Seq(JB_genes[strain][gene]).translate()])
                    else:
                        JB_out.append(['|'.join([strain, gene]), Seq(JB_genes[strain][gene])])
                                      
            except(AttributeError):
                continue
    
    out = JB_out + ayb_out
    return out
    