# -*- coding: utf-8 -*-
"""
Created on Sun Aug 15 14:54:22 2021

@author: 10979
"""

from Bio import SeqIO
from Bio import Seq
import regex as re
import pandas as pd
import numpy as np


def lncRNA_features(fasta):
    records = SeqIO.parse(fasta, 'fasta')
    orf_length = []
    orf_count = []
    orf_position = []
    ID = []
    for record in records:
        for strand, seq in (1, record.seq), (-1, record.seq.reverse_complement()):
            for frame in range(3):
                length = 3 * ((len(seq)-frame) // 3)
                for pro in seq[frame:frame+length].translate(table = 1).split("*")[:-1]:
                    if 'M' in pro:
                        orf = pro[pro.find('M'):]
                        pos = seq[frame:frame+length].translate(table=1).find(orf)*3 + frame +1
                        orf_length.append(len(orf)*3+3)
                        orf_count.append(frame)
                        orf_position.append(pos)
                        ID.append(record.id)
                    else:
                        orf_length.append(0)
                        orf_count.append(0)
                        orf_position.append(0)
                        ID.append(record.id)
    data = {
        'ID':ID,
        'orf_length':orf_length,
        'orf_count':orf_count,
        'orf_position':orf_position}
    df = pd.DataFrame(data)
    df.sort_values(by=["ID","orf_length"],ascending=[False,False],inplace=True)
    df.duplicated(['ID'])      
    df2=df.drop_duplicates(['ID'])
    ORF = df2
    print('ORF features over.')
    
    transcript_length = []
    start_codon_number = []
    end_codon_number = []
    GC_pro = []
    ID=[]

    for record in SeqIO.parse(fasta, 'fasta'):
        ID.append(record.id)
        record = record.upper()
        transcript_length.append(len(record))
        start_codon_number.append(record.seq.count('ATG'))
        end_codon_number.append(record.seq.count('TAG')+record.seq.count('TAA')+record.seq.count('TGA'))
        GC_pro.append(100*float(record.seq.count('G')+record.seq.count('C'))/len(record))

    data = {
        'ID':ID,
        'transcript_length':transcript_length,
        'start_codon_number':start_codon_number,
        'end_codon_number':end_codon_number,
        'GC_pro':GC_pro
    }
    df = pd.DataFrame(data)
    codon = df
    print('codon features over.')
    All = pd.merge(codon,ORF,how='left')
    A = All.fillna(0)
    return A


