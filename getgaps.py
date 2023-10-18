#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File    :   getgaps.py
@Time    :   2029/09/01
@Author  :   Xu Wang, Qing Zhang
@Version :   1.0
@Contact :   wangxu2018@cau.edu.cn
@License :   (C)Copyright 2021-2022, CAAS ShenZhen
@Desc    :   Get the information of gaps with gff format
@Usage   :   getgaps.py FILENAME.fasta > FILENAME.gff3
'''
import argparse
import re
from Bio import SeqIO

# Parse command-line arguments
parser = argparse.ArgumentParser()
parser.add_argument("fasta")
args = parser.parse_args()

# Open FASTA, search for masked regions, print in GFF3 format
with open(args.fasta) as handle:
    i = 0
    for record in SeqIO.parse(handle, "fasta"):
        for match in re.finditer('N+', str(record.seq)):
            i = i+1
            print (record.id, ".", "gap", match.start() + 1, match.end(), ".", ".", ".", "Name=gap" + str(i) + ";size=" + str(match.end()-match.start()), sep='\t')

# Done
