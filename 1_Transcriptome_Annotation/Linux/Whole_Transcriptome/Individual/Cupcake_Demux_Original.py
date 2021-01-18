#!/usr/bin/env python
# Szi Kay Leung
# 07/01/2020: Create fake flnc_report.csv for demultiplex https://github.com/Magdoll/cDNA_Cupcake/wiki/Tutorial:-Demultiplexing-SMRT-Link-Iso-Seq-Jobs
# script.py <inputreadstats> <outputfile>

# load packages
import sys

# Input files
inputreadstats = sys.argv[1]            # input collapsed read stats file
outputfile = sys.argv[2]                # output file (fake flnc read)

d = {'m54082_180607_173058': 'Q21_WT',
     'm54082_180605_141944': 'O18_TG',
     'm54082_190306_083150': 'L22_TG',
     'm54082_190307_045507': 'K18_TG',
     'm54082_190401_165425': 'O23_WT',
     'm54082_190403_135102': 'S23_WT',
     'm54082_190404_101400': 'S18_TG',
     'm54082_190405_063832': 'K17_WT',
     'm54082_190430_163756': 'M21_WT',
     'm54082_190524_145911': 'K23_WT',
     'm54082_190527_173356': 'Q20_TG',
     'm54082_190529_082942': 'K24_TG'
   }

h = open(outputfile, 'w')
h.write("id,primer\n")

f = open(inputreadstats)
f.readline()
for line in f:
    seqid = line.strip().split()[0]
    movie = seqid.split('/')[0]
    h.write("{0},{1}\n".format(seqid, d[movie]))
h.close()
