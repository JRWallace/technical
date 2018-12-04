#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 30 13:26:05 2018

@author: joselynn
"""
import subprocess
import argparse
#import sys

#sys.argv.extend(['-f', 'filename'])

parser = argparse.ArgumentParser()

parser.add_argument('-f', nargs=1, required=True, help='forward fastq', dest='forward')
parser.add_argument('-r', nargs=1, required=True, help='reverse fastq', dest='reverse')
parser.add_argument('-x', nargs=1, required=True, help='reference fasta', dest='reference')
parser.add_argument('-q', type=int, choices=range(0, 41), nargs=1, required=True, help='trimmomatic mean read quality cutoff', dest='quality')
parser.add_argument('-b', nargs=1, required=True, help='trimmomatic output base file name', dest='baseout')
    
#arguments=parser.parse_args()

args=parser.parse_args(['-f=SRR961514_1.fastq', '-r=SRR961514_2.fastq', '-x=sequence.fasta', '-q=30', '-b=out'])

#print(args.forward[0])

#f=arguments.forward
#r=arguments.reverse
#x=arguments.reference
#q=arguments.quality
#b=arguments.baseout

#print (["trimmomatic PE -threads 8 ", arguments.forward," ", arguments.reverse, "baseout ", arguments.baseout, " AVGQUAL:", arguments.quality])

#subprocess.call(["trimmomatic PE -threads 8 ", arguments.forward, " ", arguments.reverse, "baseout ", arguments.baseout, " AVGQUAL:", arguments.quality])
proc = subprocess.run(["trimmomatic", "PE", "-threads", "8", args.forward[0], args.reverse[0], "-baseout", args.baseout[0], "AVGQUAL:" + str(args.quality[0])])

print(proc.args)


#, stdout=stdout, shell=True)

#stdout = open("arguments.baseout.txt","wb")