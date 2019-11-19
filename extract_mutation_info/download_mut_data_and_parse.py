#!/usr/bin/env python2
# coding: utf-8


import sys, getopt, argparse
from libs.textMineDescriptions import ParseMutInfo
from libs.UniprotXmlParser import GetMutationsData



def main():
    parser = argparse.ArgumentParser(description="Script to download mutation information data from uniprot and identify regular expressions of interest.")
    parser.add_argument('-outfileA', help='outfileA: Path to outfile containing all mutations')
    parser.add_argument('-outfileP', help='outfileA: Path to outfile containing parsed mutations')
    parser.add_argument('-infile', help='Infile: list of all uniprot proteins for a given specie', type=str)
    args = parser.parse_args()

    print 'Retrieving mutation features from UNIPROT site ...'
    (name2accession,accession2name,accession2descriptions) = GetMutationsData(args.outfileA, args.infile)

    print ''
    print 'Extracting useful information from mutation description field ...'
    ParseMutInfo(name2accession,accession2name,accession2descriptions, args.outfileP)


if __name__ == '__main__':
    main()

