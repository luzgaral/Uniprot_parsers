#!/usr/bin/env python2
# coding: utf-8


import sys, getopt, argparse
from libs.textMineDescriptions import PareseBgeeExpression
from libs.UniprotXmlParser import GetBgeeExpression



def main():
    parser = argparse.ArgumentParser(description="Script to download active sites from uniprot.")
    parser.add_argument('-outfileP', help='Path to outfile containing active sites')
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-uniprots', help='Protein identifiers. Separated by ",". Example P15056,Q9P0W2.', type=str)
    group.add_argument('-infile', help='Infile: list of uniprot proteins', type=str)
    args = parser.parse_args()

    print 'Retrieving descriptions from UNIPROT site ...'
    (accession2name,accession2descriptions) = GetBgeeExpression(args.uniprots, args.infile)


    print 'Extracting Bgee expression information ...'
    PareseBgeeExpression(accession2name,accession2descriptions, args.outfileP)


if __name__ == '__main__':
    main()

