# -*- coding: utf-8 -*-
#!/usr/bin/env python

from Bio import SeqIO
from collections import defaultdict
import argparse
import re

def main(GenomeFasta, AnnotationGTF, GenesToIsolate, ReverseMinusStrandGenes = False, OutputPrefix = None, Quiet = False):

    print GenesToIsolate
    if OutputPrefix:
        fastaFilepathOut = OutputPrefix + '.fa'
        gtfFilepathOut = OutputPrefix + '.gtf'
    else:
        fastaFilepathOut = '_'.join(GenesToIsolate) + '.fa'
        gtfFilepathOut = '_'.join(GenesToIsolate) + '.gtf'
    print fastaFilepathOut, gtfFilepathOut


    SwitchSignDict = {'+':'-', '-':'+'}
    
    #First pass iterate through annotation file, locate lines that contain a match to any of the genes, and save the chromosome, start, stop, and strand to a variable. This is important to find the min and max coordinates that match the features of a gene so I know the full extent of the gene.
    ChromosomesAndCoordinatesForGenesToIsolate = defaultdict(lambda: defaultdict(set))
    SearchPattern = re.compile('|'.join(GenesToIsolate), re.IGNORECASE)
    with open(AnnotationGTF, 'rU') as AnnotationGTF_fh:
        for line in AnnotationGTF_fh:
            GeneMatch = SearchPattern.search(line)
            if GeneMatch:
                chrom, source, feature, start, end, score, strand, frame, attributes = line.strip('\n').split('\t')
                ChromosomesAndCoordinatesForGenesToIsolate[GeneMatch.group()]['chrom'].add(chrom)
                ChromosomesAndCoordinatesForGenesToIsolate[GeneMatch.group()]['strand'].add(strand)
                ChromosomesAndCoordinatesForGenesToIsolate[GeneMatch.group()]['coordinates'].update((int(start), int(end)))
    FastaEntriesToWriteOut=[]
    ReferenceGenome = SeqIO.index(GenomeFasta, "fasta")
    for gene in sorted(GenesToIsolate): #Check that each gene was found on one and only one chromosome, and only on one strand/
        if gene not in ChromosomesAndCoordinatesForGenesToIsolate:
            print('ERROR. "{}" not found in input annotation file. Check that this gene (an exact string match) is in the annotation file or run script without "{}" in the --GenesToIsolate parameter. Aborted script.').format(gene, gene)
            return None
        elif len(ChromosomesAndCoordinatesForGenesToIsolate[gene]['chrom']) != 1:
            print('ERROR. "{}" found on {} chromosomes in input annotation file. Must be on one and only one chromosome. Aborted script.').format(gene, len(ChromosomesAndCoordinatesForGenesToIsolate[gene]['chrom']))
            return None
        elif len(ChromosomesAndCoordinatesForGenesToIsolate[gene]['strand']) != 1:
            print('ERROR. "{}" found on both + and - strand in input annotation file. Must be on + or - strand only. Aborted script.').format(gene)
            return None
        else: #If all the above conditions look good, add the gene sequence to the lits of Fasta entries to write out
            chrom = next(iter(ChromosomesAndCoordinatesForGenesToIsolate[gene]['chrom']))
            strand = next(iter(ChromosomesAndCoordinatesForGenesToIsolate[gene]['strand']))
            GeneCoordinatesInReferenceGenome = (min(ChromosomesAndCoordinatesForGenesToIsolate[gene]['coordinates'])-1, max(ChromosomesAndCoordinatesForGenesToIsolate[gene]['coordinates']))
            ChromosomesAndCoordinatesForGenesToIsolate[gene]['coordinatesmin'] = GeneCoordinatesInReferenceGenome[0] + 1
            ChromosomesAndCoordinatesForGenesToIsolate[gene]['coordinatesmax'] = GeneCoordinatesInReferenceGenome[1]

        if ReverseMinusStrandGenes and strand == '-':
            ChromosomesAndCoordinatesForGenesToIsolate[gene]['CoordinateSystem'] = 'UseNegativeStrandCoordinates'
            FastaEntryOut = ReferenceGenome[chrom][GeneCoordinatesInReferenceGenome[0]:GeneCoordinatesInReferenceGenome[1]].reverse_complement()
        else:
            ChromosomesAndCoordinatesForGenesToIsolate[gene]['CoordinateSystem'] = 'UsePositiveStrandCoordinates'
            FastaEntryOut = ReferenceGenome[chrom][GeneCoordinatesInReferenceGenome[0]:GeneCoordinatesInReferenceGenome[1]]
        FastaEntryOut.id = gene
        FastaEntryOut.description = ''
        FastaEntriesToWriteOut.append(FastaEntryOut)
    SeqIO.write(FastaEntriesToWriteOut, fastaFilepathOut, "fasta")

    AnnotationGTF_fh = open(AnnotationGTF, 'rU')
    GtfLinesOut = []
    for line in AnnotationGTF_fh:
        GeneMatch = SearchPattern.search(line)
        if GeneMatch:
            gene = GeneMatch.group()
            chrom, source, feature, start, end, score, strand, frame, attributes = line.strip('\n').split('\t')
            if ChromosomesAndCoordinatesForGenesToIsolate[gene]['CoordinateSystem'] == 'UseNegativeStrandCoordinates':
                StartOut = len(ReferenceGenome[chrom]) - int(end) + 1 - (len(ReferenceGenome[chrom]) - ChromosomesAndCoordinatesForGenesToIsolate[gene]['coordinatesmax'])
                EndOut = len(ReferenceGenome[chrom]) - int(start) + 1 - (len(ReferenceGenome[chrom]) - ChromosomesAndCoordinatesForGenesToIsolate[gene]['coordinatesmax'])
                strand = '+'
            else:
                StartOut = int(start) - ChromosomesAndCoordinatesForGenesToIsolate[gene]['coordinatesmin'] + 1
                EndOut = int(end) - ChromosomesAndCoordinatesForGenesToIsolate[gene]['coordinatesmin'] + 1
            GtfLinesOut.append([gene, source, feature, str(StartOut), str(EndOut), score, strand, frame, attributes])
    AnnotationGTF_fh.close()
    GtfLinesOut.sort(key=lambda k: (k[0], k[3], k[4]))
    AnnotationGTF_out_fh = open(gtfFilepathOut, 'w')
    AnnotationGTF_out_fh.write('\n'.join(['\t'.join(line) for line in GtfLinesOut]))
    AnnotationGTF_out_fh.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-G','--GenomeFasta', metavar = "[filepath]", help="REQUIRED. The input genome file [.fasta]", required=True)
    parser.add_argument('-A', '--AnnotationGTF', metavar = "[filepath]", help="REQUIRED. The input annotation file [.gtf]", required=True)
    parser.add_argument('-genes', '--GenesToIsolate', metavar = "[string]", help="The systematic ID of the gene to isolate", required=True, nargs='+')
    parser.add_argument('-o', '--OutputPrefix', metavar = "[string]", help="The prefix of the output [prefix].fa and [prefix].gtf file. If not specified, will use the value of the '-gene/--GeneToIsolate' argument as the prefix", required=False, default=None)
    parser.add_argument("--Quiet", help="quiet the output verbosity", action="store_true", default = False)
    parser.add_argument("--ReverseMinusStrandGenes", help="the output will use minus strand coordinates for minus strand genes; The sign of the strand for minus strand genes will be positive and the region Start and End coordinates will be determined by (Start = chromSize - OldEnd) and (End = chromSize - OldStart). The output fasta for that sequence will be the reverse compliment of the input fasta. Thus, if you opened the output fasta and annotation in IGV, the gene would be displayed left to right rather than the typical right to left for minus strand genes. Plus strand genes will be unaffected", action="store_true", default=False)
    args = parser.parse_args()
    print vars(args)
    main(**vars(args))