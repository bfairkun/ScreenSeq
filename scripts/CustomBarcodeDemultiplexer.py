# -*- coding: utf-8 -*-
#!/usr/bin/env python

import os
import HTSeq
import time
import re
import sys
import itertools
import cStringIO
import argparse
import distutils.dir_util

def extend_ambiguous_dna(seq, AmbiguousDict={"N":["A", "G","C","T"],"V":["A", "G","C"],"H":["A","C","T"],"D":["A", "G","T"],"B":["G","C","T"],"M":["A","C"],"K":["G","T"],"W":["A","T"],"S":["G","C"],"Y":["C","T"],"R":["A","G"]}):
    """
    Given a DNA sequence that may contain ambiguous bases (ex N for any base, Y for pyrimidines, etc.) as a string, this function returns as list of non-ambiguous sequences with only A, G, C, or T
    """
    groups = itertools.groupby(seq, lambda char:char not in AmbiguousDict)
    splits = []
    for b,group in groups:
        if b:
            splits.extend([[g] for g in group])
        else:
            for nuc in group:
                splits.append(AmbiguousDict[nuc])
    return [''.join(p) for p in itertools.product(*splits)]

def flush_stringio_buffers_to_disk(StringIOBufferDict, FilesToAppendInsteadOfWrite=None):
    """
    Given a dictionary that maps filepaths (the keys) to stringio objects (the values), this function will write out the stringio objects to the files, making directories as necessary if they are included in the filepath
    """
    for filename, stringio_buffer in StringIOBufferDict.items():
        distutils.dir_util.mkpath(os.path.dirname(filename))
        if filename in FilesToAppendInsteadOfWrite:
            fh = open(filename, 'a')
        else:
            fh = open(filename, 'w')
        stringio_buffer.seek(0)
        fh.write(stringio_buffer.read())
        fh.close()

def restart_line():
    sys.stdout.write('\r')
    sys.stdout.flush()

def main():
    parser = argparse.ArgumentParser(usage = '''
    This script is used to demultiplex a fastq file into daughter fastq files
    based on a series of barcodes as supplied by the user in a tab-delimited-text
    file.
    
    It is assumed that each fastq read will have a 500 series barcode (first index read) and a 700 series barcode (second index read) in the header line of the read. Also it is assumed that each insert read will have a plate-specific barcode and a primer-specific sequence as the first bases sequenced in the insert read. (See the example below:)
    
    Given a user supplied barcode key text file that looks like this:
    
    Filepath(Absolute or Relative)	500 series barcode	700 series barcode	plate barcode	primer barcode
    my_folder_for_GCTGCCCATAAA_primer/my_sample1.fastq	TAAGGCGA	ATAGAGAG	GTCAG	GCTGCCCATAAA
    my_folder_for_GCTGCCCATAAA_primer/my_sample2.fastq	CGTACTAG	ATAGAGAG	GTCAG	GCTGCCCATAAA
    my_folder_for_GGAAGGAAATTG_primer/my_sample3.fastq	TCGACGTC	TCGCATAA	NCCGTGAGA	GGAAGGAAATTG
    
     ...and a fastq file that looks like this:
    
    @NB500947:241:HTVMFBGXY:4:13608:10075:9804 1:N:0:TAAGGCGA+ATAGAGAG
    GTCAGGCTGCCCATAAACTCTTTGCCCTCAAAGTCATTTACGATATCACGAGCATCGCGAGCGTCTTCAACCTCAA
    +
    AAAAAEEEEEEEEEEEEEEEEEEEEEEE/EEEEEEEAEEEEEEEEEEEEEEEEEEEEEAEEEEEEAEAEEEEAEEE
    @NB500947:241:HTVMFBGXY:1:11305:11868:13191 1:N:0:CGTACTAG+ATAGAGAG
    GTCAGGCTGCCCATAAACTCTTTGCCCTGAAAGTCATTTACGATATCACGAGCATCGCGAGCGTCTTCAACCTCAA
    +
    AA/AAEEEEEEEEEEEEEEEEEEEEEEAEEEEEEEEEEEEAEAAEEEEE</EEEAEEA<AEAEEEEEEEEEEAAEE
    @NB500947:241:HTVMFBGXY:1:12109:5082:17813 1:N:0:CGTACTAG+ATAGAGAG
    GTCAGGCTGCCCATAAACTCTTTGCCCTGAAAGTCATTTACGATATCACGAGCATCGCGAGCGTCTTCAACCTCAA
    +
    A/A/AEEEEEAEEEEEEEEEEEEEEEAEEEEEAEEA/E6AEAEEEEEEEE/E<AEEEEA/EE/E//EA//EE/AEA
    @NB500947:241:HTVMFBGXY:1:11112:15794:18211 1:N:0:TCGACGTC+TCGCATAA
    CCCGTGAGAGGAAGGAAATTGACCATATTAAGTTCGTGAAACTTTGGGTCTGTCTCTTATACACATCTCCAAGCCC
    +
    AAAAAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE/EEEEE
    @NB500947:241:HTVMFBGXY:1:21210:6296:17810 1:N:0:TCGACGTC+TCGCATAA
    TCCGTGAGAGGAAGGAAATTGACCATATTAAGTTCGTGAAACTTTGGGTCTGTCTCTTATACACATCTCCGAGCCC
    +
    AAAAAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE/EEEEEEEEEEEEEEEEEEEEEEEAEEEAEE/EEEE
    
    Then the first read will be written to my_folder_for_GCTGCCCATAAA_primer/my_sample1.fastq, the next two reads to my_folder_for_GCTGCCCATAAA_primer/my_sample2.fastq, and the last two reads to my_folder_for_GGAAGGAAATTG_primer/my_sample3.fastq. Directories will be made as necessary if they are in the filepath in the barcodekey file. Also, notice that plate barcodes in the barcokekey file may have "N" for any nucleotide in the barcode.
    
     
    DEPENDENCIES:
    This python script requires the python module 'HTSeq' to be installed.
    '''
    )
    parser.add_argument('-I','--InputFile', metavar = "[filepath]", help="REQUIRED. The input fastq file. Use 'stdin' or '-' to read from stdin", required=True)
    parser.add_argument('-Key', '--BarcodeKey', metavar = "[filepath]", help="REQUIRED. File that describes barcodes and filepaths to separate into", required=True)
    parser.add_argument('-U', '--UnmatchedReads', metavar = "[filepath]", help="OPTIONAL. Filepath for fastq reads which do not match any barcode combinations", required=False, default=False)
    parser.add_argument("--quiet", help="OPTIONAL. quiet the output verbosity", action="store_true")
    parser.add_argument("--MaxReadsInMemory", metavar = "[integer]", help="OPTIONAL. Max number of reads in file buffers memory at once before flushing to files on disk. Default = 50000", type=int, default=50000)
    parser.add_argument("--MaxFileBuffersInMemory", metavar = "[integer]", help="OPTIONAL. Max number of file buffers in memory at once, Default = 500", type=int, default=500)
    parser.add_argument("--CreateEmptyFastqFiles", help="OPTIONAL. Create a fastq file for every filepath in the BarcodeKey file, regardless of if it ends up being an empty file without any barcode matched reads", action="store_true")
    args = parser.parse_args()
    
    FilesOpened = set()
    BarcodeFile = open(args.BarcodeKey, 'rU')
    BarcodeFile.readline()
    BarcodeDict={} #will hold a tuple of the four barcodes as keys to access the corresponding destination filepath
    for line in BarcodeFile:
        filepath,barcode_500,barcode_700,barcode_plate,barcode_primer = line.strip('\n').split('\t')
        if args.CreateEmptyFastqFiles:
            distutils.dir_util.mkpath(os.path.dirname(filepath))
            open(filepath, 'w').close()
            FilesOpened.add(filepath)
        for UnambiguousPlateBarcode in extend_ambiguous_dna(barcode_plate):
            BarcodeDict[(barcode_500,barcode_700,UnambiguousPlateBarcode,barcode_primer)] = filepath
    
    
    BarcodeFile.close()
    barcode_primer_StartIndex = min([len(barcode[2]) for barcode in BarcodeDict.keys()]) #The first possible position a primer_barcode may start in the read
    barcode_primer_StopIndex = max([len(barcode[2]) for barcode in BarcodeDict.keys()]) + max([len(barcode[3]) for barcode in BarcodeDict.keys()]) #The last possible position a primer_barcode may end in the read
    
    PrimerBarcodes = re.compile('|'.join(set([barcode[3] for barcode in BarcodeDict.keys()])))
    tic = time.clock()
    
    MaxOpenFiles = 500
    FileBufferDict = {}
    ReadsInMemory = 0
    BarcodeMatchedReadCount = 0
    BarcodeUnmatchedReadCount = 0
    if args.InputFile == 'stdin' or args.InputFile == '-':
        FastqIn = sys.stdin
    else:
        FastqIn = args.InputFile
    for i,read in enumerate(HTSeq.FastqReader(FastqIn)):
        if i>=0:
            OutFile = None
            barcode_500,barcode_700 = read.name.split('N:0:')[-1].split('+')
            RegexSearchObject = PrimerBarcodes.search(str(read[barcode_primer_StartIndex:barcode_primer_StopIndex])) #re.search returns None if there was no match
            if RegexSearchObject:
                barcode_plate = str(read)[0:RegexSearchObject.span()[0]+barcode_primer_StartIndex]
                barcode_primer = RegexSearchObject.group()
                try:
                    OutFile = BarcodeDict[(barcode_500,barcode_700,barcode_plate,barcode_primer)]
                    if OutFile not in FileBufferDict:
                        FileBufferDict[OutFile] = cStringIO.StringIO()
                    read.write_to_fastq_file(FileBufferDict[OutFile])
                    ReadsInMemory += 1
                    BarcodeMatchedReadCount += 1
                except KeyError:
                    pass
            if args.UnmatchedReads and OutFile is None:
                if args.UnmatchedReads not in FileBufferDict:
                    FileBufferDict[args.UnmatchedReads] = cStringIO.StringIO()
                read.write_to_fastq_file(FileBufferDict[args.UnmatchedReads])
                ReadsInMemory += 1
                BarcodeUnmatchedReadCount += 1 
            if len(FileBufferDict) >= args.MaxFileBuffersInMemory or ReadsInMemory >= args.MaxReadsInMemory:
                flush_stringio_buffers_to_disk(FileBufferDict,FilesToAppendInsteadOfWrite=FilesOpened)
                FilesOpened.update(FileBufferDict.keys())
                FileBufferDict = {}
                ReadsInMemory = 0
            if not args.quiet:
                if i % 100000 == 0:
                    print(time.strftime('%X') + ' | ' + format((i)/float(1000000), '0.1f') + 'million reads processed')
    flush_stringio_buffers_to_disk(FileBufferDict,FilesToAppendInsteadOfWrite=FilesOpened)
    toc = time.clock()
    if not args.quiet:
        print('DONE')
        print('Elapsed time: ' + str(toc-tic))
        print(str(BarcodeMatchedReadCount) + '/' + str(i) + ' reads had identifiable barcodes, which demultiplexed into ' + str(len(FilesOpened)) + ' files')
        if args.UnmatchedReads:
            print(str(BarcodeUnmatchedReadCount) + '/' + str(i) + ' reads had unidentifiable barcodes, were written to ' + args.UnmatchedReads)
main()
