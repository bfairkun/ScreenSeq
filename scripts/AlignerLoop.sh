#!/usr/bin/env bash

#This script will align files on a screen-seq demultiplexed fastq files on loop and output spliced and unspliced counts for each file into ScreenSeqResults.tab To account for potential alternative splicing, the script first aligns all of the reads from all files, counts annotated and unannotated splice junctions with >10000 reads, and only reports counts for those junctions when aligning each sample individually.

RefGenomeFasta="MyPath/To/NewGenomeFiles.fa" #Change accordingly
RefGenomeGTF="MyPath/To/NewGenomeFiles.gtf" #Change accordingly
FastqDemultiplexedGlobPattern="FastqDemultiplex/*.fastq" #Change accordingly

samtools faidx $RefGenomeFasta #Index the genome. This is required for some bedtools commands
awk -v OFS='\t' '{print $1,$2}' $RefGenomeFasta.fai > $RefGenomeFasta.chrome.sizes #Make a file containing chromosome sizes. This is required for some bedtools commands
mkdir -p GenomeDir #Make a directory to hold the files that will be made by STAR genomeGenerate
STAR --runMode genomeGenerate --genomeDir ./GenomeDir/ --genomeFastaFiles $RefGenomeFasta --sjdbGTFfile $RefGenomeGTF --sjdbOverhang 100
cat ./FastqDemultiplexedOutAllGenes/* > AllScreenSeqReadsCombined.fastq #Combine all the files that were demultiplexed into a single fastq file. We will align this separately to identify which annotated and novel junctions are present at a high enough frequency that they are worth keeping track of when we do the alignments for each sample in a loop
mkdir -p CombinedAlignmentsNoMM #Make a new directory to hold the alignment files from the Combined fastq file
STAR --genomeDir ./GenomeDir --outFileNamePrefix ./CombinedAlignmentsNoMM/ --readFilesIn ./fastqfiles/AllDemultiplexableReads.fastq --genomeSAindexNbases 5 --outSAMtype BAM SortedByCoordinate --outSAMattributes All --clip5pNbases 23 --outReadsUnmapped Fastx --alignEndsType EndToEnd --clip3pAdapterSeq TCGAGCATATGATTTAACAAAGCGACCTGTCTCTTATACACATCTCCGAGCCCACGAGAC --runThreadN 8 --limitBAMsortRAM 20000000000 --alignIntronMax 290 --outFilterMismatchNmax 15
printf "FastqFile\tMappedReadCount\tUnmappedReadCount\t" > ScreenSeqResults.tab #Start by making a header for the ScreenSeqResults.tab file

#Grab annotated introns that have a Combined junction coverage >= 10000 and save them as regions in a bed file. These regions will be assessed for coverage later for each sample to quantify intron retention. Also append these junctions to the header of the ScreenSeqResults.tab file
awk -v OFS='\t' -F '\t' '$4==2 && $6==1 && $7>=10000 {print $1, $2-1, $3, ".", $7, "-"} $4==1 && $6==1 && $7>=10000 {print $1, $2-1, $3, ".", $7, "+"}' CombinedAlignmentsNoMM/SJ.out.tab | tee intronCovered.bed | awk -v OFS='\t' -v ORS='\t' -F '\t' '{print $1"_"$2+1"_"$3"_"$6"_Coverage", $1"_"$2+1"_"$3"_"$6"_Junctions"}' >> ScreenSeqResults.tab

#Grab unannotated introns that have a Combined coverage >= 10000 and save them as a single field to use later with the join command on each sample to grab the number of spliced reads supporting that junction. Also append these junctions to the header of the ScreenSeqResults.tab file
awk -v OFS='_' -F '\t' '$6==0 && $7>=10000 {print $1, $2, $3, $4}' CombinedAlignmentsNoMM/SJ.out.tab | tee CombinedAlignmentsNoMM/CommonUnnanotedEvents.txt | awk -v OFS='\t' -v ORS='\t' -F '_' '$4==2 {print "Unannotated_"$1"_"$2"_"$3"_-_Junctions"} $4==1 {print "Unannotated_"$1"_"$2"_"$3"_+_Junctions"}' >> ScreenSeqResults.tab
printf "OtherUnannotatedJunctionsCount\tGeneticVariants(Chr,Pos,Ref,Alt,Qual)" >> ScreenSeqResults.tab

#Align each sample individually on a loop. For each sample, extract the pertinent information and write it to a new row in the ScreenSeqResults.tab file. While the script is running, one could then continuously monitor progress through the loop using the tail command on the ScreenSeqResults.tab file
for myfilepath in $FastqDemultiplexedGlobPattern
do
    echo $myfilepath
    STAR --genomeDir ./GenomeDir --readFilesIn $myfilepath --genomeSAindexNbases 5 --outSAMtype BAM SortedByCoordinate --outSAMattributes All --clip5pNbases 23 --alignEndsType EndToEnd --clip3pAdapterSeq TCGAGCATATGATTTAACAAAGCGACCTGTCTCTTATACACATCTCCGAGCCCACGAGAC --genomeLoad LoadAndKeep --limitBAMsortRAM 2000000000 --alignIntronMax 290 --outReadsUnmapped Fastx --runThreadN 8 --outFilterMismatchNmax 15 > /dev/null
    samtools flagstat Aligned.sortedByCoord.out.bam | awk -v myvar="$myfilepath" -v ORS='\t' 'BEGIN {print "\n"'myvar'} NR==1 {print $1}' >> ScreenSeqResults.tab #Write the sample name and the number of mapped reads to output
    wc -l Unmapped.out.mate1 | awk -v myvar="$myfilepath" -v ORS='\t' '{print $1/4}' >> ScreenSeqResults.tab #Write number of unmapped reads to output
    awk -F '\t' '$4==2 && $6==1 && $7>=10 {print $1"_"$2"_"$3"_-", $7} $4==1 && $6==1 && $7>=10 {print $1"_"$2"_"$3"_+",$7}' SJ.out.tab > TempAnnotated_SJ.out.tab #Parse the read count of annotated junctions from STAR's SJ.out file to a temporary file
    
    #Use bedtools to count the coverage over the annotated intron regions that were previously stored as a bedfile. Match that intron coverage count with the corresponding junction count stored in the temporary file in the step above and write that information out to ScreenSeqResults.tab
    bedtools coverage -b ./Aligned.sortedByCoord.out.bam -a ./intronCovered.bed -sorted -split | awk -F '\t' '{print $1"_"$2+1"_"$3"_"$6,$7}' | join - TempAnnotated_SJ.out.tab -a 1 -e 0 -o 0,1.2,2.2 --nocheck-order | awk -v ORS='\t' -v OFS='\t' '{print $2, $3}' >> ScreenSeqResults.tab
    
    #Parse the read count of unannotated junctions from STAR's SJ.out file and filter it for only those that match unannotated introns previously stored as a bed file. Write the junction count to ScreenSeqResults.tab
    awk -v OFS='_' -F '\t' '{print $1, $2, $3, $4" "$7}' SJ.out.tab | join CombinedAlignmentsNoMM/CommonUnnanotedEvents.txt - --nocheck-order -a 1 -e 0 -o 2.2 | awk -v ORS='\t' {print} >> ScreenSeqResults.tab
    #Parse the read count of unannotated junctions from STAR's SJ.out file and filter it for only those that don't match unannotated introns previously stored as a bed file. Sum up that number (which represents count of rare unannoted junctions) and write it to ScreenSeqResults.tab
    awk -F '\t' -v OFS='\t' '$6==0 {print $1"_"$2"_"$3"_"$4, $7}' SJ.out.tab | awk -v ORS='\t' 'NR==FNR{a[$0];next} !($1 in a) {sum += $2} END {print sum}' CombinedAlignmentsNoMM/CommonUnnanotedEvents.txt - >> ScreenSeqResults.tab
    
    #Use samtools and bcftools to call SNPs, filtering out those with a quality score below 40. Write the SNPs to ScreenSeqResults.tab
    samtools mpileup -uf $RefGenomeFasta -E -d 100000 Aligned.sortedByCoord.out.bam 2> /dev/null | bcftools call -cv --ploidy 1 - | awk -v OFS=':' -v ORS=';' -F '\t' '/^[^#]/ && $6>=40 {print $1, $2, $4, $5, $6}' >> ScreenSeqResults.tab
done

rm -R GenomeDir