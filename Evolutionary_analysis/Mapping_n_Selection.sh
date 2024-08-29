########################################################################################
##############          Mapping of sample reads with pangenome          ###################
###########################################################################################


## 1.  Blacklisted_reads_removed Samples were mapping to reference genome using bwa-mem 

### On terminal window

module load bwa/0.7.12  ##load module
bwa index  Xp_Xeu.fna #also with Xp_Xeu.fna


## indexing of reference file used AL65 and AL22 non-redundant pangenome
##in queue

#!/bin/bash
source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load bwa/0.7.12

#input fastq files
for fq1 in /path/to/clean_reads/*.blk_filter.1.fq   ## made by the script mentioned in the https://github.com/Potnislab/AtDep_2021_metagenome/blob/main/Raw%20read%20processing/quality-trimming_script
do
echo "working with file $fq1"
base=$(basename $fq1 .blk.1.fq)
echo "base name is $base"
fq1=/path/to/1/${base}.blk.1.fq
fq2=/path/to/2/${base}.blk.2.fq
bwa mem /path/to/indexed_reference/Xp_Xeu.fna  $fq1 $fq2 > ${base}.sam    
done


## 2. Picard Sortsam  (https://gatk.broadinstitute.org/hc/en-us/articles/360046788792-SortSam-Picard-#--VALIDATION_STRINGENCY )

# downloaded picard
wget http://sourceforge.net/projects/picard/files/picard-tools/1.119/picard-tools-1.119.zip

# and gave the path to file

#!/bin/bash
source activate picard
for file in ./*.sam;    ## using the output files from previous step as input
    do tag=${file%.sam};
java -Xmx2g -jar /Path_to/picard-tools-1.119/SortSam.jar I="$file" O=$tag.sort.sam SORT_ORDER=coordinate
done


## 3. remove low quality alignments 
 
#!/bin/bash
module load samtools/1.11
for file in /path/to/*.sort.sam;   ##using output of sort sam as input in this step
    do tag=${file%.sort.sam}; 
samtools view -q 20 -f 0x0002 -F 0x0004 -F 0x0008 -Sb "$file" > $tag.q20.bam 
done 

## 4.Picard remove duplicates

#!/bin/bash
module load picard
for file in *.Xp_Xeu.q20.bam;
    do tag=${file%.Xp_Xeu.q20.bam};
picard MarkDuplicates -I "$file" -O $tag.Xp_Xeu.rmd.sort.bam -M $tag.dupstat.txt --REMOVE_DUPLICATES true --REMOVE_SEQUENCING_DUPLICATES true
done




########################################################################################
##############     Polymorphism  &   Genes_under_selection        ###################
###########################################################################################

# Input file requirement for METAPOP package
# to count no. of reads within input files
module load samtools
for each in *.bam; do data1=$(samtools view -c -F 4 "${each}");   data2=$(samtools view -c -f 4 "${each}"); printf '%s\t%s\t%s\n' "$each" "$data1" "$data2"; done > Reads_counts_Xp_Xeu.txt


#installation of Metapop with conda

conda install bioconda::metapop

#!/bin/bash
source activate metapop2
metapop --input_samples ./Bam --reference  reference \
--norm Reads_counts_Xp_Xeu.txt  --output Metapop_Output \
--plot_all --snp_scale both

# ./Bam is the folder where all the Bam files are located
# reference is the folder where pangenome as a refeerence file is located
# 