
################################################################################################################################################
################################################################################################################################################
################################################         SEASONAL VARIANTS          ############################################################
################################################################################################################################################
################################################################################################################################################

#1. used MIDAS package for calling the SNV sites

#!/bin/bash
source activate Midas
for fq1 in /path_to_input_reads/*.blk.1.fq    # here we used samples without blacklisted genes
do
echo "working with file $fq1"
base=$(basename $fq1 .blk.1.fq)
echo "base name is $base"
fq1=./${base}.blk.1.fq
fq2=./${base}.blk.2.fq

run_midas.py snps MIDAS_SNP/${base}_Snp_bk --species_id Xanthomonas_perforans \
-1 $fq1   -2 $fq2 -d Xp_Xeu_custom_db  --remove_temp
done


# Xp_Xeu_custom_db is the database formated using MIDAS before identifying the genes similar to phylosphere

# now we merged all the output files based on the single sampling point (one sampling point = Single population) here is the example for one.

# we removed the samples where we did not have any Xanthomonas presence in the population.

source activate Midas
merge_midas.py snps ./Snps_SUS_A_End --species_id Xanthomonas_perforans -i A_ALE_Snp_bk,B_ALE_Snp_bk,D3_SCM_Snp_bk,D_SCE_Snp_bk,S_SCE_Snp_bk,L_ALE_Snp_bk,S_ALE_Snp_bk,\
L_NCE_Snp_bk,U_SCE_Snp_bk -t list -d Xp_Xeu_custom_db --all_sites --all_samples  --snp_type any

######################################################################################################################################################################################################
######################################################################################################################################################################################################
      # for finding the average frequency we used the following command used directly on the output of merged files from single sampling point
######################################################################################################################################################################################################
######################################################################################################################################################################################################


awk -F',' -v OFS=',' 'NR==1{ 
    print $0, "Average depth"; 
    next 
} 
{ 
    s=0; 
    numFields=0; 
    for(i=2; i<=NF;i++){ 
        if(length($i)){ 
            s+=$i; 
            numFields++ 
        } 
    } 
    print $0, (numFields ? s/numFields : 0) 
}' snps_freq_20_mid.csv > Av_freq_20_mid.csv

# we did same for averaging the depth 

# then we combined the two files together and subset the files based on the read depth (removed the sites where we had less than 10 average read depth)

# further to explore the changes in the frequency of alleles during seasons, we used R to subset and ploting the figure.

#######################################################################################################################################################################################################
#######################################################################################################################################################################################################


## 2.  Now finding the variants of interest finding the alleles with freq less than 0.2 and became greater than 0.8 next season with changes in allele from minor to major

#we have done for all possible combinations such as :
# from mid 2020: to end 2020, to fall 2020, to mid 2021, to end 2021, to mid 2022, to end 2022, to fall 2022
# from end 2020: to fall 2020, to mid 2021, to end 2021, to mid 2022, to end 2022, to fall 2022
# from fall 2020: to mid 2021, to end 2021, to mid 2022, to end 2022, to fall 2022
# from mid 2021: to end 2021, to mid 2022, to end 2022, to fall 2022
# from end 2021: to mid 2022, to end 2022, to fall 2022
# from mid 2022: to end 2022, to fall 2022
# from end 2022: to fall 2022

######################################################################################################################################################################################################

#here is the example for one of the combination:

# we got four different output files from last step of merging the samples with single sampling point

# snp_depth.txt; snp_info.txt, snp_freq.txt, readme.txt


#I. combining the  snp_info files together based on the combination, from mid_2020 to end 2020 in this case.

awk -F, -vOFS=, '(NR==FNR){a[$1]=$0; next}  { if(a[$1]){print $0,a[$1]} else{print $0,"no match"}}' snps_info_20_end.csv snps_info_20_mid.csv > Snp_info_20_Mid_End.csv

######################################################################################################################################################################################################

#II. now find the major minor changes during the season if major allele changed or not

# removed the column using "cut" function and kept which are required. here are the column names of the final file

# 1 site_id 
# 2 ref_id  
# 3 ref_pos 
# 4 ref_allele  
# 5 major_allele  
# 6 minor_allele  
# 7 count_samples 
# 8 major_allele  
# 9 minor_allele  
# 10 count_samples

# to find if there are same major from mid to end season: creating new column, by checking if each row of column 5 is equal to column 7, then type true
awk -F, -vOFS=,  'NR==1 { $0 = $0 OFS  "Same_major"; print; next } {if ($5==$8) {print $0 OFS "True"} else {print $0 OFS "False"}}' OFS=, Snp_info_20_Mid_End_1.csv > tmp

# to find if there are same minor allele from mid to end season: creating new column
awk  -F, -vOFS=,  'NR==1 { $0 = $0 OFS  "Same_minor"; print; next } {if ($6==$9) {print $0 OFS "True"} else {print $0 OFS "False"}}' OFS=, tmp > tmp1

# to find allele is shifting from minor to major from mid to end season: creating new column
awk  -F, -vOFS=, 'NR==1 { $0 = $0 OFS  "Min_major"; print; next } {if ($6==$8) {print $0 OFS "True"} else {print $0 OFS "False"}}' OFS=, tmp1 > tmp2

# to find allele is shifting from major to minor from mid to end season: creating new column
awk -F, -vOFS=,  'NR==1 { $0 = $0 OFS  "maj_minor"; print; next } {if ($5==$9) {print $0 OFS "True"} else {print $0 OFS "False"}}' OFS=, tmp2 > Snp_info_20_Mid_End_Final.csv

######################################################################################################################################################################################################


#III.  subsetting the sites : keeping if the minor allele changed and then from those if the minor from previous season became major

# now the last output file has these column names

# 1 site_id 
# 2 ref_id  
# 3 ref_pos 
# 4 ref_allele  
# 5 major_allele  
# 6 minor_allele  
# 7 count_samples 
# 8 major_allele  
# 9 minor_allele  
# 10 count_samples
# 11 Same_major 
# 12 Same_minor
# 13  Min_major 
# 14 maj_minor

#Snp_info_20_Mid_End_Final.csv

awk -F , '$12=="False" { print }' ../Snp_info_20_Mid_End_Final.csv > Snp_info_20_Mid_End_Final_minor.csv
awk -F , '$13=="True" { print }' Snp_info_20_Mid_End_Final_minor.csv > Changed_Snp_info_20_Mid_End.csv

#######################################################################################################################################################################################################

# IV. removing alleles based  on snp freq and read depth

# for that we need depth and freq from both seasons, so combined the files
awk -F, -vOFS=, '(NR==FNR){a[$1]=$0; next}  { if(a[$1]){print $0,a[$1]} else{print $0,"no match"}}' snps_depth_20_mid.csv snps_freq_20_mid.csv  > FD_20_mid.csv 
awk -F, -vOFS=, '(NR==FNR){a[$1]=$0; next}  { if(a[$1]){print $0,a[$1]} else{print $0,"no match"}}' snps_depth_20_end.csv snps_freq_20_end.csv  > FD_20_end.csv 


# NOW combine these above files to the changes file ( from line 115 )
awk -F, -vOFS=, '(NR==FNR){a[$1]=$0; next}  { if(a[$1]){print $0,a[$1]} else{print $0,"no match"}}' FD_20_mid.csv Changed_Snp_info_20_Mid_End.csv  > tmp
awk -F, -vOFS=, '(NR==FNR){a[$1]=$0; next}  { if(a[$1]){print $0,a[$1]} else{print $0,"no match"}}' FD_20_end.csv tmp5  > FD_Changed_Snp_info_20_Mid_End.csv 


# Next we used R to remove the sites with depth less than 10 and kept if they are present in more than 10 read depth in at least 50% of the samples.
# we looked for the frequency changes if the allele was less than 0.2 and became > 0.8 in the next season

# Then, similarly, we find if these allele which became minor to major are staying as major during the next seasons (next sampling points).

