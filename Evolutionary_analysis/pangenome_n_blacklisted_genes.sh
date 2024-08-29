
#######################################################################################################################################################################################################################################
########################################################################################################################################################################################################################################
##################################################################     building pangenome using SuperPang    ####################################################################################################################################
#######################################################################################################################################################################################################################################
########################################################################################################################################################################################################################################

# For SuperPang, one of the input file required is checkm report of contamination and  completeness of genomes

#!/bin/bash
module load drep/3.2.2
checkm lineage_wf -t 14 --reduced_tree -x fna ./genome_files ./checkm

# fna is for the extension format
# ./genome_files folder where all the input/genome files located
#./ output folder name


# Building pangenome 

#!/bin/bash
SuperPang.py --fasta GCF_003993095.1_ASM399309v1_genomic.fna  GCF_006979915.1_ASM697991v1_genomic.fna   GCF_020879555.1_ASM2087955v1_genomic.fna GCF_003993105.1_ASM399310v1_genomic.fna  \
GCF_006979945.1_ASM697994v1_genomic.fna GCF_020879675.1_ASM2087967v1_genomic.fna GCF_003993115.1_ASM399311v1_genomic.fna  GCF_006979955.1_ASM697995v1_genomic.fna   \
GCF_020879695.1_ASM2087969v1_genomic.fna GCF_003993135.1_ASM399313v1_genomic.fna GCF_006979965.1_ASM697996v1_genomic.fna   GCF_020879715.1_ASM2087971v1_genomic.fna \
GCF_003993535.1_ASM399353v1_genomic.fna GCF_006979975.1_ASM697997v1_genomic.fna GCF_020879735.1_ASM2087973v1_genomic.fna GCF_003993575.1_ASM399357v1_genomic.fna \
GCF_006979985.1_ASM697998v1_genomic.fna GCF_020879955.1_ASM2087995v1_genomic.fna GCF_004102065.1_ASM410206v1_genomic.fna GCF_006980045.1_ASM698004v1_genomic.fna   \
GCF_020880155.1_ASM2088015v1_genomic.fna GCF_004102075.1_ASM410207v1_genomic.fna GCF_006980055.1_ASM698005v1_genomic.fna GCF_020880335.1_ASM2088033v1_genomic.fna \
GCF_004102085.1_ASM410208v1_genomic.fna GCF_006980075.1_ASM698007v1_genomic.fna GCF_020881595.1_ASM2088159v1_genomic.fna GCF_004102095.1_ASM410209v1_genomic.fna  \
GCF_006980085.1_ASM698008v1_genomic.fna GCF_021607285.1_ASM2160728v1_genomic.fna GCF_004102165.1_ASM410216v1_genomic.fna GCF_006980095.1_ASM698009v1_genomic.fna \
GCF_021607305.1_ASM2160730v1_genomic.fna GCF_000192045.2_ASM19204v3_genomic.fna GCF_004102175.1_ASM410217v1_genomic.fna GCF_006980115.1_ASM698011v1_genomic.fna \
GCF_021607325.1_ASM2160732v1_genomic.fna GCF_001009365.1_v1_genomic.fna GCF_004102205.1_ASM410220v1_genomic.fna GCF_006980165.1_ASM698016v1_genomic.fna GCF_021607345.1_ASM2160734v1_genomic.fna \
GCF_001009385.1_ASM100938v1_genomic.fna  GCF_004102215.1_ASM410221v1_genomic.fna  GCF_006980175.1_ASM698017v1_genomic.fna GCF_021607365.1_ASM2160736v1_genomic.fna GCF_001009395.1_v1_genomic.fna  \
 GCF_004102225.1_ASM410222v1_genomic.fna GCF_006980195.1_ASM698019v1_genomic.fna GCF_021607395.1_ASM2160739v1_genomic.fna GCF_001009405.1_v1_genomic.fna GCF_004102235.1_ASM410223v1_genomic.fna \
 GCF_006980215.1_ASM698021v1_genomic.fna GCF_021607415.1_ASM2160741v1_genomic.fna GCF_001009445.1_v1_genomic.fna GCF_004102275.1_ASM410227v1_genomic.fna GCF_006980225.1_ASM698022v1_genomic.fna \
 GCF_021607425.1_ASM2160742v1_genomic.fna GCF_001009465.1_v1_genomic.fna GCF_004102305.1_ASM410230v1_genomic.fna  GCF_006980245.1_ASM698024v1_genomic.fna   GCF_021607465.1_ASM2160746v1_genomic.fna \
 GCF_001009475.1_v1_genomic.fna GCF_004102315.1_ASM410231v1_genomic.fna GCF_006980275.1_ASM698027v1_genomic.fna GCF_021607485.1_ASM2160748v1_genomic.fna GCF_001009485.1_v1_genomic.fna \
 GCF_004102325.1_ASM410232v1_genomic.fna GCF_006980285.1_ASM698028v1_genomic.fna GCF_021607495.1_ASM2160749v1_genomic.fna GCF_001009545.1_ASM100954v1_genomic.fna  \
 GCF_004102335.1_ASM410233v1_genomic.fna  GCF_006980305.1_ASM698030v1_genomic.fna GCF_021607525.1_ASM2160752v1_genomic.fna GCF_001009665.1_v1_genomic.fna GCF_004102345.1_ASM410234v1_genomic.fna \
 GCF_006980315.1_ASM698031v1_genomic.fna GCF_021607545.1_ASM2160754v1_genomic.fna GCF_001009675.1_v1_genomic.fna GCF_004102405.1_ASM410240v1_genomic.fna GCF_006980365.1_ASM698036v1_genomic.fna \
 GCF_021607565.1_ASM2160756v1_genomic.fna GCF_001009685.1_v1_genomic.fna GCF_004102415.1_ASM410241v1_genomic.fna GCF_006980375.1_ASM698037v1_genomic.fna GCF_021607575.1_ASM2160757v1_genomic.fna \
 GCF_001009705.1_ASM100970v1_genomic.fna GCF_004102425.1_ASM410242v1_genomic.fna GCF_006980385.1_ASM698038v1_genomic.fna GCF_021607605.1_ASM2160760v1_genomic.fna GCF_001009745.1_v1_genomic.fna \
 GCF_004102435.1_ASM410243v1_genomic.fna GCF_006980395.1_ASM698039v1_genomic.fna GCF_021607625.1_ASM2160762v1_genomic.fna GCF_001009765.1_v1_genomic.fna GCF_004102455.1_ASM410245v1_genomic.fna \
 GCF_006980405.1_ASM698040v1_genomic.fna GCF_021607645.1_ASM2160764v1_genomic.fna GCF_001009795.1_v1_genomic.fna GCF_004102505.1_ASM410250v1_genomic.fna GCF_006980465.1_ASM698046v1_genomic.fna \
 GCF_021607655.1_ASM2160765v1_genomic.fna GCF_001009825.1_ASM100982v1_genomic.fna GCF_004102515.1_ASM410251v1_genomic.fna GCF_006980475.1_ASM698047v1_genomic.fna \
 GCF_021607685.1_ASM2160768v1_genomic.fna GCF_001009845.1_ASM100984v1_genomic.fna GCF_004102525.1_ASM410252v1_genomic.fna GCF_006980495.1_ASM698049v1_genomic.fna \
 GCF_021607705.1_ASM2160770v1_genomic.fna GCF_001009855.1_ASM100985v1_genomic.fna GCF_004102535.1_ASM410253v1_genomic.fna  GCF_006980505.1_ASM698050v1_genomic.fna \
 GCF_021607725.1_ASM2160772v1_genomic.fna GCF_001009865.1_ASM100986v1_genomic.fna GCF_006979525.1_ASM697952v1_genomic.fna GCF_006980525.1_ASM698052v1_genomic.fna \
 GCF_021607745.1_ASM2160774v1_genomic.fna GCF_001009885.1_ASM100988v1_genomic.fna  GCF_006979535.1_ASM697953v1_genomic.fna  GCF_006980565.1_ASM698056v1_genomic.fna \
 GCF_021607765.1_ASM2160776v1_genomic.fna GCF_001009925.1_ASM100992v1_genomic.fna GCF_006979545.1_ASM697954v1_genomic.fna GCF_006980575.1_ASM698057v1_genomic.fna \
 GCF_021607775.1_ASM2160777v1_genomic.fna GCF_001009935.1_ASM100993v1_genomic.fna GCF_006979555.1_ASM697955v1_genomic.fna GCF_006980585.1_ASM698058v1_genomic.fna \
 GCF_021607805.1_ASM2160780v1_genomic.fna GCF_001009945.1_ASM100994v1_genomic.fna GCF_006979565.1_ASM697956v1_genomic.fna GCF_006980615.1_ASM698061v1_genomic.fna \
 GCF_021607815.1_ASM2160781v1_genomic.fna GCF_001009955.1_ASM100995v1_genomic.fna GCF_006979625.1_ASM697962v1_genomic.fna GCF_006980625.1_ASM698062v1_genomic.fna \
 GCF_021607845.1_ASM2160784v1_genomic.fna GCF_001010005.1_ASM101000v1_genomic.fna GCF_006979635.1_ASM697963v1_genomic.fna GCF_006980655.1_ASM698065v1_genomic.fna \
 GCF_021607865.1_ASM2160786v1_genomic.fna GCF_001010015.1_ASM101001v1_genomic.fna GCF_006979645.1_ASM697964v1_genomic.fna GCF_006980665.1_ASM698066v1_genomic.fna \
 GCF_021607885.1_ASM2160788v1_genomic.fna GCF_001010025.1_ASM101002v1_genomic.fna GCF_006979655.1_ASM697965v1_genomic.fna GCF_007713955.1_ASM771395v1_genomic.fna \
 GCF_021607895.1_ASM2160789v1_genomic.fna GCF_001010035.1_ASM101003v1_genomic.fna GCF_006979675.1_ASM697967v1_genomic.fna GCF_007713965.1_ASM771396v1_genomic.fna \
 GCF_021607925.1_ASM2160792v1_genomic.fna GCF_001010085.1_ASM101008v1_genomic.fna GCF_006979715.1_ASM697971v1_genomic.fna GCF_007713985.1_ASM771398v1_genomic.fna \
 GCF_021607945.1_ASM2160794v1_genomic.fna GCF_001010105.1_ASM101010v1_genomic.fna GCF_006979725.1_ASM697972v1_genomic.fna GCF_007714045.1_ASM771404v1_genomic.fna \
 GCF_021607955.1_ASM2160795v1_genomic.fna GCF_001908855.1_ASM190885v1_genomic.fna GCF_006979735.1_ASM697973v1_genomic.fna GCF_007714065.1_ASM771406v1_genomic.fna \
 GCF_021607985.1_ASM2160798v1_genomic.fna GCF_001976075.1_ASM197607v1_genomic.fna GCF_006979785.1_ASM697978v1_genomic.fna GCF_007714075.1_ASM771407v1_genomic.fna \
 GCF_021608005.1_ASM2160800v1_genomic.fna GCF_003136155.1_ASM313615v1_genomic.fna GCF_006979805.1_ASM697980v1_genomic.fna GCF_007714105.1_ASM771410v1_genomic.fna \
 GCF_021608025.1_ASM2160802v1_genomic.fna GCF_003992975.1_ASM399297v1_genomic.fna GCF_006979815.1_ASM697981v1_genomic.fna GCF_007714115.1_ASM771411v1_genomic.fna \
 GCF_021608055.1_ASM2160805v1_genomic.fna GCF_003993015.1_ASM399301v1_genomic.fna GCF_006979835.1_ASM697983v1_genomic.fna GCF_009733625.1_ASM973362v1_genomic.fna LH3.fna \
 GCF_003993025.1_ASM399302v1_genomic.fna GCF_006979855.1_ASM697985v1_genomic.fna GCF_009733635.1_ASM973363v1_genomic.fna GCF_003993035.1_ASM399303v1_genomic.fna \
 GCF_006979875.1_ASM697987v1_genomic.fna GCF_013112235.1_ASM1311223v1_genomic.fna GCF_003993055.1_ASM399305v1_genomic.fna GCF_006979895.1_ASM697989v1_genomic.fna \
 GCF_020879295.1_ASM2087929v1_genomic.fna GCF_000488955.1_Xaf-CFBP3836-G1_genomic.fna GCF_003993445.1_ASM399344v1_genomic.fna GCF_000730305.1_G1_genomic.fna GCF_003993595.1_ASM399359v1_genomic.fna \
 GCF_000802325.1_XEU_66b_E1.0_genomic.fna GCF_003993605.1_ASM399360v1_genomic.fna GCF_000802345.1_XEU_83M_E1.0_genomic.fna GCF_003993615.1_ASM399361v1_genomic.fna GCF_001008805.1_v1_genomic.fna \
 GCF_003993655.1_ASM399365v1_genomic.fna GCF_001008815.1_v1_genomic.fna GCF_003993675.1_ASM399367v1_genomic.fna GCF_001008825.1_v1_genomic.fna GCF_003993685.1_ASM399368v1_genomic.fna \
 GCF_001008835.1_v1_genomic.fna GCF_003993725.1_ASM399372v1_genomic.fna GCF_001008885.1_v1_genomic.fna GCF_005059795.1_ASM505979v1_genomic.fna GCF_001008895.1_v1_genomic.fna \
 GCF_014198935.1_ASM1419893v1_genomic.fna GCF_001008905.1_v1_genomic.fna GCF_017724035.1_ASM1772403v1_genomic.fna GCF_001008915.1_v1_genomic.fna GCF_019192985.1_ASM1919298v1_genomic.fna \
 GCF_001008965.1_v1_genomic.fna GCF_019193005.1_ASM1919300v1_genomic.fna GCF_001008975.1_v1_genomic.fna GCF_020879135.1_ASM2087913v1_genomic.fna GCF_001008985.1_v1_genomic.fna \
 GCF_020879155.1_ASM2087915v1_genomic.fna GCF_001008995.1_v1_genomic.fna GCF_020879265.1_ASM2087926v1_genomic.fna GCF_001009045.1_v1_genomic.fna GCF_020879415.1_ASM2087941v1_genomic.fna \
 GCF_001009055.1_v1_genomic.fna GCF_020879485.1_ASM2087948v1_genomic.fna GCF_001009075.1_v1_genomic.fna GCF_020879635.1_ASM2087963v1_genomic.fna GCF_001009095.1_ASM100909v1_genomic.fna \
 GCF_020879655.1_ASM2087965v1_genomic.fna GCF_001009125.1_ASM100912v1_genomic.fna GCF_020879745.1_ASM2087974v1_genomic.fna GCF_001009135.1_ASM100913v1_genomic.fna \
 GCF_020879775.1_ASM2087977v1_genomic.fna GCF_001009165.1_ASM100916v1_genomic.fna GCF_020879795.1_ASM2087979v1_genomic.fna GCF_001009175.1_ASM100917v1_genomic.fna \
 GCF_020879815.1_ASM2087981v1_genomic.fna GCF_001009205.1_ASM100920v1_genomic.fna GCF_020879825.1_ASM2087982v1_genomic.fna GCF_001009215.1_ASM100921v1_genomic.fna \
 GCF_020879845.1_ASM2087984v1_genomic.fna GCF_001009245.1_ASM100924v1_genomic.fna GCF_020879875.1_ASM2087987v1_genomic.fna GCF_001009255.1_ASM100925v1_genomic.fna \
 GCF_020880015.1_ASM2088001v1_genomic.fna GCF_001010095.1_ASM101009v1_genomic.fna GCF_020880035.1_ASM2088003v1_genomic.fna GCF_001401555.1_Xeu_LMG27970_genomic.fna \
 GCF_020880045.1_ASM2088004v1_genomic.fna GCF_001401625.1_Xaa_LMG495_genomic.fna GCF_020880115.1_ASM2088011v1_genomic.fna GCF_001401675.2_ASM140167v2_genomic.fna \
 GCF_020880135.1_ASM2088013v1_genomic.fna GCF_001691315.1_XEU_LMG667_E1.0_genomic.fna GCF_020880165.1_ASM2088016v1_genomic.fna GCF_001691325.1_XEU_LMG905_E1.0_genomic.fna \
 GCF_020880235.1_ASM2088023v1_genomic.fna GCF_001691345.1_XEU_LMG909_E1.0_genomic.fna GCF_020880275.1_ASM2088027v1_genomic.fna GCF_001691375.1_XEU_LMG933_E1.0_genomic.fna \
 GCF_020880315.1_ASM2088031v1_genomic.fna GCF_001691385.1_XEU_LMG918_E1.0_genomic.fna GCF_020880345.1_ASM2088034v1_genomic.fna GCF_001908795.1_ASM190879v1_genomic.fna \
 GCF_020880375.1_ASM2088037v1_genomic.fna GCF_002939715.1_XaclCFBP3371_genomic.fna GCF_020880395.1_ASM2088039v1_genomic.fna GCF_003136175.1_ASM313617v1_genomic.fna \
 GCF_020880415.1_ASM2088041v1_genomic.fna GCF_003992785.1_ASM399278v1_genomic.fna GCF_020880435.1_ASM2088043v1_genomic.fna GCF_003992805.1_ASM399280v1_genomic.fna \
 GCF_020880575.1_ASM2088057v1_genomic.fna GCF_003993175.1_ASM399317v1_genomic.fna GCF_020880635.1_ASM2088063v1_genomic.fna GCF_003993185.1_ASM399318v1_genomic.fna \
 GCF_020880835.1_ASM2088083v1_genomic.fna GCF_003993195.1_ASM399319v1_genomic.fna GCF_020880855.1_ASM2088085v1_genomic.fna GCF_003993225.1_ASM399322v1_genomic.fna \
 GCF_020880895.1_ASM2088089v1_genomic.fna GCF_003993255.1_ASM399325v1_genomic.fna GCF_003993265.1_ASM399326v1_genomic.fna GCF_003993275.1_ASM399327v1_genomic.fna \
 GCF_003993315.1_ASM399331v1_genomic.fna GCF_003993335.1_ASM399333v1_genomic.fna GCF_003993345.1_ASM399334v1_genomic.fna --checkm checkm_result_xp_xeu \
 --output-dir Xp_Xeu_Pangenome --force-overwrite --verbose-mOTUpan -t 10


# all *.fna are the genome files as input
# Xp_Xeu_Pangenome is the output directory 


#######################################################################################################################################################################################################################################
########################################################################################################################################################################################################################################
##################################################################    blacklisting of genes and removal of reads    ######################################################################################################################
#######################################################################################################################################################################################################################################
########################################################################################################################################################################################################################################


# finding blacklisiting genes

#1. step is to find the abundant species present in the phyllosphere

#installation of MIDAS software with conda 

conda install bioconda::midas


# downloading the midas default database of all species


git clone https://github.com/snayfach/MIDAS
cd MIDAS
wget http://lighthouse.ucsf.edu/MIDAS/midas_db_v1.2.tar.gz 
tar -zxvf midas_db_v1.2.tar.gz


# now run the command for species identification within the samples

#!/bin/bash
source activate Midas
for fq1 in ./Input_clean_reads/*_1.fq;
do
echo "working with file $fq1"

base=$(basename $fq1 _1.fq)

echo "base name is $base"

fq1=./Input_clean_reads/${base}_1.fq
fq2=./Input_clean_reads/${base}_2.fq

run_midas.py species Output_${base} -1 $fq1 -2 $fq2 -d midas_db_v1.2

done


#midas_db_v1.2 is the default database


#2. combing all the output files from different samples using merge command in the midas. Here is one example for single season samples 

conda  activate Midas
merge_midas.py species ./Merge_FFAR_2022_species -i \
ffar2022_AALE22,ffar2022_CHGA22,ffar2022_EVS22,ffar2022_NGAM22,ffar2022_SALM22,ffar2022_SSCM22,ffar2022_AALM22,ffar2022_DSCE22,ffar2022_LALE22,ffar2022_PUSC22,ffar2022_STSC22,ffar2022_BALE22,\
ffar2022_DSCM22,ffar2022_LALM22,ffar2022_RPGA22,ffar2022_WSC22,ffar2022_BALM22,ffar2022_ESC22,ffar2022_NGAE22,ffar2022_SALE22 -t list -d midas_db_v1.2

# ./Merge_FFAR_2022_species  is the output file

#below are 
#then, we took species prevalence file from all the three year samples (from compiled MIDAS species results of all the three years), and took out the genomes showing greater than zero mean prevalence.

# those are (tags according to MIDAS default database) :
#Pseudomonas_oleovorans_57108
#Pseudomonas_sp_59673
#Xanthomonas_perforans_55843
#Xanthomonas_axonopodis_56719
#Xanthomonas_cassavae_59587
#Xanthomonas_axonopodis_61257
#Xanthomonas_axonopodis_57683
#Pantoea_vagans_57743
#Pseudomonas_fulva_58092
#Enterobacter_cloacae_60366
#Sphingomonas_parapaucimobilis_59581
#Pantoea_sp_61911
#Pseudomonas_viridiflava_57558
#Xanthomonas_arboricola_57436
#Pantoea_sp_60701
#Pantoea_ananatis_58251
#Pantoea_agglomerans_54643
#Pantoea_agglomerans_56951
#Enterobacter_cancerogenus_61658
#Enterobacter_cloacae_58148
#Enterobacteriaceae_bacterium_60920
#Escherichia_vulneris_58627
#Pseudomonas_cichorii_60533
#Sphingomonas_sp_60678
#Enterobacter_cloacae_55011
#Enterobacter_cloacae_60571
#Pseudomonas_syringae_54805
#Pseudomonas_syringae_59460
#Acinetobacter_baumannii_57014
#Methylobacterium_populi_61518

##we took all the listed reference genomes into single directroy from midas database and concatenated all the listed genes files into single fasta file:

zcat ./Pan_genome/*/centroids.ffn.gz  >> all_species_pan_genomes_genes.fasta  



# 3. In parallel, we formatted our pangenome reference file according to MIDAS. 

#Build custom database

#input file for custom database = 
#According to github page, we need following files: indir: Path to directory of genomes. Each subdirectory should be named with a genome identifier. 
#Each subdirectory should contain the following files; 
#1. <genome_id>.fna: Genomic DNA sequence in FASTA format; 
# 2.<genome_id>.faa: Protein sequences in FASTA format; 
# 3. <genome_id>.ffn: Gene sequences in FASTA format; 
# 4. <genome_id>.genes: Tab delimited file with genomic coordinates of genes. 
#The file should be tab-delimited file with a header and the following fields:  gene_id (CHAR); scaffold_id (CHAR); start (INT); end (INT); strand (+ or -);    gene_type (CDS or RNA); mapfile: Path to mapping file that specifies which genomes belonging to the same species.
#5. The file should be tab-delimited file with a header and 3 fields:   genome_id (CHAR): corresponds to subdirectory within INDIR; species_id (CHAR): : species identifier for genome_id; rep_genome (0 or 1): indicator if genome_id should be used for SNP calling


#The first three files was built using prokka

prokka --prefix Xp_Xeu --locustag Xp_Xeu --increment 10 --outdir Xp_Xeu_Prokka --force --addgenes --genus Xanthomonas --gcode 11 assembly.fasta

#And the fourth one was made using gff  file from prokka using the following software.

#It was installed using conda. https://github.com/shenwei356/csvtk#installation 


# we got 1,2,3, file from above list
# to get .genes file we will convert the .gff file using csvtk package to cut the specific columns  from file 
#To select fields/columns:
conda activate csvtk
csvtk cut -f 1,3-5,7,9 --ignore-illegal-row -t Xp_Xeu.gff > Xp_Xeu.genes

#To rename fields/columns in genes file:
csvtk rename -f 1-6 -t -n scaffold_id,gene_type,start,end,strand,gene_id Xp_Xeu.genes > Xp_Xeu.genes


#genome mapfile looked like this
#genome_id   species_id  rep_genome
#Xp_LH3  Xanthomonas_perforans   1

#custom database =
#hmmer pacakge is required

source activate Midas
source /apps/profiles/modules_asax.sh.dyn
module load vsearch/2.22.1
module load hmmer

build_midas_db.py genomes  genomes.mapfile  Xp_Xeu_custom_db  --compress



# 4. finding genes: Now blast the both species database and Xp_Xeu database to find if there are any genes are similar
##build database for blast
source /apps/profiles/modules_dmc.sh.dyn
module load blast+
makeblastdb -in all_species_pan_genomes_genes.fasta -out all_species_pan_genomes_genes_db -dbtype nucl

# all_species_pan_genomes_genes.fasta ( from 193 line) is the input file consisting concatenated gene database from species present in the phyllosphere


#blasted the database with the Xanthomonas perforans database from Midas

blastn -db all_species_pan_genomes_genes_db -query Xp_Xeu_genes.ffn  -outfmt 6 -out blast_output_Xp_Xeu.txt 

#separted the genes with more than 97% identity

#Then, we took all the gene with more than 97% identity from the blast result and extracted the list of those genes out of our databases, we called these genes as blacklisted genes
git clone https://github.com/lh3/seqtk.git; 
cd seqtk; make 
seqtk subseq  ../Xp_Xeu_genes.ffn  ../genes_Xp_Xeu > ./Xp_Xeu-blk.fasta

#./genes_Xp_Xeu contains gene names

#4. removal of reads mapping with blacklisted genes


#!/bin/bash
source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load bbmap

for fq1 in /path_to_input_clean_reads/*.paired.1.fq
    do
    echo "working with file $fq1"

    base=$(basename $fq1 .paired.1.fq)
    echo "base name is $base"

    fq1=/path_to_input_sample_files/${base}.paired.1.fq
    fq2=/path_to_input_sample_files/${base}.paired.2.fq


bbduk.sh -Xmx30g in1=$fq1 in2=$fq2 out1=${base}_blk_filter.1.fq out2=${base}_blk_filter.2.fq outm1=${base}_matched.1.fq outm2=${base}_matched.2.fq \
ref=./Xp_Xeu-blk.fasta k=31 hdist=1 stats=${base}.stats.txt

done




