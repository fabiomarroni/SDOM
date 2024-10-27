
#!/bin/bash

##############################################
#
# Variables you will need to set based on your file system and/or local properties
# Set variables to the desired values and then remove the comment from beginning of line
#
# FUNCDIR=/path/to/folder/with/R/functions
# BASEDIR=/base/project/folder
# READDIR=${BASEDIR}/00_raw_reads
# FILTDIR=${BASEDIR}/00a_filt_reads
# CONTAMDIR=${BASEDIR}/00b_contam_reads
# BBMAPDIR=/if/not/in/path/directory/to/bbmap/executable
# OMREF=path/to/Oncorhynchus_mykiss/reference/genome
# KRES=${BASEDIR}/02a_filt_kraken2_16s_RefSeq
# KRAKBIN=/if/not/in/path/directory/to/kraken2/executable
# KDB=/path/to/kraken2/databases/
# BRACK_BIN=/if/not/in/path/directory/to/Bracken2/executable
#
#
# READPATTERN=_R1_001.fastq.gz
#
# End of variable specification
#
###############################################


############################################
#  ___ ___     _      _   
# | _ ) _ ) __| |_  _| |__
# | _ \ _ \/ _` | || | / /
# |___/___/\__,_|\_,_|_\_\
#
############################################

#Remove trout contamination with bbduk


READPATTERN=_R1_001.fastq.gz
#End of parameters that you likely need to change
mem=32g
CPU=8
#Run kraken on the selected subreads
# PRE Execution
#s1 - classify
for aaa in $READDIR/*${READPATTERN}
do
read1=$(basename $aaa)
read2=${read1/_R1_001.fastq.gz/_R2_001.fastq.gz}
echo $read1
echo $read2
${BBMAPDIR}/bbduk.sh in=$READDIR/${read1} in2=$READDIR/${read2} out=${FILTDIR}/${read1} out2=${FILTDIR}/${read2} \
    -Xmx${mem} thread=${CPU} speed=0 outm=${CONTAMDIR}/${read1} outm2=${CONTAMDIR}/${read2} ref=${OMREF}
done




############################################
#
#  _              _              ___  
# | |            | |            |__ \ 
# | | ___ __ __ _| | _____ _ __    ) |
# | |/ / '__/ _` | |/ / _ \ '_ \  / / 
# |   <| | | (_| |   <  __/ | | |/ /_ 
# |_|\_\_|  \__,_|_|\_\___|_| |_|____|
#                                     
############################################
                                     

######################
#
# kraken2 on refseq
#
######################



#Create a kraken database from the 16s refseq database
mkdir -p $KRES

READPATTERN=_R1_001.fastq.gz
#End of parameters that you likely need to change

#Run kraken 
for aaa in $FILTDIR/*${READPATTERN}
do
read1=$(basename $aaa)
read2=${read1/_R1_001.fastq.gz/_R2_001.fastq.gz}
pref=${read1/_R1_001.fastq.gz/}
echo $read1
echo $read2
mkdir -p ${KRES}/logs
if [ -s ${KRES}/${pref}.kraken.report.txt ]
	then
	echo "Kraken already run, i will skip it"	
else
cd ${KRES} 
export TMPDIR=${KRES}
${KRAKBIN}/kraken2 --threads 12 --paired --gzip-compressed --db $KDB $FILTDIR/${read1} $FILTDIR/${read2} --output ${KRES}/${pref}.kraken --classified-out ${KRES}/${pref}#.fastq --use-names --report ${KRES}/${pref}.kraken.report.txt
fi
done



###########################################                                      
#  _                    _              
# | |                  | |             
# | |__  _ __ __ _  ___| | _____ _ __  
# | '_ \| '__/ _` |/ __| |/ / _ \ '_ \ 
# | |_) | | | (_| | (__|   <  __/ | | |
# |_.__/|_|  \__,_|\___|_|\_\___|_| |_|
#                                      
###########################################                                      

#Build bracken DB for 16s RefSeq
THRESHOLD=10
READ_LEN=250

${BRACK_BIN}/bracken-build -d $KDB -t 12 -k 35 -l 250 
for CLASSIFICATION_LEVEL in S G
do
for READ1 in ${READDIR}/*${READPATTERN}
do
read1=$(basename $READ1)
read1=${read1/_R1_001.fastq.gz/}
echo $read1
#Run bracken on kraken results
echo $THRESHOLD
#rm ${KRES}/*_${CLASSIFICATION_LEVEL}.bracken.txt
${BRACK_BIN}/bracken -d ${KDB} -i ${KRES}/${read1}.kraken.report.txt -o ${KRES}/${read1}_${CLASSIFICATION_LEVEL}.bracken.txt -r ${READ_LEN} -l ${CLASSIFICATION_LEVEL} -t ${THRESHOLD}
done
done




################################################
#  _____  _       _           _                     _                      
# |  __ \| |     | |         | |                   | |                     
# | |__) | | ___ | |_    __ _| |__  _   _ _ __   __| | __ _ _ __   ___ ___ 
# |  ___/| |/ _ \| __|  / _` | '_ \| | | | '_ \ / _` |/ _` | '_ \ / __/ _ \
# | |    | | (_) | |_  | (_| | |_) | |_| | | | | (_| | (_| | | | | (_|  __/
# |_|    |_|\___/ \__|  \__,_|_.__/ \__,_|_| |_|\__,_|\__,_|_| |_|\___\___|
#                                                                          
################################################
  

mkdir ${BASEDIR}/04a_filt_plots_RefSeq
mkdir ${BASEDIR}/05a_filt_tables_RefSeq
#Using RefSeq16s db
Rscript ${FUNCDIR}/01_summarize_bracken.r \
    --indir ${KRES} \
    -N 15 --removeme '' \
    --outfile ${BASEDIR}/05a_filt_tables_RefSeq/bracken.txt \
    --outgraph ${BASEDIR}/04a_filt_plots_RefSeq/species_bracken_barplot.pdf

Rscript ${FUNCDIR}/01b_summarize_bracken_genus.r \
    --indir ${KRES} \
    -N 15 --removeme '' \
    --outfile ${BASEDIR}/05a_filt_tables_RefSeq/genus_bracken.txt \
    --outgraph ${BASEDIR}/04a_filt_plots_RefSeq/genus_bracken_barplot.pdf



###############################################################
#  _____               _____  ______  _____            ___  
# |  __ \             |  __ \|  ____|/ ____|          |__ \ 
# | |__) |   _ _ __   | |  | | |__  | (___   ___  __ _   ) |
# |  _  / | | | '_ \  | |  | |  __|  \___ \ / _ \/ _` | / / 
# | | \ \ |_| | | | | | |__| | |____ ____) |  __/ (_| |/ /_ 
# |_|  \_\__,_|_| |_| |_____/|______|_____/ \___|\__, |____|
#                                                   | |     
#                                                   |_|     
#
##################################################################
COND=Status
#Standard

# on 16s_RefSeq
mkdir -p ${BASEDIR}/06a_filt_DA_16s_RefSeq
Rscript ${FUNCDIR}/02_DEseq.r \
--abundance ${BASEDIR}/05a_filt_tables_RefSeq/bracken_raw.txt \
-c $COND -r 200 -O ${BASEDIR}/06a_filt_DA_16s_RefSeq/

mkdir -p ${BASEDIR}/06c_filt_DA_16s_RefSeq_genus
Rscript ${FUNCDIR}/02_DEseq.r \
--abundance ${BASEDIR}/05a_filt_tables_RefSeq/genus_bracken_raw.txt \
-c $COND -r 200 -O ${BASEDIR}/06c_filt_DA_16s_RefSeq_genus/



