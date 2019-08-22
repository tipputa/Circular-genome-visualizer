#$ -S /bin/bash
#$ -cwd
#$ -pe def_slot 1 
#$ -o /dev/null
#$ -e /dev/null
#$ -pe def_slot 4

binDir=$1
targetId=$2
FastaDir=$3
DBDir=$4
OutputDir=$5
AllInput=$6
NSLOTS=4

module load singularity

for id in `echo $AllInput| tr ',' '\n'`
do
singularity exec /usr/local/biotools/b/blast\:2.7.1--boost1.64_1 blastp -evalue 1e-05  -outfmt 5 -num_alignments 1 \
-db $DBDir/$id.fa \
-query $FastaDir/$targetId.fa \
-out $OutputDir/${targetId}vs${id}.xml
done

