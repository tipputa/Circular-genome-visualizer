#$ -S /bin/bash
#$ -cwd
#$ -pe def_slot 4
#$ -l s_vmem=4G -l mem_req=4G 

# please activate python environment
source ~/tool/py37/bin/activate 

python runAfterBlastProcess.py $1 $2
