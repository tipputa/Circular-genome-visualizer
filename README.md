## About
Circular visualizer for complete genomes

Input: GenBank files (.gb)

Usage: 
```sh
python runAllProcess.py <output directory> <input directory (GenBank files)>
```

## Requirements
 - BLAST+ (ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
 - Circos (http://circos.ca/)
 - python 3 & modules
   - `pandas`, `numpy`, `biopython`, `tqdm`	
 
#### PATH
You need to add the blast+/bin and the circos/bin to your PATH.
Please check below.
```sh
echo $PATH
blastn -h
circos -h
```
## Run
#### Run all processes
```sh
python runAllProcess.py <output directory> <input directory (GenBank files)>
```
Two required arguments are as follows:
<ol>
<li> path to the output directory </li>
<li> path to the directory containing GenBank files (.gb) </li>
</ol>

#### Run after the blastp processes
```sh
python runAfterBlastProcess.py <output directory> <input directory (GenBank files)>
```
Two required arguments are as follows:
<ol>
<li> path to the output directory </li>
<li> path to the directory containing GenBank files (.gb) </li>
</ol>

#### Run visualization
```sh
python runVisualize.py  <output directory> <configuration file> <option; key word for output; default:"test"> <option; the minimum number of genes in each cluster; default: 1> <option; sorting column name; default: None>
```
Five arguments are as follows:
<ol>
<li> path to the output directory </li>
<li> path to the configuration file. Please see below.</li>
<li> optional: suffix for the output image file. Default is "test". </li>
<li> optional: the minimum number of genes in each cluster for visualization. Default is 1. </li>
<li> optional: the column name in the configuration file for sorting. Default is None. </li>
</ol>

#### Configuration file
This file will be outputed as "RingOrder_*_df.tsv" by runAllProcess.py and runAfterBlastProcess.py. Please see ./testResult/RingOrder_aligned_df.tsv. and ./testResult/changed_setting.tsv for examples.
<table>
<tr>
  <td>AccNo</td>
  <td>Genome_size</td>
  <td>Strand</td>  
  <td>Angle</td>  
  <td>Deviation (Aligned)</td>  
  <td>Deviation (Original)</td>  
  <td>optional</td>  
</tr>
<tr>
  <td>NC_000915.1</td>
  <td>1667867</td>
  <td>0</td>  
  <td>0</td>  
  <td>53.311</td>  
  <td>47.061</td>  
  <td>...</td>  
</tr>
<tr>
  <td>NC_014256.1</td>
  <td>1673997</td>
  <td>1</td>  
  <td>342</td>  
  <td>53.07</td>  
  <td>177.67</td>  
  <td>...</td>  
</tr>
<tr>
  <td>...</td>
  <td>...</td>
  <td>...</td>  
  <td>...</td>  
  <td>...</td>  
  <td>...</td>  
  <td>...</td>  
</tr>
</table> 

You can edit the visualization result, such as the number of genomes and the ring order, by deleting / reordering rows in this file.
 - AccNo: GenBank accession number
 - Genome_size: the genome size (bp)
 - Strand: choice of the strand for alignment. 0 = original, 1 = complement
 - Angle: rotation in degree for alignment. (-20)–359
 - Deviation (Aligned): The deviation to the consensus after alignment
 - Deviation (Original): The deviation to the consensus before alignment
 - optional: you can add the new column for sorting the ring order.


## Test run
#### Test data: 5 Helicobacter pylori genomes
#### all processes
It takes 10 minutes (BLASTP 8 min, other 2 min) on a standalone desktop server of 32GB memory.
```sh
cd
git clone git@github.com:tipputa/Circular-genome-visualizer.git
python ~/Circular-genome-visualizer/bin/runAllProcess.py ~/Circular-genome-visualizer/test/ ~/Circular-genome-visualizer/test/gb/
```


#### RunVisualize. Configuring the visualization.
In this example, "changed_setting.tsv" is a modified configuring file, where the first row was deleted from /test/data/RingOrder_aligned_df.tsv.
```sh
python ~/Circular-genome-visualizer/bin/runVisualize.py ~/Circular-genome-visualizer/test/ ~/Circular-genome-visualizer/test/changed_setting.tsv "rm1genome" 4
```
"rm1genome" is a suffix for the output file (e.g. circos_rm1genome.png). The genes conserved in >= 4 genomes are visualized.


## Citation
 - Tada I, Tanizawa Y, Arita M. [Visualization of consensus genome structure without using a reference genome.](http://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-017-3499-7) BMC Genomics. 2017;18(Suppl 2):208. doi:10.1186/s12864-017-3499-7.
 
