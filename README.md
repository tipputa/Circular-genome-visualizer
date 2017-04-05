## About
Complete circular genome visualizer.

Input: GenBank files (.gb)

Usage: 
```sh
python runAllProcess.py <output directory> <input directory (GenBank files)>
python runAfterBlastProcess.py <output directory> <input directory (GenBank files)>
python runVisualize.py  <output directory> <data info file> <option; key word for output; default:"test"> <option; the minimum number of genes in each cluster; default: 1> <option; sorting column name; default: None>
```

## Requirements
 - BLAST+ (ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
 - Circos (http://circos.ca/)
 - python 3 & modules
   - `pandas`, `numpy`, `biopython`, `tqdm`	
 
#### PATH
You need to add the blast+/bin and the circus/bin to your PATH.
Please check below.
```sh
echo $PATH
blastn -h
circos -h
```
## Run
#### Run all process
```sh
python runAllProcess.py <output directory> <input directory (GenBank files)>
```
There are two arguments as follows:
<ol>
<li> path to output directory </li>
<li> path to directory containing GenBank files (.gb) </li>
</ol>

#### Run after blastp process
```sh
python runAfterBlastProcess.py <output directory> <input directory (GenBank files)>
```
There are two arguments as follows:
<ol>
<li> path to output directory </li>
<li> path to directory containing GenBank files (.gb) </li>
</ol>

#### Run visualize
```sh
python runVisualize.py  <output directory> <setting file> <option; key word for output; default:"test"> <option; the minimum number of genes in each cluster; default: 1> <option; sorting column name; default: None>
```
There are five arguments as follows:
<ol>
<li> path to output directory </li>
<li> path to the setting file; please check seetting file description</li>
<li> optional: key words for output files; default: "test" </li>
<li> optional: the  minimum number of genes in each cluster for visualization; default: 1 </li>
<li> optional: the column name in the setting file for sorting. default: None </li>
</ol>

#### Setting file
This file will be outputed as "RingOrder_*_df.tsv". Please check ./testResult/RingOrder_aligned_df.tsv.
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

You can modify the setting for visualization, such as the number of genomes and the ring order.
Please delete or reorder any rows where you want.
 - AccNo: GenBank accession number
 - Genome_size: the genome size (bp)
 - Strand: choice of strand for alignment;  0 = original, 1 = complement
 - Angle: rotation degree for alignment; (-20)–359
 - Deviation (Aligned): The deviation to the consensus after alignment
 - Deviation (Original): The deviation to the consensus before alignment
 - optional: you can add the new column, and set the fifth argument.

## Test run
### - Five Helicobacter pylori genomes
#### all process
It will take 10 minutes (BLASTP 8 min, other 2 min).
```sh
cd
git clone git@github.com:tipputa/Circular-genome-visualizer.git
cd Circular-genome-visualizer
python bin/runAllProcess.py ~/Circular-genome-visualizer/test/ ~/Circular-genome-visualizer/test/gb/
```

#### RunVisualize. Changing the visualization settings.
We can change the setting file.
In this test case, we can use "changed_setting.tsv" file (first row was deleted from test/data/RingOrder_aligned_df.tsv).
```sh
python3 bin/runVisualize.py ~/Circular-genome-visualizer/test/ ~/Circular-genome-visualizer/test/changed_setting.tsv "rm1genome" 4
```
"rm1genome" is a key word for output files (e.g. circos_rm1genome.png). 
The genes conserved in over 4 genomes are visualized.


## Citation:
 - Tada I, Tanizawa Y, Arita M. [Visualization of consensus genome structure without using a reference genome.](http://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-017-3499-7) BMC Genomics. 2017;18(Suppl 2):208. doi:10.1186/s12864-017-3499-7.
 
