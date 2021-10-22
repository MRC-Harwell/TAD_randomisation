# TAD_randomisation_strategy

This repository contains the scripts required to generate random TADs as outlined in:

<b>Making sense of the linear genome, gene function and TADs</b><br>
Helen S Long, Simon Greenaway, George Powell, Ann-Marie Mallon, Cecilia M Lindgren, Michelle M Simon<br>
doi: https://doi.org/10.1101/2020.09.28.316786
<p>&nbsp;</p>

<b>Scripts:</b><br>
<b>Randomgenome.R:</b> <br>
A script to generate random genomes for mm10.

This r script takes a file containing the coordinates of genes in mm10 ("mm10_proteincodinggenes_biomart_sorted.bed") with the columns (not labelled):
```
<chr> <start> <end> <strand> <gene ID>
```
This file was downloaded from ensembl biomart (V96) and sorted using bedtools. 

The script will generate 100 "random genomes" by randomising the gene ID column within each chromosome. 

The created random genomes can then be overlapped with TADs using bed tools e.g.

```
bedtools intersect -a <TADs.bed> -b <randomgenome.bed> -wao -F 1 > <output>
```
<p>&nbsp;</p>

<b>randomTADsmethod_Arrowhead.py:</b> <br>
A script to generate random TADs with properties similar to Arrowhead TADs as outlined in the paper. It can be run for either human or mouse. 

This python script takes the file "datasetname"_"TADcaller""binsize"kb.txt" with the columns (labelled):
```
<TAD ID> <No. genes in TAD> <Median length of genes in TAD> <chr> <start> <end>
```
It also requires hg19.chrom.sizes/mm10.chrom.sizes as downloaded from UCSC and "mm10_proteincodinggenes_biomart_sorted.bed" formated as above or the equivalent for hg19.

The script takes five arguments: <br>
* Dataset name <br>
* TAD caller <br>
* Species <br>
* Chromosome <br>
* Bin size (kb) <br>

The script is designed to be run as an array job, parallelising per chromsome and the chromsome must be specified such that 1=chr1, 2=chr2 etc. for human 23=chrX and 24=chrY, for mouse 20=chrX and 21=chrY. 


e.g. To run for Bonev ESC at 10kb chr1:

In a folder containing the file "BonevESC_Arrowhead10kb.txt"
```
python randomTADsmethod_Arrowhead.py BonevESC Arrowhead Mouse 1 10
```

Requires python 3, pybedtools, bedtools, random, pandas and numpy.

<p>&nbsp;</p>

<b>randomTADsmethod_TopDom.py:</b> <br>
A script to generate random TADs with properties similar to TopDom TADs as outlined in the paper. It can be run for either human or mouse. 

This python script takes the file "datasetname"_"TADcaller""binsize"kb.txt" with the columns (labelled):
```
<TAD ID> <No. genes in TAD> <Median length of genes in TAD> <chr> <start> <end>
```
It also requires hg19.chrom.sizes/mm10.chrom.sizes as downloaded from UCSC and "mm10_proteincodinggenes_biomart_sorted.bed" formated as above or the equivalent for hg19.

The script takes five arguments: <br>
* Dataset name <br>
* TAD caller <br>
* Species <br>
* Chromosome <br>
* Bin size (kb) <br>

The script is designed to be run as an array job, parallelising per chromsome and the chromsome must be specified such that 1=chr1, 2=chr2 etc. for human 23=chrX and 24=chrY, for mouse 20=chrX and 21=chrY. 


e.g. To run for Bonev ESC at 10kb chr1:

In a folder containing the file "BonevESC_TopDom10kb.txt"
```
python randomTADsmethod_Topdom.py BonevESC TopDom Mouse 1 10
```

Requires python 3, pybedtools, bedtools, random, pandas and numpy.

