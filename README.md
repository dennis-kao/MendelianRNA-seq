# MendelianRNA-seq

#### Modification of Beryl Cummings scripts for discovering novel splicing events through RNA-seq

MendelianRNA-seq is a collection of scripts to help to discover splice sites in a list of bam files. 

[MendelianRNA-seq-DB](https://github.com/dennis-kao/MendelianRNA-seq-DB) is a version of this tool which stores junction information in a database. The main benefit of this is that results can be reused and do not have to be recomputed for previously processed BAM files. The tool also has a few additional features like a column for total read counts seen in control or patients bams for a specific junction and exon skipping detection. 

## Pipeline

SpliceJunctionDiscovery.py calls upon samtools to report the presence of introns in a list of regions of interest, summarizes their read counts, and writes this to a text file. NormalizeSpliceJunctionValues.py reads this text file, normalizes the read counts of each site based on read support from annotated junctions, and adds this as an additional column to the file. FilterSpliceJunctions.py can then be used to added OMIM annotations, filter out sites that are of low quality and/or are present in a given number of samples, and much more.

SpliceJunctionDiscovery.py is the script that takes the longest to execute. It calls upon samtools based on the number of BAMs * the number of regions of interest. This step is parallelized and the number of worker processes can specified in the torque file or as an arguement to the standalone script. This number should be equal to or less than the number of cores on your system.

## Required files

1. BAM (and BAI) files produced from an RNA-seq pipeline - You need a sufficient number of high quality control BAMs so that you can filter out more splice junctions and discover those that are specific to a diseased sample. The [GTEx project](https://www.gtexportal.org/home/) is a good resource for control BAMs. These BAM files should all be from the same tissue due to tissue specific expression. A way to test for contaminated tissue samples has been described in the study cited below. You can also generate .bai files from .bam files using this line: ```parallel  samtools index ::: *.bam```

2. gene_list - A text file containing a list of genes and their spanning chromosome positions that you want to discover junctions in:
	```
	GENE	ENSG	STRAND	CHROM	START	STOP	GENE_TYPE
	```
You can use [genes.R](https://github.com/naumenko-sa/bioscripts/blob/master/genes.R) for that, or convert an existing .bed file using this bash line:
```
cat kidney.glomerular.genes.bed | awk '{print $4"\t"$4"\t+\t"$1"\t"$2"\t"$3"\tNEXONS"}' >> gene.list
```

3. bamlist - A text file containing the names of all the bams you want to discover junctions in. The file should quite simply be:

		```
		G65693.GTEX.8TY6-5R5T.2.bam
		G55612.GTEX.7G54-3GS8.1.bam
		G09321.GTEX.0EYJ-9E12.3.bam
		PATIENT.bam
		```
An easy way to generate this file would be to navigate to a directory containing the .bam files you want to use and running this line: ```ls *.bam | grep '' > bamlist.list```

4. transcript_model - A text file containing a list of junctions. These will be used to evaluate a junction's annotation (none, one, both) and it's normalization calculation. You can use your own, or use the [included file](gencode.comprehensive.splice.junctions.txt). This file is based off of junctions from gencode v19.

## Steps

1. Get all the files you need as listed above and put them in a seperate directory. Navigate to it. NOTE: there should not be any .txt files present beforehand in order for SpliceJunctionDiscovery.py to run correctly.

2. For [Torque](http://www.adaptivecomputing.com/products/open-source/torque/) users there is a [PBS file](Analysis/rnaseq.novel_splice_junction_discovery.pbs) containing all the commands you need to run. Just change the "home" directory in the file to match where you placed the MendelianRNA-seq folder and run: 

	```qsub MendelianRNA-seq/Analysis/rnaseq.novel_splice_junction_discovery.pbs -v transcriptFile=kidney.glomerular.genes.list,bamList=bamlist.list,sample=sampName```

	Parameters:
	1. transcriptFile, path to file produced in step 2
	2. bam_list
	3. sample, the name of the bam file you want to find novel junctions in, without the ".bam" extension. For example, if your file name is "findNovel.bam", then write "sample=findNovel"

	Optional parameters:
	1. minread, the minimum number of reads a junction needs to have (default=10)
	2. threshold, the minimum normalized read count a junction needs to have (default=0.5)
	3. [transcript_model](https://github.com/dennis-kao/MendelianRNA-seq/blob/master/gencode.comprehensive.splice.junctions.txt), the absolute path to a text file containing a list of known canonical splice junctions (default=/home/dennis.kao/tools/MendelianRNA-seq/gencode.comprehensive.splice.junctions.txt). This file is used in NormalizeSpliceJunctionValues.py.
	4. processes, the number of worker processes running in the background calling samtools. This the slowest step in the program. This number should be equal to or less than the number of cores on your machine. 
	
		For torque users: This number should also be equal to or less than the number specified for ppn in rnaseq.novel_splice_junction_discovery.pbs:

		
			#PBS -l walltime=10:00:00,nodes=1:ppn=10

4. (Optional) By default the pipeline automatically generates a file with filtered results based on the following torque parameters: sample, threshold and minread. The filtered file only has splice sites which are specific to one given sample (sample), have a minimum normalized read count (threshold) and have a minimum read count (minread). However, as a researcher, you may look to filter your dataset in multiple ways according to coverage, diagnosis, etc. 

	For a more comprehensive tool, use FilterSpliceJunctions.py on the file produced by the normalization script. Great documentation on how to use the script can be found near the end of this [blog post](https://macarthurlab.org/2017/05/31/improving-genetic-diagnosis-in-mendelian-disease-with-transcriptome-sequencing-a-walk-through/). 

## Output

The scripts output 2 files:

norm-All.**kidney.glomerular.genes**.list, (where kidney.glomerular.genes is the name of your transcriptFile) 

which contains all splice site information pertaining to all samples with a column for normalized read counts,

and

threshold**X.XX**\_novel\_**sampleName**\_norm\_**All.kidney.glomerular.genes.list**.splicing.txt, (where X.XX is the threshold value, sampleName is the sample you want to discover novel junctions in, and All.kidney.glomerular.genes.list is the name of the input file) 

which contains splice sites only seen in sampleName that have a read count > minRead and a normalized read count > threshold.

The "threshold" file contains text information in the format:

```GENE	GENE_TYPE	CHROM:START-STOP	READ_COUNT	SAMPLES_SEEN	READ_COUNT:SAMPLE	SITES_ANNOTATED	NORM_READ_COUNT:SAMPLE```

Here is a sample output:

```
WDPCP	NEXON	2:64040820-64054790	23	1	23:SAMPLE	Neither annotated	-
YIPF7	NEXON	4:44631580-44637960	21	1	21:SAMPLE	Neither annotated	-
ZNF331	NEXON	19:54041790-54042480	15	1	15:SAMPLE	One annotated	15.0:SAMPLE
ZNF404	NEXON	19:44384341-44388442	40	1	40:SAMPLE	Neither annotated	-
```
## Differences between this MendelianRNA-seq and Beryl Cumming's original MendelianRNA-seq

**SpliceJunctionDiscovery.sh and SpliceJunctionSummary.py have been rewritten and condensed in a single script called SpliceJunctionDiscovery.py. A small bug in NormalizeSpliceJunctionValues.py has been fixed. No other scripts have been modified. SpliceJunctionDiscovery still outputs a text file that works with Beryl's NormalizeSpliceJunctionValues.py and FilterSpliceJunctions.py scripts.**

- SpliceJunctionDiscovery has been rewritten in Python and parallelized - decreasing processing time by a factor proprotional to the number of worker processes
- CIGAR string parsing is handled by a function called parseCIGARForIntrons() whereas before CIGAR strings were handled by piping through multiple bash tools. As a result of improper parsing using bash tools, junction start and/or stop positions were not reported properly (i.e. 1:100-200*1D30 represents an alignment that should really be 1:100-230 or 1:100-231)
- Junction flanking in NormalizeSpliceJunctionValues.py has been fixed and now works. When flanking junctions were added to the set in make_annotated_junction_set(), individual characters in the string were added as opposed to the entire string itself (i.e. 1:100-200 gets added as '1', ':', '0', '2', '-')

## Citations

[Improving genetic diagnosis in Mendelian disease with transcriptome sequencing](http://stm.sciencemag.org/content/9/386/eaal5209)

Beryl Cumming's original scripts: [MendelianRNA-seq](https://github.com/berylc/MendelianRNA-seq)

## Footnotes

The included transcript_model file [_gencode.comprehensive.splice.junctions.txt_](https://github.com/dennis-kao/MendelianRNA-seq/blob/master/gencode.comprehensive.splice.junctions.txt) is based off of gencode v19.
