﻿# MendelianRNA-seq

#### Modification of Beryl Cummings scripts for discovering novel splicing events through RNA-seq

MendelianRNA-seq helps to discover splice sites in a sample given a list of bam files. SpliceJunctionDiscovery.py calls upon samtools to report the presence of introns in a list of regions of interest and summarizes these results for read count. NormalizeSpliceJunctionValues.py normalizes the read count of each site based on read support from nearby junctions. FilterSpliceJunctions.py then can be used to added OMIM annotations, filter out sites that are of low quality and/or are present in a given number of samples, and much more.

SpliceJunctionDiscovery.py usually takes the longest to execute because it calls upon samtools based on the number of samples * the number of regions of interest. This step is parallelized and the number of worker processes can specified in the torque file or as an arguement to the standalone script. This number should be equal to or less than the number of cores on your system.

[MendelianRNA-seq-DB](https://github.com/dennis-kao/MendelianRNA-seq-DB) is a version of this tool which stores junction information in a database. The main benefit of this is that results can be reused and do not have to be recomputed for previously processed BAM files. The tool also has a few additional features like a column for total read counts seen in control or patients bams for a specific junction and exon skipping detection. 

## Steps

1. Run bcbio RNA-seq pipeline to get bam files

2. Create a list of genes of interest (muscular or kidney), in the format:
	
	```GENE	ENSG	STRAND	CHROM	START	STOP	NEXONS```

	use [genes.R](https://github.com/naumenko-sa/bioscripts/blob/master/genes.R) for that.

	Some ready list are in data folder.

	```cat kidney.glomerular.genes.bed | awk '{print $4"\t"$4"\t+\t"$1"\t"$2"\t"$3"\tNEXONS"}' >> kidney.glomerular.genes.list```

3. Make and/or navigate to a directory containing all your .bam and corresponding .bai files. Run the [novel splice junction discovery script](https://github.com/dennis-kao/MendelianRNA-seq/blob/master/Analysis/rnaseq.novel_splice_junction_discovery.pbs). NOTE: there should not be any .txt files present beforehand in order for SpliceJunctionDiscovery.py to run correctly.

	```qsub MendelianRNA-seq/Analysis/rnaseq.novel_splice_junction_discovery.pbs -v transcriptFile=kidney.glomerular.genes.list,bamList=bamlist.list,sample=sampName```

	Mandatory parameters:
	1. transcriptFile, path to file produced in step 2
	2. bamList, a text file containing the names of all bam files used in the analysis, each on a seperate line. For example:

		```
		control1.bam
		control2.bam
		control3.bam
		findNovel.bam
		```
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

** SpliceJunctionDiscovery has been entirely rewritten in Python and a small bug in NormalizeSpliceJunctionValues.py has been fixed. No other scripts have been modified. SpliceJunctionDiscovery still outputs a text file that works with Beryl's NormalizeSpliceJunctionValues.py and FilterSpliceJunctions.py scripts. **

- SpliceJunctionDiscovery has been rewritten in Python and parallelized - decreasing processing time by a factor proprotional to the number of worker processes
- CIGAR string parsing is handled by a function called parseCIGARForIntrons() whereas before CIGAR strings were handled by piping through multiple bash tools. As a result of improper parsing using bash tools, junction start and/or stop positions were not reported properly (i.e. 1:100-200*1D30 represents an alignment that should really be 1:100-230 or 1:100-231)
- Junction flanking in NormalizeSpliceJunctionValues.py has been fixed and now works. When flanking junctions were added to the set in make_annotated_junction_set(), individual characters in the string were added as opposed to the entire string itself (i.e. 1:100-200 gets added as '1', ':', '0', '2', '-')

## Footnotes

The included transcript_model file [_gencode.comprehensive.splice.junctions.txt_](https://github.com/dennis-kao/MendelianRNA-seq/blob/master/gencode.comprehensive.splice.junctions.txt) is based off of gencode v19.
