# MendelianRNA-seq

#### Modification of Beryl Cummings' scripts for discovering novel splicing events through RNA-seq

MendelianRNA-seq is a collection of scripts to help to discover splice sites in a list of bam files. It is a rewrite of Beryl Cummings' original [MendelianRNA-seq](https://github.com/berylc/MendelianRNA-seq) to support parallel processing. It should be used if you want to process 10 bams or less.

[MendelianRNA-seq-DB](https://github.com/dennis-kao/MendelianRNA-seq-DB) is a version of this tool which stores junction information in a database. The main benefit of this is that results from previously processed bam files can be reused and there is a much lower ram usage. The tool also has a few additional features such as fields like total_patient_read_count and n_gtex_seen and exon skipping detection.

## Pipeline

SpliceJunctionDiscovery.py calls upon samtools to report the presence of introns in a list of regions of interest, summarizes their read counts, and writes this to a text file. NormalizeSpliceJunctionValues.py reads this text file, normalizes the read counts of each site based on read support from annotated junctions, and adds this as an additional column to the file. FilterSpliceJunctions.py can then be used to added OMIM annotations, filter out sites that are of low quality and/or are present in a given number of samples, and much more.

SpliceJunctionDiscovery.py is the script that takes the longest to execute. It calls upon samtools based on the number of BAMs * the number of regions of interest. This step is parallelized and the number of worker processes can specified in the torque file or as an arguement to the standalone script. This number should be equal to or less than the number of cores on your system.

NormalizeSpliceJunctions.py and FilterSpliceJunctions.py execute pretty much instantaneously with minimal RAM use.

## Required files

1. .bam (and .bai) files produced from an RNA-seq pipeline - You need a sufficient number of high quality control BAMs so that you can filter out more splice junctions and discover those that are specific to a diseased sample. The [GTEx project](https://www.gtexportal.org/home/) is a good resource for control BAMs. These BAM files should all be from the same tissue due to tissue specific expression. A way to test for contaminated tissue samples has been described in the study cited below. Note that you can generate .bai files from .bam files using this line: ```parallel  samtools index ::: *.bam```

2. transcript_file - A text file containing a list of genes and their spanning chromosome positions that you want to discover junctions in:
	```
	GENE	ENSG	STRAND	CHROM	START	STOP	GENE_TYPE
	```
	You can use [genes.R](https://github.com/naumenko-sa/bioscripts/blob/master/genes.R) for that, or convert an existing .bed file using this bash line:
```
cat kidney.glomerular.genes.bed | awk '{print $4"\t"$4"\t+\t"$1"\t"$2"\t"$3"\tNEXONS"}' >> gene.list
```
3. bamlist.list - A file containing the names of all the bams you want to discover junctions in. The file should quite simply be:
	
	
		G65693.GTEX.8TY6-5R5T.2.bam
		G55612.GTEX.7G54-3GS8.1.bam
		G09321.GTEX.0EYJ-9E12.3.bam
		PATIENT.bam
	
	
	An easy way to generate this file would be to navigate to a directory containing the .bam files you want to use and running this line: ```ls *.bam | grep '' > bamlist.list```

4. transcript_model - A text file containing a list of known canonical splice junctions. These will be used to evaluate a junction's annotation (none, one, both) and it's annotated normalization calculation. You can use your own, or use the [included file](gencode.comprehensive.splice.junctions.txt). This file contains junctions from gencode v19.

## Steps

1. Put bamlist.list, .bam files, .bai files in a seperate directory. Navigate to it. 
	NOTE: there should not be any .txt files present beforehand in order for SpliceJunctionDiscovery.py to run correctly.

2. For [Torque](http://www.adaptivecomputing.com/products/open-source/torque/) users there is a [PBS file](Analysis/rnaseq.novel_splice_junction_discovery.pbs) containing all the commands you need to run. Just change the "home" directory in the file to match where you placed the MendelianRNA-seq folder and run: 

	```qsub MendelianRNA-seq/Analysis/rnaseq.novel_splice_junction_discovery.pbs -v transcriptFile=transcript_file,bamList=bamlist.list,processes=10```
	
	For non-Torque users, the scripts can be run from terminal:
	
	```python3 MendelianRNA-seq/Analysis/SpliceJunctionDiscovery.py -transcriptFile=$transcriptFile -bamList=$bamList -processes=$processes```
	
	```python MendelianRNA-seq/Analysis/NormalizeSpliceJunctionValues.py -transcript_model=$transcript_model -splice_file=$sjdOutput --normalize > $normOutput ```

	Parameters:
	1. transcriptFile, path to file #2
	2. bamList, path to file #3
	3. processes, the number of worker processes running in the background calling samtools. This the slowest step in the program. This number should be equal to or less than the number of cores on your machine. **Keep in mind that increasing the number of worker processes consumes more ram.** Read further on for an explanation on RAM use.
	
		For torque users: This number should also be equal to or less than the number specified for ppn in rnaseq.novel_splice_junction_discovery.pbs:

		
			#PBS -l walltime=10:00:00,nodes=1:ppn=10
	4. sjdOutput, the file produced from SpliceJunctionDiscovery.py. The name of this file is always 'All.' + transcriptFile + '.splicing.list'
	5. transcript_model, path to file #4
	6. normOutput, the name of a file in which NormalizeSpliceJunctionValues.py writes its output to. You can remove this to have the script write to stdout

4. Use FilterSpliceJunctions.py on the file produced by the normalization script. The methododology on how to filter out splice junctions is entirely up to you. You can filter junctions based on read count, specificity to a sample, normalized read count, gene and more. For my research purposes, I often filtered out junctions with a read count less than 5, a normalized read count less than 0.05 and were present in GTEx samples. Having a gene panel, that is a list of genes or regions you suspect has the causative mutation, is very helpful in cutting down the number of splice sites.

	Great documentation on how to use FilterSpliceJunctions.py can be found near the end of this [blog post](https://macarthurlab.org/2017/05/31/improving-genetic-diagnosis-in-mendelian-disease-with-transcriptome-sequencing-a-walk-through/) by Beryl Cumming.

## RAM use in SpliceJunctionDiscovery and runtime

RAM use in SpliceJunctionDiscovery is largely a function of the number of bam files, the number of unqiue splice sites pertaining to a specified gene region in transcript_file (somewhat correlated with coverage) and the number of worker processes. In my experience, with 4 bam files (3 GTEx and 1 patient) and 10 worker processes, ram use leveled out at a constant 35 GB. SpliceJunctionDiscovery completed within 30 minutes. 

The reason for high ram use is that junctions across all samples in a specified gene region are being stored in a single Python dictionary, spliceDict. If you look at the function each worker process goes through, intronDiscovery(), you will notice the loop through each BAM file and the samtools view call to a chromosome region specified in transcript_file.

If RAM use is an issue and/or you are going to process more than 10 BAM files, you should consider using a smaller number of worker processes or use [MendelianRNA-seq-DB](https://github.com/dennis-kao/MendelianRNA-seq-DB) instead.

## Output

The scripts output a single file:

norm-All.**transcript_file**.list, (where transcript_file is the name of require file #2) 

which contains all splice site information pertaining to all samples with a column for normalized read counts.

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
- CIGAR string parsing is handled by a function called parseCIGARForIntrons() whereas before CIGAR strings were handled by piping through multiple bash tools. As a result of improper parsing using bash tools, junction start and/or stop positions were not reported properly (e.x. 1:100-200*1D30 represents an alignment that should really be 1:100-230 or 1:100-231)
- Junction flanking in NormalizeSpliceJunctionValues.py has been fixed and now works. When flanking junctions were added to the set in the original make_annotated_junction_set(), individual characters in the string were added as opposed to the entire string itself (e.x. 1:100-200 gets added as '1', ':', '0', '2', '-')

## Citations

[Improving genetic diagnosis in Mendelian disease with transcriptome sequencing](http://stm.sciencemag.org/content/9/386/eaal5209)

Beryl Cummings' original scripts: [MendelianRNA-seq](https://github.com/berylc/MendelianRNA-seq)

## Footnotes

The included transcript_model file [_gencode.comprehensive.splice.junctions.txt_](https://github.com/dennis-kao/MendelianRNA-seq/blob/master/gencode.comprehensive.splice.junctions.txt) is based off of gencode v19.
