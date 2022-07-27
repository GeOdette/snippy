# snippy
![image](https://user-images.githubusercontent.com/96872843/176043639-359d461a-2185-4f68-ae09-be005285ce16.png)
[Launch the workflow at the Latch console](https://console.latch.bio/explore/61353/info)
 [![Twitter URL](https://img.shields.io/twitter/url/https/twitter.com/ge_odette.svg?style=social&label=Follow%20%40ge_odette)](https://twitter.com/ge_odette)
[![Build Status](https://travis-ci.org/tseemann/snippy.svg?branch=master)](https://travis-ci.org/tseemann/snippy)
[![License: MIT](https://img.shields.io/badge/License-GPL%20v2-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)
![Don't judge me](https://img.shields.io/badge/Language-latch-steelblue.svg)

# Snippy
Rapid haploid variant calling and core genome alignment

## Synopsis

Snippy finds SNPs between a haploid reference genome and your NGS sequence
reads.  It will find both substitutions (snps) and insertions/deletions
(indels).  It will use as many CPUs as you can give it on a single computer
(tested to 64 cores).  It is designed with speed in mind, and produces a
consistent set of output files in a single folder.  It can then take a set
of Snippy results using the same reference and generate a core SNP alignment
(and ultimately a phylogenomic tree).

# Calling SNPs

## Input Requirements
* a reference genome in FASTA or GENBANK format (can be in multiple contigs)
* sequence read file(s) in FASTQ or FASTA format (can be .gz compressed) format
* a folder to put the results in

## Output Files

Extension | Description
----------|--------------
.tab | A simple [tab-separated](http://en.wikipedia.org/wiki/Tab-separated_values) summary of all the variants
.csv | A [comma-separated](http://en.wikipedia.org/wiki/Comma-separated_values) version of the .tab file
.html | A [HTML](http://en.wikipedia.org/wiki/HTML) version of the .tab file
.vcf | The final annotated variants in [VCF](http://en.wikipedia.org/wiki/Variant_Call_Format) format
.bed | The variants in [BED](http://genome.ucsc.edu/FAQ/FAQformat.html#format1) format
.gff | The variants in [GFF3](http://www.sequenceontology.org/gff3.shtml) format
.bam | The alignments in [BAM](http://en.wikipedia.org/wiki/SAMtools) format. Includes unmapped, multimapping reads. Excludes duplicates.
.bam.bai | Index for the .bam file
.log | A log file with the commands run and their outputs
.aligned.fa | A version of the reference but with `-` at position with `depth=0` and `N` for `0 < depth < --mincov` (**does not have variants**)
.consensus.fa | A version of the reference genome with *all* variants instantiated
.consensus.subs.fa | A version of the reference genome with *only substitution* variants instantiated
.raw.vcf | The unfiltered variant calls from Freebayes
.filt.vcf | The filtered variant calls from Freebayes
.vcf.gz | Compressed .vcf file via [BGZIP](http://blastedbio.blogspot.com.au/2011/11/bgzf-blocked-bigger-better-gzip.html)
.vcf.gz.csi | Index for the .vcf.gz via `bcftools index`)

:warning: :x: Snippy 4.x does **NOT** produce the following files that Snippy 3.x did

Extension | Description
----------|--------------
.vcf.gz.tbi | Index for the .vcf.gz via [TABIX](http://bioinformatics.oxfordjournals.org/content/27/5/718.full)
.depth.gz | Output of `samtools depth -aa` for the `.bam` file
.depth.gz.tbi | Index for the `.depth.gz` file

## Columns in the TAB/CSV/HTML formats

Name | Description
-----|------------
CHROM | The sequence the variant was found in eg. the name after the ```>``` in the FASTA reference
POS | Position in the sequence, counting from 1
TYPE | The variant type: snp msp ins del complex
REF | The nucleotide(s) in the reference
ALT | The alternate nucleotide(s) supported by the reads
EVIDENCE | Frequency counts for REF and ALT

If you supply a Genbank file as the `--reference` rather than a FASTA
file, Snippy will fill in these extra columns by using the genome annotation
to tell you which feature was affected by the variant:

Name | Description
-----|------------
FTYPE | Class of feature affected: CDS tRNA rRNA ...
STRAND | Strand the feature was on: + - .
NT_POS | Nucleotide position of the variant withinthe feature / Length in nt
AA_POS | Residue position / Length in aa (only if FTYPE is CDS)
LOCUS_TAG | The `/locus_tag` of the feature (if it existed)
GENE | The `/gene` tag of the feature (if it existed)
PRODUCT | The `/product` tag of the feature (if it existed)
EFFECT | The `snpEff` annotated consequence of this variant (ANN tag in .vcf)

## Variant Types

Type | Name | Example
-----|------|-------------
snp  | Single Nucleotide Polymorphism |  A => T
mnp  | Multiple Nuclotide Polymorphism | GC => AT
ins  | Insertion | ATT => AGTT
del  | Deletion | ACGG => ACG
complex | Combination of snp/mnp | ATTC => GTTA

## The variant caller

The variant calling is done by
[Freebayes](https://github.com/ekg/freebayes).
The key parameters under user control are:

* `--mincov` - the minimum number of reads covering a site to be considered (default=10)
* `--minfrac` - the minimum proportion of those reads which must differ from the reference
* `--minqual` - the minimum VCF variant call "quality" (default=100)

## Looking at variants in detail with `snippy-vcf_report`

If you prompt the workflow to generate a report,you will get a snps.report.txt that looks like below:
```
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
>LBB_contig000001:10332 snp A=>T DP=7 Q=66.3052 [7]

         10301     10311     10321     10331     10341     10351     10361
tcttctccgagaagggaatataatttaaaaaaattcttaaataattcccttccctcccgttataaaaattcttcgcttat
........................................T.......................................
,,,,,,  ,,,,,,,,,,,,,,,,,,,,,t,,,,,,,,,,t,,t,,,,,,,,,,,,,,,,g,,,,,,,g,,,,,,,,,t,
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,, .......T..................A............A.......
.........................A........A.....T...........    .........C..............
.....A.....................C..C........CT.................TA.............
,a,,,,,a,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,t,t,,,g,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
,,,,,ga,,,,,,,c,,,,,,,t,,,,,,,,,,g,,,,,,t,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
                            ............T.C..............G...............G......
                                                    ,,,,,,,g,,,,,,,,g,,,,,,,,,,,
                                                           g,,,,,,,,,,,,,,,,,,,,
```
Character | Meaning
----------|-----------
`ATGC`    | Same as the reference
`atgc`    | Different from the reference
`-`       | Zero coverage in this sample **or** a deletion relative to the reference
`N`       | Low coverage in this sample (based on `--mincov`)
`X`       | Masked region of reference (from `--mask`)
`n`       | Heterozygous or poor quality genotype  (has `GT=0/1` or `QUAL < --minqual` in `snps.raw.vcf`)


## Finding SNPs between contigs

The workflow has been configured to work with contig files. Check the contigs options to run. 
Sometimes one of your samples is only available as contigs, without
corresponding FASTQ reads. It does this by shredding the contigs
into 250 bp single-end reads at `2 &times; --mincov` uniform coverage.

## Example output
![snp](https://user-images.githubusercontent.com/96872843/181271654-d71fb361-ac14-4fc6-8f66-b760d161369b.png)

