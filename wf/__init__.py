"""
Rapid haploid variant calling and core genome alignment
"""

import subprocess
from pathlib import Path
from latch import small_task, workflow, large_gpu_task
import os
from enum import Enum
from typing import Optional
from latch.types import LatchFile, LatchDir, LatchParameter, LatchAuthor, LatchMetadata
from latch.resources.launch_plan import LaunchPlan


metadata = LatchMetadata(
    display_name="Variant call with snippy",
    documentation="https://github.com/tseemann/snippy",
    author=LatchAuthor(
        name="GeOdette",
        email="steveodettegeorge@gmail.com",
        github="https://github.com/GeOdette",
    ),
    repository="https://github.com/GeOdette/snippy",
    license="MIT"
)
metadata.parameters = {
    "read_type": LatchParameter(
        display_name="Which files do you wish to process",
        description="Click on the reads you wish to process, eighter paired end, single end, of a contigs file"
    ),
    "ref_gen": LatchParameter(
        display_name="Reference genome",
        description="Reference genome. E.g Listeria.gbk",
    ),
    "seq_read1": LatchParameter(
        display_name="Read 1",
        description="First read sequence.",
    ),
    "seq_read2": LatchParameter(
        display_name="Read 2",
        description="The second read sequence.",
    ),
    "out_dir": LatchParameter(
        display_name="Output directory",
        description="Directory for outputs. Tip * Create a directory at the latch console",
    ),
    "cpus": LatchParameter(
        display_name="CPUS",
        description="Number of cpus for processing. Tip * Default at 10"
    ),
    "gen_report": LatchParameter(
        display_name="Generate SNP report",
        description="Check this option to generate snippy-vcf_report. Checking this option increase run time to about 2 hours. Just check the option and go for lunch"
    ),
    "contigs_file": LatchParameter(
        display_name="Contigs file",
        description="The contigs file"
    ),
    "cores": LatchParameter(
        display_name="CPUS",
        description="CPUS to run the samples. User higher CPUS for many samples"
    )

}


class ReadType(Enum):
    SE = "single end read"
    PE = "Paired end reads"
    CONT = "contigs"


def malformedrun():
    pass


def NoContigsFileFound(Exception):
    pass


@large_gpu_task
def snippy_task(ref_gen: LatchFile,
                seq_read1: Optional[LatchFile],
                seq_read2: Optional[LatchFile],
                contigs_file: Optional[LatchFile],
                out_dir: LatchDir,
                read_type: ReadType = ReadType.PE,
                gen_report: bool = False,
                cores: str = "16") -> LatchDir:

    # Defining outputs

    local_dir = "/root/snippy_output/"
    local_prefix = os.path.join(local_dir, "Latch_snippy_files")

    # Lets define command options
    if read_type.value == "Paired end reads":
        if gen_report == True:
            reads = ["--report", "--outdir", str(local_prefix), "--ref",
                     ref_gen.local_path, "--R1", seq_read1.local_path, "--R2", seq_read2.local_path, ]
        elif gen_report != True:
            reads = ["--outdir", str(local_prefix), "--ref", ref_gen.local_path,
                     "--R1", seq_read1.local_path, "--R2", seq_read2.local_path, ]
    elif read_type.value == "single end read":
        if gen_report == True:
            reads = ["--report", "--outdir",
                     str(local_prefix), "--ref", ref_gen.local_path, "--R1", seq_read1.local_path, ]
        elif gen_report != True:
            reads = ["--outdir",
                     str(local_prefix), "--ref", ref_gen.local_path, "--R1",  seq_read1.local_path, ]
    elif read_type.value == "contigs":
        if contigs_file is None:
            raise NoContigsFileFound(
                print(
                    "No contigs file found. Please input a contigs file and try again")
            )
        elif contigs_file is not None:
            if gen_report == True:
                reads = ["--report", "--outdir", str(local_prefix),
                         "--ref", ref_gen.local_path, "--contigs", contigs_file.local_path]
            elif gen_report != True:
                reads = ["--outdir", str(local_prefix), "--ref",
                         ref_gen.local_path, "--contigs", contigs_file.local_path]

    # Define the command
    _snippy_cmd = [
        "snippy",
        "--cpus",
        str(cores),
        *reads

    ]

    # Run subprocess
    subprocess.run(_snippy_cmd, check=True)
    return LatchDir(local_dir, out_dir.remote_path)


@ workflow(metadata)
def snippy(ref_gen: LatchFile,
           seq_read1: Optional[LatchFile],
           seq_read2: Optional[LatchFile],
           contigs_file: Optional[LatchFile],
           out_dir: LatchDir,
           read_type: ReadType = ReadType.PE,
           gen_report: bool = False,
           cores: str = "16") -> LatchDir:
    """Rapid haploid variant calling and core genome alignment

    Snippy
    ----

    Latch implementation of Snippy, Rapid haploid variant calling and core genome alignment

    # Basic usage:

    The workflow will require the following inputs:

    > A reference genome in FASTA or GENBANK format (can be in multiple contigs)<br>
    > A sequence read file(s) in FASTQ or FASTA format (can be .gz compressed) format<br>
    > A folder to put the results in

    # Outputs

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
    .vcf.gz.tbi | Index for the .vcf.gz via [TABIX](http://bioinformatics.oxfordjournals.org/content/27/5/718.full)
    .depth.gz | Output of `samtools depth -aa` for the `.bam` file
    .depth.gz.tbi | Index for the `.depth.gz` file

    # Columns in the TAB/CSV/HTML formats

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

    # Variant Types

    Type | Name | Example
    -----|------|-------------
    snp  | Single Nucleotide Polymorphism |  A => T
    mnp  | Multiple Nuclotide Polymorphism | GC => AT
    ins  | Insertion | ATT => AGTT
    del  | Deletion | ACGG => ACG
    complex | Combination of snp/mnp | ATTC => GTTA

    # Example output
    ![snp](https://user-images.githubusercontent.com/96872843/181271654-d71fb361-ac14-4fc6-8f66-b760d161369b.png)




    """

    return snippy_task(ref_gen=ref_gen,
                       seq_read1=seq_read1,
                       seq_read2=seq_read2,
                       contigs_file=contigs_file,
                       read_type=read_type,
                       out_dir=out_dir,
                       gen_report=gen_report,
                       cores=cores)
