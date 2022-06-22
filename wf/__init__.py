"""
Assemble and sort some COVID reads...
"""

import subprocess
from pathlib import Path

from latch import small_task, workflow
from latch.types import LatchFile, LatchDir
from typing import Optional
import os


@small_task
def snippy_task(ref_gen: LatchFile, seq_read1: LatchFile, seq_read2: LatchFile, out_name: str, out_dir: LatchDir, cores: int = 16) -> LatchDir:

    # Defining outputs

    local_dir = "/root/snippy_output/"
    local_prefix = os.path.join(local_dir, out_name)

    # Define the command
    _snippy_cmd = [
        "snippy",
        "--cpus",
        str(cores),
        "--outdir",
        str(local_prefix),
        "--ref",
        ref_gen.local_path,
        "--R1",
        seq_read1.local_path,
        "--R2",
        seq_read2.local_path,
    ]
    # Run subprocess
    subprocess.run(_snippy_cmd, check=True)

    # Return out_dir
    return LatchDir(local_dir, out_dir.remote_path)


@workflow
def snippy(ref_gen: LatchFile, seq_read1: LatchFile, seq_read2: LatchFile, out_name: str, out_dir: LatchDir, cores: int = 16) -> LatchDir:
    """

    _Snippy_
    ----

    Latch implementation of Snippy, Rapid haploid variant calling and core genome alignment

    __metadata__:
        display_name: Rapid haploid variant calling and core genome alignment

        author:
            name: SGodette

            email: steveodettegeorge@gmail.com

            github:
        repository: https://github.com/tseemann/snippy.git

        license:
            id: MIT

    Args:

        ref_gen:
          Reference genome. E.g Listeria.gbk

          __metadata__:
            display_name: Reference genome

        seq_read1:
          First read sequence. Tip* Input reference genome in FASTA or GENBANK format

          __metadata__:
            display_name: Read 1

        seq_read2:
          The second read sequence. Tip* Input reference genome in FASTA or GENBANK format

          __metadata__:
            display_name: Read 2

        out_dir:
          Directory for outputs. Tip* Create a directory at the latch console

          __metadata__:
            display_name: Outout directory


        out_name:
          Preffered output prefix

          __metadata__:
            display_name: Output Name

        cpus:
          Number of cpus for processing. Tip* Default at 10

          __metadata__:
            display_name: CPUS

    """

    return snippy_task(ref_gen=ref_gen, seq_read1=seq_read1, seq_read2=seq_read2, out_name=out_name, out_dir=out_dir, cores=cores)
