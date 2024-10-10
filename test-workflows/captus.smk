#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import tempfile


# demonstrate how to use a param to find the input directory for e.g. iqtree
def get_alignment_dir(wildcards):
    captus_align_outdir = rules.captus_align.output["outdir"]
    captus_align_output = [x for x in Path(captus_align_outdir).glob("**")]
    trimpath = Path(
        "03_trimmed",
        "06_informed",
        "03_coding_MIT",
        f"{wildcards.alignment_type}",
    )
    untrimpath = Path(
        "02_untrimmed",
        "06_informed",
        "03_coding_MIT",
        f"{wildcards.alignment_type}",
    )
    for path in captus_align_output:
        if path.match(Path("**", trimpath).as_posix()) and path.is_dir():
            return path
    for path in captus_align_output:
        if path.match(Path("**", untrimpath).as_posix()) and path.is_dir():
            return path
    return ""


sample_data = Path("test-data", "hybpiper", "samples.csv")
target_file = Path("test-data", "hybpiper", "combined_targetfiles.fixed.fasta")
read_directory = Path("test-data", "captus", "reads")
output_directory = Path(
    "test-output",
    "captus",
)

# just run on whatever read files we have
all_samples = sorted(
    set(x.stem.split(".")[0] for x in read_directory.glob("*r1.fastq.gz"))
)

# You can use any sample(s) from the dataset as an outgroup to root the
# alignments, but not samples from the target file.
outgroup_samples = all_samples[-2:]

# You can pass external FASTA files to captus extract, e.g. to make them into
# an outgroup. It's much slower (e.g. this sample takes 4h 41m 9.2s on my
# computer) so request a longer allocation for the extract step.
external_fasta_files = [
    Path("test-data", "captus", "GCF_008831285.2_ASM883128v2_genomic.fna")
]

# WHAT HAPPENS IF LOCI ARE DUPLICATED IN THE TARGET FILE?
# NC_045138.2 has the most hits so could use that. First try NC_045147.1, which has a hit for 5858, which is the first locus in the target file
single_chr_fasta_files = [
    # Path("test-data", "captus", "NC_045138.2.fasta")
    Path("test-data", "captus", "NC_045147.1.fasta")
]
target_file_with_duplicated_loci = Path(
    "test-data", "captus", "target_file_with_duplicated_loci.fasta"
)

# You can add miscellaneous DNA markers at the extract stage. Testing with the
# rhizanthella ITS sequences from genbank (retrieved from the nucleotide db by
# searching 'txid158356[Organism:exp]')
misc_dna = [Path("test-data", "captus", "rhizanthella_loci.fasta")]

# captus_snakefile = "../modules/captus/Snakefile"
captus_snakefile = github(
    "tomharrop/smk-modules",
    path="modules/captus/Snakefile",
    tag="0.6.5",
)


captus_alignments = ["01_AA", "02_NT", "03_genes"]


# target is at the end


rule post_captus:
    input:
        Path(output_directory, "plain_captus", "04_alignments"),
    output:
        directory(Path(output_directory, "99_post-captus", "{alignment_type}")),
    params:
        alignment_dir=get_alignment_dir,
    shell:
        "cp -r {params.alignment_dir} {output}"


# NOTE. this does NOT work well with `prefix:` because of the way captus parses
# directory arguments. better to use the outdir where possible.
# Another NOTE. It's very hard to make this work without a checkpoint. If you
# pass the samples as a list, you have to define the module once for each list.
# If you pass the samples in a namelist, you can't start the workflow until
# you've read the list.
module captus:
    snakefile:
        captus_snakefile
    config:
        {
            "namelist": Path("test-data", "captus", "namelist.txt"),
            "read_directory": Path(output_directory, "inputs"),
            "target_file": target_file,
            "outdir": Path(output_directory, "plain_captus"),
        }


use rule * from captus as captus_*


# do this one first, otherwise the wildcard here clobbers the checkpoints in
# other modules. Probably a bug?
module captus_with_wscore_cutoff:
    snakefile:
        captus_snakefile
    config:
        {
            "namelist": Path("test-data", "captus", "namelist.txt"),
            "outdir": Path(
                output_directory, "with_wscore_cutoff.{minimum_sample_wscore}"
            ),
            "read_directory": Path(output_directory, "inputs"),
            "target_file": target_file,
            "minimum_sample_wscore": "{minimum_sample_wscore}",
            "cluster_leftovers": False,
        }


use rule * from captus_with_wscore_cutoff as captus_with_wscore_cutoff_*


module captus_with_outgroup:
    snakefile:
        captus_snakefile
    config:
        {
            "namelist": Path("test-data", "captus", "namelist.txt"),
            "outgroup": outgroup_samples,  # a list of sample names
            "read_directory": Path(output_directory, "inputs"),
            "target_file": target_file,
            "outdir": Path(output_directory, "with_outgroup"),
        }


use rule * from captus_with_outgroup as captus_with_outgroup_*


module captus_with_external:
    snakefile:
        captus_snakefile
    config:
        {
            "namelist": Path("test-data", "captus", "namelist.txt"),
            "outgroup": [x.with_suffix("").name for x in external_fasta_files],
            "external_fasta_files": external_fasta_files,
            "read_directory": Path(output_directory, "inputs"),
            "target_file": target_file,
            "outdir": Path(output_directory, "with_external"),
        }


use rule * from captus_with_external as captus_with_external_*


module captus_with_duplicate_loci:
    snakefile:
        captus_snakefile
    config:
        {
            "namelist": Path("test-data", "captus", "namelist.txt"),
            "external_fasta_files": single_chr_fasta_files,
            "read_directory": Path(output_directory, "inputs"),
            "target_file": target_file_with_duplicated_loci,
            "outdir": Path(output_directory, "with_duplicate_loci"),
        }


use rule * from captus_with_duplicate_loci as captus_with_duplicate_loci_*


module captus_with_misc:
    snakefile:
        captus_snakefile
    config:
        {
            "namelist": Path("test-data", "captus", "namelist.txt"),
            "read_directory": Path(output_directory, "inputs"),
            "target_file": target_file,
            "misc_markers": misc_dna,
            "outdir": Path(output_directory, "with_misc"),
        }


use rule * from captus_with_misc as captus_with_misc_*


rule set_up_read_files:
    input:
        read_directory=Path(read_directory, "{sample}.r{r}.fastq.gz"),
    output:
        temp(Path(output_directory, "inputs", "{sample}.r{r}.fastq.gz")),
    shell:
        "ln -s $( readlink -f {input} ) $( readlink -f {output} )"


rule generate_namelist:
    output:
        Path("test-data", "captus", "namelist.txt"),
    run:
        with open(output[0], "wt") as f:
            f.write("\n".join(all_samples))


rule target:
    default_target: True
    input:
        expand(
            Path(output_directory, "99_post-captus", "{alignment_type}"),
            alignment_type=captus_alignments,
        ),
        rules.captus_with_duplicate_loci_target.input,
        rules.captus_with_external_target.input,
        rules.captus_with_misc_target.input,
        rules.captus_with_outgroup_target.input,
        expand(
            rules.captus_with_wscore_cutoff_target.input,
            minimum_sample_wscore=["0", "0.2", "0.3"],
        ),
