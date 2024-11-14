# smk-modules captus

Run `captus_assembly assemble`, `captus_assembly extract` and `"captus_assembly
align`, and get some summary stats.

```python
module captus:
    snakefile:
        github(
       "tomharrop/smk-modules",
        path="modules/captus/Snakefile",
        tag="0.7.0",
        )
    config:
        {
            "namelist": Path("test-data", "captus", "namelist.txt"),
            "read_directory": Path(output_directory, "inputs"),
            "target_file": target_file,
            "outdir": Path(output_directory, "plain_captus"),
        }


use rule * from captus as captus_*
```

## Options

```python3
"minimum_sample_wscore": "{minimum_sample_wscore}",
"outgroup": outgroup_samples,  # a list of sample names,
"external_fasta_files": external_fasta_files,
"misc_markers": misc_dna,
```

## Workflow

![`snakemake --rulegraph -s
test_workflows/captus.smk`](../../assets/captus_graph.svg)
