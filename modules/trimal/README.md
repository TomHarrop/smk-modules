# smk-modules trimal

```python
module trimal:
    snakefile:
        github(
            "tomharrop/smk-modules",
            path="modules/trimal/Snakefile",
            tag="0.7.1",
        )
    config:
        {
            "alignment_directory": Path(test_alignments, "{dataset}"),
            "outdir": Path(output_directory, "{dataset}"),
        }


use rule * from trimal as trimal_*
```

## Workflow

![`snakemake --rulegraph -s
test_workflows/trimal.smk`](../../assets/trimal_graph.svg)

