# smk-modules iqtree

```python
module iqtree:
    snakefile:
        iqtree_snakefile = github(
            "tomharrop/smk-modules",
            path="modules/iqtree/Snakefile",
            tag="0.2.13",
        )
    config:
        {
            "alignment_directory": Path(test_alignments, "{dataset}"),
            "outdir": Path(output_directory, "{dataset}"),
        }

use rule * from iqtree as iqtree_*
```

## Workflow

![`snakemake --rulegraph -s
test_workflows/iqtree.smk`](../../assets/iqtree_graph.svg)

