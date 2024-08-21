# Changelog

## 0.5.1 (2024-08-21)

### New

* Run trimal on a folder of alignments.

## 0.4.6 (2024-08-15)

### New

* Generate sample stats from Captus extractions.

## 0.4.5 (2024-08-13)

### Fix

* Bugfix in decruft.

## 0.4.4 (2024-08-12)

### New

* Include miscellaneous markers for captus extract.

## 0.4.3 (2024-08-09)

### New

* Add external reference sequence to captus extract.

## 0.4.02 (2024-08-08)

### New

* Aggressively decruft paragone input.

## 0.4.01 (2024-08-08)

### New

* Add outgroups for Captus align.

## 0.3.04 (2024-08-07)

### Other

* Detect files with two or fewer sequences.

* Changelog.

## 0.3.03 (2024-08-05)

### Fix

* Check for empty files in orient_sequences.

### Other

* Changelog.

## 0.3.02 (2024-07-27)

### Fix

* Find cluster file inside temporary path script.

* Find cluster file inside temporary path script.

### Other

* Update changelog.

* Changelog.

## 0.3.01 (2024-07-25)

### New

* Recover plastid markers with Captus.

## 0.2.18 (2024-07-25)

### Other

* Log what happened during orientation.

## 0.2.17 (2024-07-24)

### Other

* Multiprocess fasta files.

## 0.2.13 (2024-07-24)

### Other

* Tags.

* It's impossible to run modules on separate namelists without a checkpoint inthe modules.

* Remove checkpoints from captus and demonstrate how to find output.

* Checkpoints removed from iqtree.

* Initial commit.

## 0.2.12 (2024-07-21)

### Other

* Mark tempfile as ancient.

## 0.2.11 (2024-07-19)

### Other

* Switch paragone to clustalo.

## 0.2.10 (2024-07-16)

### Other

* Paragone version bump.

## 0.2.09 (2024-07-15)

### Other

* Update iqtree nt flag.

## 0.2.08 (2024-07-14)

### Other

* Avoid temporary input directory.

## 0.2.07 (2024-07-14)

### Other

* Use ancient to prevent reruns.

* Debug.

## 0.2.06 (2024-07-12)

### Other

* Refactor paragone to use archives for out and input.

* Petrichor is making recursive symlinks.

## 0.2.05 (2024-07-10)

### Other

* Handle archives.

## 0.2.04 (2024-07-10)

### Other

* Revert captus changes.

* Workaround lack of tmpdir for hybpiper tar.

## 0.2.03 (2024-07-09)

### Other

* Retrieve captus clusters.

* Merge pull request #2 from TomHarrop/expand_captus_wildcards.

  Expand captus wildcards

## 0.2.02 (2024-07-09)

### Other

* Expand workaround.

## 0.2.01 (2024-07-09)

### Other

* Workaround for expand???

* Update test workflows.

## 0.2.0 (2024-07-09)

### Other

* Refactor captus.

* Attempt to use tmpdir.

* Attempt to use tmpdir.

* Attempt to use tmpdir.

* Attempt to use tmpdir.

* Attempt to use tmpdir.

* Attempt to use tmpdir.

* Attempt to use tmpdir.

* Attempt to use tmpdir.

* Attempt to use tmpdir.

* Attempt to use tmpdir.

* Attempt to use tmpdir.

* Attempt to use tmpdir.

* Use shadow rather than tmpdir.

* Use shadow rather than tmpdir.

* Use shadow rather than tmpdir.

* Use shadow rather than tmpdir.

* Attempt to use tmpdir.

* Use shadow rather than tmpdir.

* Use shadow rather than tmpdir.

* Use shadow rather than tmpdir.

* Attempt to use tmpdir.

* Attempt to use tmpdir.

* Attempt to fix tmp issues in modules.

## 0.1.07 (2024-07-04)

### Other

* Iqtree update.

* Update iqtree container.

## 0.1.06 (2024-06-28)

### Other

* Update captus target.

## 0.1.05 (2024-06-27)

### Other

* Update captus target.

* Merge remote-tracking branch 'origin'

## 0.1.03 (2024-06-27)

### Other

* Update captus target.

## 0.1.04 (2024-06-25)

### Other

* Typo in path.

* Use a shadow directory for fasterq-dump.

## 0.1.02 (2024-06-20)

### Other

* Iqtree.

* Iqtree.

## 0.1.01 (2024-06-20)

### Other

* Something wrong with wildcards.

## 0.1.00 (2024-06-19)

### Other

* Exclude groups with <3 seqs from tree.

* Note on problematic alignments.

* Use checkpoints to resolve alignments.

## 0.0.59 (2024-06-19)

### Other

* Format params in iqtree.

## 0.0.58 (2024-06-13)

### Other

* Make directory an input.

## 0.0.57 (2024-06-13)

### Other

* Use wildcards in iqtree.

## 0.0.56 (2024-06-13)

### Other

* Add iqtree.

## 0.0.55 (2024-04-19)

### Other

* Catch extract stats.

## 0.0.54 (2024-04-12)

### Other

* Mark captus directories temp.

## 0.0.53 (2024-04-12)

### Other

* Catch more captus output.

## 0.0.52 (2024-04-05)

### Other

* Captus typo:

## 0.0.51 (2024-04-05)

### Other

* Working captus v1.

* Assemble works but extract fails.

## 0.0.50 (2024-03-28)

### Other

* Shadow for purge.

## 0.0.49 (2024-03-27)

### Other

* Use shadow directories to avoid CD.

## 0.0.48 (2024-03-21)

### Other

* Constrict sample names.

## 0.0.47 (2024-03-15)

### Other

* Tag.

* Basic purge haplotigs function.

## 0.0.46 (2024-03-08)

### Other

* Rm typo.

## 0.0.45 (2024-03-08)

### Other

* Clean up RM.

## 0.0.44 (2024-03-01)

### Other

* Use module.

* Module to get reads from hybpiper output.

* Use read_usage_target to extract the reads used in assembly.

## 0.0.43 (2024-02-29)

### Other

* Hybpiper now reads samples from a namelist.

## 0.0.42 (2024-02-23)

### Other

* Bpdownload: cache csv and xml file to avoid repeated parsing.

* Skip SE SRR files.

* Test bpdownload module.

* Cache csv files in bpdownload.

## 0.0.41 (2024-02-07)

### Other

* Braker3 pick up gtf and gff.

## 0.0.40 (2024-02-01)

### Other

* Add header length parameter.

## 0.0.39 (2024-02-01)

### Other

* Tag.

## 0.0.38 (2024-02-01)

### Other

* Split bamfiles.

* Fa readme.

## 0.0.37 (2024-02-01)

### Other

* Tag.

* Add optional parameters.

* Update readme.

## 0.0.36 (2024-01-31)

### Other

* Tag.

* Remove comments.

* Conditional interproscan.

* Add note about interproscan'

* Tidy up env vars.

* Move build script to container recipe repo.

## 0.0.35 (2024-01-24)

### Other

* Braker version bump.

* Wrong path.

## 0.0.34 (2024-01-19)

### Other

* Tag.

* Merge pull request #1 from TomHarrop/funannotate.

  Funannotate

* Merge branch 'main' into funannotate.

* Merge but add interproscan warning.

* Update vars in singularity build script.

* Add SLURM build.

* Working annotation but without interproscan.

* Download ipr database.

* Working up to interproscan"

* Draft interproscan rule.

* Gitignlore.

* Running predict and emapper.

* Working up to train.

* Funannotate db test.

* Work out db extraction / download.

## 0.0.33 (2024-01-11)

### Other

* Turn off strandedness.

## 0.0.32 (2024-01-11)

### Other

* Hisat2.

## 0.0.31 (2023-12-29)

### Other

* Wrong flag in samtools.

## 0.0.30 (2023-12-27)

### Other

* Write index to outdir.

## 0.0.29 (2023-12-21)

### Other

* Relax sample names.

## 0.0.28 (2023-12-21)

### Other

* Allow running without annothation.

## 0.0.27 (2023-12-21)

### Other

* Removed drafts.

* Draft files.

* Imports.

* Forgot container.

## 0.0.26 (2023-12-14)

### Other

* Forgot container.

## 0.0.25 (2023-12-14)

### Other

* No need to combine samples in bbmap.

## 0.0.24 (2023-12-14)

### Other

* Add bbmap.

* Star pipeline.

## 0.0.23 (2023-12-12)

### Other

* Check the bamfile before running braker3.

## 0.0.22 (2023-12-01)

### Other

* Biocontainer.

## 0.0.21 (2023-12-01)

### Other

* Add resources for paragone.

## 0.0.20 (2023-11-30)

### Other

* Merge branch 'paragone'

* Use module.

* Simplify paragone module.

* Too complicated, simplify.

* Make internal and external outgroups optional.

* Targetfile note.

* Working tutorial.

* Qc_trees_and_extract_fasta.

* Initial paragone step.

## 0.0.19 (2023-11-30)

### Other

* Exponent time for repair.

## 0.0.18 (2023-10-27)

### Other

* Add some logging info for URLs.

* Add logging info for URL checking.

* Don't need to configure.

* Configure ncbi.

* Initial download working from accession alone.

## 0.0.17 (2023-09-28)

### Other

* Be stricter on targetfiles.

## 0.0.16 (2023-09-28)

### Other

* Collect targetfile in hybpiper.

## 0.0.15 (2023-09-27)

### Other

* Less ram usage for HP.

## 0.0.14 (2023-09-21)

### Other

* Deal with whitespace in proteins and genome.

## 0.0.13 (2023-09-21)

### Other

* Don't try to reformat proteins.

## 0.0.12 (2023-09-21)

### Other

* Simplify output.

## 0.0.11 (2023-09-21)

### Other

* Wrong arg.

## 0.0.10 (2023-09-21)

### Other

* Wrong arg.

## 0.0.9 (2023-09-21)

### Other

* Logger.

## 0.0.8 (2023-09-21)

### Other

* Make braker evidence optional.

* Specify braker output.

## 0.0.7 (2023-09-19)

### Other

* Why re.

* Add braker.

* Remove bamfile.

* Test commit.

* Add braker.

## 0.0.6 (2023-09-18)

### Other

* Remote module.

* Repeatmaskter.

## 0.0.5 (2023-09-14)

### Other

* Retrieve paralogs.

* Add paralogs.

## 0.0.4 (2023-09-14)

### Other

* Don't search for files, name them appropriately.

## 0.0.3 (2023-09-14)

### Other

* Tag.

* Hook into snakemake logger.

* Handle missing files.

## 0.0.2 (2023-09-13)

### Other

* Add targetfile check to hybpiper.

* Check steps.

* Add check targetfile.

## 0.0.1 (2023-09-13)

### Other

* Try from GH.

* Working hybpiper module.

* Working hybpiper?

* Initial commit.

* Initial commit.
