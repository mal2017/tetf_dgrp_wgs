# TETF_DGRP_WGS: processing of WGS data for the Ellison lab TETF project

## dependencies

- `snakemake` (tested with v7.21.0 on Ubuntu 22 LTS and 7.18.2 on CentOS 7)
- `mamba` (tested with version 1.0.0 on Ubuntu 22 LTS and 0.9.1 on CentOS 7). `conda` alone should work just fine.
- `singularity` (tested with v3.1.0-1 on CentOS 7) or `apptainer` (tested with v1.1.3 on Ubuntu 22 LTS)

## usage

```
snakemake --use-conda --use-singularity --cores <threads>
```

A test run (just DGRP_100 and DGRP_101) can be performed by adding `--config RUN_TYPE=test`

## Description

This workflow performs basic DNAseq alignment for a large dataset of WGS on the DGRP.

In its current iteration, the primary output is `results/copies/copies.tsv`.

```
Strain	| sample_name	sequence	length	bases	median.cov	est.copies
DGRP_313	DGRP_313_SAMN00014255	R1A1-element	5356	852468	159.16	64.96326530612244
DGRP_313	DGRP_313_SAMN00014255	roo	9092	1420089	156.19	63.75102040816326
...
...
```

