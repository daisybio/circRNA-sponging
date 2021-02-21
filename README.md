# ![nf-core/circrnasponging](docs/images/nf-core-circrnasponging_logo.png)

**Analysis of circRNA and miRNA sponging**.

<!--[![GitHub Actions CI Status](https://github.com/nf-core/circrnasponging/workflows/nf-core%20CI/badge.svg)](https://github.com/nf-core/circrnasponging/actions)
[![GitHub Actions Linting Status](https://github.com/nf-core/circrnasponging/workflows/nf-core%20linting/badge.svg)](https://github.com/nf-core/circrnasponging/actions)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A520.04.0-brightgreen.svg)](https://www.nextflow.io/)

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](https://bioconda.github.io/)
[![Docker](https://img.shields.io/docker/automated/nfcore/circrnasponging.svg)](https://hub.docker.com/r/nfcore/circrnasponging)
[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23circrnasponging-4A154B?logo=slack)](https://nfcore.slack.com/channels/circrnasponging)-->

This pipeline was written by Octavia Ciora as part of her Advanced Lab Course Bioinformatics under the supervision of Dr. Markus List.

## Table of Contents

* [Introduction](#introduction)
* [Pipeline Summary](#pipeline-summary)
* [Documentation](#documentation)
  + [Input](#input)
  + [Configuration](#configuration)
  + [Output](#output)
* [Running the Pipeline](#running-the-pipeline)

## Introduction

**nf-core/circrnasponging** is a pipeline for the systematical analysis of circRNAs that act as sponges for miRNAs. It requires total RNA and small RNA sequencing data and is based on the hypothesis that the negatively correlating expression of miRNAs and circRNAs (having miRNA binding sites) is an indicator of sponging.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker containers making installation trivial and results highly reproducible.

## Pipeline Summary

By default, the pipeline currently performs the following:
<!-- TODO nf-core: Fill in short bullet-pointed list of default steps of pipeline -->
* circRNA analysis:
  + totalRNA read mapping (using `STAR`)
  + circRNA detection and quantification (using `circExplorer2`)
  + aggregation of detected circRNAs over all samples and filtering
  + extraction of fasta sequences for filtered circRNAs
* miRNA analysis:
  + smRNA read mapping and quantification (using `miRDeep2`)
  + aggregation of identified miRNAs over all samples and filtering
* Binding sites:
  + miRNA binding sites detection on quantified circRNAs (`miranda`)
  + binding sites processing and filtering
* Sponging analysis between circRNAs and miRNAs
* Overall results summaries and plots

## Documentation
<!--The nf-core/circrnasponging pipeline comes with documentation about the pipeline: [usage](https://nf-co.re/circrnasponging/usage) and [output](https://nf-co.re/circrnasponging/output). -->
This analysis is based on both totalRNA (or rRNA-depleted) and smRNA data coming from the same samples. It is recommended to run this analysis with a minimum of 5 samples. The exact input format and how to get the needed reference files is described below.

### Input
These pipeline can be applied in two ways:
* Option 1: Input fastq files only for totalRNA data and tabulated raw read counts for miRNAs quantified previously. It is recommended to use the [`nf-core/smarnaseq`](https://github.com/nf-core/smrnaseq) pipeline for the quantification of miRNAs.
* Option 2: Input fastq files for both totalRNA and smRNA data. In this case the pipeline perform the read mapping and quantification of circRNAs and miRNAs.

#### Samplesheet
The pipeline requires a tab-separated samplesheet file with information about the samples and the paths to the fastq.gz read files.

##### single-end
If your totalRNA sequencing data is single-end, the samplesheet file should have the following structure:
* Column 1: sample name
* Column 2: path to totalRNA fastq.gz file corresponding to the sample mentioned in column 1
* Column 3: path to the miRNA fastq.gz file corresponding to the sample mentioned in column 1. In case you are running the pipeline without miRNA fastq files (see Option 1 in [Input](#input)), replace the path with NA.
```
      sample     |             totalRNA1               |             smallRNA
---------------- | ------------------------------------|------------------------------------
 cerebellum_rep1 | path/to/<totalRNA_sample1>.fastq.gz | path/to/<smallRNA_sample1>.fastq.gz
 cerebellum_rep2 | path/to/<totalRNA_sample2>.fastq.gz | path/to/<smallRNA_sample2>.fastq.gz
hippocampus_rep1 | path/to/<totalRNA_sample3>.fastq.gz | path/to/<smallRNA_sample3>.fastq.gz
      ...        |                ...                  |               ...
```

##### paired-end
If your totalRNA sequencing data is single-end, the ```dataset.tsv``` file should look like this:
Column 1: sample name
Column 2: path to totalRNA read1 fastq.gz file corresponding to the sample mentioned in column 1
Column 3: path to totalRNA read2 fastq.gz file corresponding to the sample mentioned in column 1
Column 4: path to the miRNA fastq.gz file corresponding to the sample mentioned in column 1. In case you are running the pipeline without miRNA fastq files (see Option 1 in [Input](#input)), replace the path with NA.

```
      sample     |                totalRNA1               |               totalRNA2               |             smallRNA
---------------- | ---------------------------------------|---------------------------------------|-------------------------------------
 cerebellum_rep1 | path/to/<totalRNA_sample1_R1>.fastq.gz | path/to/<totalRNA_sample1_R2>.fastq.gz| path/to/<smallRNA_sample1>.fastq.gz
 cerebellum_rep2 | path/to/<totalRNA_sample2_R1>.fastq.gz | path/to/<totalRNA_sample2_R2>.fastq.gz| path/to/<smallRNA_sample2>.fastq.gz
hippocampus_rep1 | path/to/<totalRNA_sample3_R1>.fastq.gz | path/to/<totalRNA_sample3_R2>.fastq.gz| path/to/<smallRNA_sample3>.fastq.gz
      ...        |                  ...                   |                  ...                  |                ...
```

#### Parameters
For instructions about how to specify the parameters, see [Configuration](#configuration).



### Reference Files
#### Mandatory
Reference files needed for this analysis are:
* GTF file
* genome fasta file
* miR-Base reference files for hairpin and mature miRNA sequences
* gene-pred file 
#### Optional
A STAR index is required for the detection of circRNAs and a bowtie index is needed for the read mapping of smRNAs. If you skip the miRNA quantification step, and input the tabulated miRNA read counts, the bowtie index is not needed. Specifying already built STAR index or bowtie index is supported, but this is optional, since the pipeline generates them once from the fasta file, if they are not provided.

### Configuration
### Output

## Running the Pipeline
After preparing the input files and setting up the configuration file, you can run the analysis using the following command: 
 ```bash
 nextflow run nf-core-circrnasponging/ -c my.config -profile cluster
 ```
 Add the `-resume` tag to continue a previously stopped pipeline execution without losing already finished tasks.
 
 
 
 
 
 
 
 
 
 
<!--
1. Install [`nextflow`](https://nf-co.re/usage/installation)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) or [`Podman`](https://podman.io/) for full pipeline reproducibility _(please only use [`Conda`](https://conda.io/miniconda.html) as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))_

3. Download the pipeline and test it on a minimal dataset with a single command:

    ```bash
    nextflow run nf-core/circrnasponging -profile test,<docker/singularity/podman/conda/institute>
    ```

    > Please check [nf-core/configs](https://github.com/nf-core/configs#documentation) to see if a custom config file to run nf-core pipelines already exists for your Institute. If so, you can simply use `-profile <institute>` in your command. This will enable either `docker` or `singularity` and set the appropriate execution settings for your local compute environment.

4. Start running your own analysis!

    TODO nf-core: Update the example "typical command" below used to run the pipeline

    ```bash
    nextflow run nf-core/circrnasponging -profile <docker/singularity/podman/conda/institute> --input '*_R{1,2}.fastq.gz' --genome GRCh37
    ```

See [usage docs](https://nf-co.re/circrnasponging/usage) for all of the available options when running the pipeline.
-->

<!--## Credits

First scripts were originally written by Octavia Ciora as part of her Bachelor's Thesis at the Chair of Experimental Bioinformatics under the supervision of Prof. Dr. Jan Baumbach and Dr. MArkus List
The pipeline was re-written in Nextflow 
We thank the following people for their extensive assistance in the development
of this pipeline: Dr. Markus List
TODO nf-core: If applicable, make list of people who have also contributed
-->

<!--## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#circrnasponging` channel](https://nfcore.slack.com/channels/circrnasponging) (you can join with [this invite](https://nf-co.re/join/slack)).-->

<!--## Citations

 TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi. -->
<!-- If you use  nf-core/circrnasponging for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!--You can cite the `nf-core` publication as follows: -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->
