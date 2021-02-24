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
  + [Basic Workflow](#basic-workflow)
  + [Input Files](#input-files)
  + [Additional Features and Advanced Options](#additional-features-and-advanced-options)
  + [Output](#output)
* [Pipeline Execution and Configuration](#pipeline-execution-and-configuration)

## Introduction

**nf-core/circrnasponging** is a pipeline for the systematical analysis of circRNAs that act as sponges for miRNAs. It requires total RNA and small RNA sequencing data and is based on the hypothesis that the negatively correlating expression of miRNAs and circRNAs (having miRNA binding sites) is an indicator of sponging.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker containers making installation trivial and results highly reproducible.

## Pipeline Summary

By default, the pipeline currently performs the following:
<!-- TODO nf-core: Fill in short bullet-pointed list of default steps of pipeline -->
* circRNA analysis:
  + totalRNA read mapping (using `STAR`)
  + circRNA detection and quantification (using `circExplorer2`)
  + aggregation of detected circRNAs over all samples, normalization and filtering
  + extraction of fasta sequences for filtered circRNAs
* miRNA analysis:
  + smRNA read mapping and quantification (using `miRDeep2`)
  + aggregation of identified miRNAs over all samples, normalization and filtering
* Binding sites:
  + miRNA binding sites detection on quantified circRNAs (`miranda`)
  + binding sites processing and filtering
* Correlation analysis between circRNAs and miRNAs expression
* Overall results summaries and plots

## Documentation
<!--The nf-core/circrnasponging pipeline comes with documentation about the pipeline: [usage](https://nf-co.re/circrnasponging/usage) and [output](https://nf-co.re/circrnasponging/output). -->
This analysis is based on both totalRNA (or rRNA-depleted) and smRNA data coming from the same samples. It is recommended to run this analysis with a minimum of 5 samples. The exact input format and how to get the needed reference files is described below.
In order to run the sponging analysis on a dataset using our pipeline, the data has to meet strict requirements.  At least 5 samples from the same organism are needed, because computing the correlation between 4 or less samples is unreasonable. For each sample, both total RNA (or rRNA-depleted) and small RNA data should be available.  Further instructions regarding the exact input format, reference files and configuration are explained below.

### Basic Workflow
The following options are mandatory for executing the basic workflow of the pipeline:


```
BASIC OPTIONS:
  --samplesheet path/to/sampleseet.tsv
  --out_dir path/to/output_directory
  --species species_in_3_letter_code # hsa for human, mmu for mouse etc.
  --miRNA_adapter adapter_sequence # miRNA adapter used for trimming
    
REFERENCE FILES
  --fasta path/to/genome.fasta
  --gtf path/to/gtf_file
  --gene_pred path/to/gene_annotation
  --mature_fasta path/to/mature_fasta
  --mature_other_fasta path/to/mature_other_fasta
  --hairpin_fasta path/to/hairpin_fasta
}
```

### Input Files
#### Samplesheet
The pipeline requires a tab-separated samplesheet file containing the sample names and the paths to the corresponding read files in fastq.gz format. By default, the  totalRNA sequencing data is considered to be paired-end. If the total RNA data is single-end, the instructions below should be followed. The samplesheet should be a tab-separated file following the structure:

```
   sample  |                totalRNA1               |               totalRNA2               |             smallRNA
-----------|----------------------------------------|---------------------------------------|-------------------------------------
  sample1  | path/to/<totalRNA_sample1_R1>.fastq.gz | path/to/<totalRNA_sample1_R2>.fastq.gz| path/to/<smallRNA_sample1>.fastq.gz
  sample2  | path/to/<totalRNA_sample2_R1>.fastq.gz | path/to/<totalRNA_sample2_R2>.fastq.gz| path/to/<smallRNA_sample2>.fastq.gz
  sample3  | path/to/<totalRNA_sample3_R1>.fastq.gz | path/to/<totalRNA_sample3_R2>.fastq.gz| path/to/<smallRNA_sample3>.fastq.gz
   ...     |                  ...                   |                  ...                  |                ...
```

#### Reference Files
Instructions for generating the miRNA reference files can be found in the [`miRDeep2 tutorial`](https://drmirdeep.github.io/mirdeep2_tutorial.html). Required are the mature and hairpin miRNA sequences for the organism used in the analysis, e.g. mouse(mmu). Additionally, mature miRNA sequences from related species are needed, e.g. rat(rno) and human (hsa).

### Additional Features and Advanced Options
#### Skip miRNA Quantification
The pipeline offers the option to skip the miRNA quantification step, if this has been done in advance and the raw read counts are available. In this case, the pipeline performs the read mapping and quantification only for circRNAs. It is recommended to use [`nf-core/smarnaseq`](https://github.com/nf-core/smrnaseq) pipeline for the quantification of miRNAs. The tabulated raw read counts can be passed directly to the pipeline using the following parameter:

```
  --miRNA_raw_counts path/to/miRNA_raw_counts.tsv
```
The file should be tab-separated and follow the structure shown below. The header should contain the same sample names as defined in the samplesheet. The last column called `smallRNA` in the samplesheet should have the value ’NA’ among all samples.  When running the pipeline in this scenario, the following parameters are no longer required:
```
  --species, --miRNA_adapter, --mature_other_fasta, --hairpin_fasta
```
This additional feature gives the user the freedom to choose the miRNA quantificationprocedure of his preference, while also saving time, in case the quantified data is alreadyavailable.

```
     miRNA    |sample1|sample2|   ...
--------------|-------|-------|---------
mmu-let-7a-5p |   53  |   0   |   ...
mmu-let-7b-3p |   37  |   93  |   ...
      ...     |  ...  |  ...  |   ...
```
#### Single-end total RNA data
By default, the pipeline works with paired-end total RNA sequencing data.  Nevertheless,it also supports single-end data, in which case the following additional parameter shouldbe used:

```
  --single_end 
      default: false
        true    total RNA sequencing data is single-end
	false   total RNA sequencing data is paired-end
```
In addition, the samplesheet should be adapted to match the following structure:
```
   sample  |             totalRNA1               |             smallRNA
-----------|-------------------------------------|------------------------------------
  sample1  | path/to/<totalRNA_sample1>.fastq.gz | path/to/<smallRNA_sample1>.fastq.gz
  sample2  | path/to/<totalRNA_sample2>.fastq.gz | path/to/<smallRNA_sample2>.fastq.gz
  sample3  | path/to/<totalRNA_sample3>.fastq.gz | path/to/<smallRNA_sample3>.fastq.gz
    ...    |                ...                  |               ...
```

#### Sample Grouping for Better Visualization
This feature can be used if the samples can be partitioned into different conditions or groups. If the sample grouping is specified, the samples will be colored according to their grouping in the plots showing sponging candidates, facilitating a better visualization andresults interpretation. The file containing the sample grouping should be structured as shown below and specified using the parameter:
```
  --sample_group path/to/sample.group.tsv
```
```
  sample | group   
---------|-------
 sample1 | group1
 sample2 | group1
 sample3 | group2
 sample4 | group3
   ...   |  ...  
```
#### Advanced Filtering
After normalizing the raw read counts for both circRNAs and miRNAs, the pipelinefilters out entries which have a low expression level. By default, only entries having at least 5 reads in at least 20% of samples are used in the downstream analysis. Stricter filtering conditions might lead to better identification of sponging candidates, as long as the parameters are not too high, in which case interesting candidates might be missed. The parameter values can be changed with the options:
```
  --read_threshold
      default: 5
	real >= 0    read counts under this threshold are considered to be low expressed
  --sample_percentage
      default: 0.2
        0<= real <=1    minimum percentage of samples that should haveno low expression
```
### Output
The output folder is structured as shown below. The circRNA/miRNA results for each sample are stored in folder `samples`. The tabulated circRNA and miRNA counts summarized over all samples are `results/circRNA` and `results/miRNA`. The results of the sponging analysis are stored in the subfolder `results/sponging`. 

```
├─── output_folder
│   ├─── samples
|   │   ├─── sample_1
|   |   |   |─── circRNA_detection
|   |   |   └─── miRNA_detection
|   │   ├─── sample_2
|   |   |   |─── circRNA_detection
|   |   |   └─── miRNA_detection
|   │   ├─── ...
│   ├─── results
|   |   |─── circRNA
|   |   |   |─── circRNA_counts_raw.tsv
|   |   |   └─── circRNA_counts_filtered.tsv
|   |   |─── miRNA
|   |   |   |─── miRNA_counts_raw.tsv
|   |   |   └─── miRNA_counts_filtered.tsv
|   |   |─── binding_sites
|   |   |─── sponging
|   |   |   |─── sponging_statistics.txt
|   |   |   |─── filtered_circRNA_miRNA_correlation.tsv
|   |   |   └─── plots
└── └── └── 
```

## Pipeline Execution and Configuration
To execute the pipeline, run the following command:
```nextflow run nf-core-circrnasponging/ -c my.config -profile cluster -resume```
where `my.config` is a configuration file specifying parameters and execution settings. An example configuration file is shown below:
```
params {
        samplesheet = "path/to/sampleseet.tsv"
        out_dir = "path/to/out_dir"
        species = "mmu"
        miRNA_adapter = "TGGAATTCTCGGGTGCCAAGG"
        single_end = true
        fasta = "path/to/genome.fasta"
        gtf = "path/to/genome.gtf"
        gene_pred = "path/to/genome.annot"
        mature_fasta = "path/to/mature_mmu.fa"
        mature_other_fasta = "path/to/mature_rno_hsa.fa"
        hairpin_fasta = "path/to/hairpin_mmu.fa"
    }
    profiles {
         standard {
             process.executor = 'local'
         }
         cluster {
             executor.queueSize = 20
             process.executor = 'slurm'
             process.cpu = '8'
    	 process.memory = '50 GB'
         }
    }
```





 
 
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
We thank the following people for their extensive assistance in the development
of this pipeline:
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
