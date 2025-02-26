![Build Docker (env/Dockerfile)](https://github.com/AndersenLab/isotype-nf/workflows/Build%20Docker%20(env/align.Dockerfile)/badge.svg)


# isotype-nf

The [isotype-nf](https://github.com/AndersenLab/isotype-nf) pipeline performs isotype group calls for wild isolate variant data __at the strain level__, and outputs isotype assignments and related information. Those isotypes can be used for downstream analysis including isotype reference variant calling, [wi-gatk-nf (variant calling)](http://andersenlab.org/dry-guide/pipelines/pipeline-wi/) and other analyses.

This page details how to run the pipeline. You can also find more information on the Andersen Lab [dry guide](http://andersenlab.org/dry-guide/latest/pipelines/pipeline-isotype/).


# Pipeline overview

```
_____   ______     ___    _________  ____  ____  _______  ________       ____  _____  ________  
|_   _|.' ____ \  .'   `. |  _   _  ||_  _||_  _||_   __ \|_   __  |     |_   \|_   _||_   __  | 
  | |  | (___ \_|/  .-.  \|_/ | | \_|  \ \  / /    | |__) | | |_ \_|______ |   \ | |    | |_ \_| 
  | |   _.____`. | |   | |    | |       \ \/ /     |  ___/  |  _| _|______|| |\ \| |    |  _|    
 _| |_ | \____) |\  `-'  /   _| |_      _|  |_    _| |_    _| |__/ |      _| |_\   |_  _| |_     
|_____| \______.' `.___.'   |_____|    |______|  |_____|  |________|     |_____|\____||_____|    
                                                                                               

To run the pipeline:

nextflow main.nf --help
nextflow main.nf --debug
nextflow main.nf --vcf_file=/path/to/vcf_file --species c_elegans -output-dir=/path/to/output
nextflow main.nf --vcf_file=/path/to/vcf_file --bam_location=/path/to/bams --previous_isotypes=/path/to/previous_isotype_file -output-dir=/path/to/output

    parameters                 description                           Set/Default
    ==========                 ===========                           ========================
    --debug                    Use --debug to indicate debug mode    false
    --vcf_file                 All strains VCF file                  null

    --species                  Species to call isotypes from         null
    and / or
    --cutoff                   Concordance cutoff for isotype calls  null
    --bam_location             Directory of BAM files                null
    --previous_isotypes        File containing previous isotypes     null
    
    username                                                         null

	HELP: http://andersenlab.org/dry-guide/pipelines/pipeline-isotype/
```

## Software requirements

* Nextflow v24+ (see the dry guide on Nextflow [here](http://andersenlab.org/dry-guide/rockfish/rf-nextflow/) or the Nextflow documentation [here](https://www.nextflow.io/docs/latest/getstarted.html)). On Rockfish, you can access this version by loading the `nf24_env` conda environment prior to running the pipeline command:

```
module load python/anaconda
source activate /data/eande106/software/conda_envs/nf24_env
```

* Singularity. On Rockfish, you can get this with `module load singularity` before running


# Usage

## Testing on Rockfish

*This command uses a test dataset*

```
nextflow run -latest andersenlab/isotype-nf --debug
```

## Running on Rockfish

You should run this in a screen or tmux session.

*Note: if you are having issues running Nextflow or need reminders, check out the [Nextflow](http://andersenlab.org/dry-guide/rockfish/rf-nextflow/) page.*

```
nextflow main.nf --vcf_file=/path/to/vcf_file --species c_elegans -output-dir=/path/to/output
```
or
```
nextflow main.nf --vcf_file=/path/to/vcf_file --bam_location=/path/to/bams --previous_isotypes=/path/to/previous_isotype_file -output-dir=/path/to/output
```

# Parameters

## -profile

There are three configuration profiles for this pipeline.

* `rockfish` - Used for running on Rockfish (default)

>[!Note]
>If you forget to add a `-profile`, the `rockfish` profile will be chosen as default

## --vcf_file

The `vcf file` for isotype calling is the output from the [wi-gatk](https://github.com/AndersenLab/wi-gatk) pipeline. The `vcf file` **must be gzipped**, is the **full path to the vcf file ** (even if it is in your current directory)

>[!Note]
>Remember that in `--debug` mode the pipeline will use the vcf file located in `test_data/vcf.gz`.

## --debug (optional)

You should use `--debug` for testing/debugging purposes. This will run the debug test set (located in the `test_data` folder) using your specified configuration profile (e.g. rockfish).

For example:

```
nextflow run -latest andersenlab/isotype-nf --debug -resume
```

Using `--debug` will automatically set the vcf file to `test_data/vcf.gz`

### --species (optional if --bam_location, --cutoff, and --previous_isotypes are specified)

Must be "c_elegans", "c_briggsae" or "c_tropicalis" to select automatically select cutoff, bam folder, and previous isotype groups. Manually specifying any of these will override the default path for a given species. If all of them are manually specified, --species will not affect the workflow.

### --bam_location (optional if --species is specified)

Location of bam folder containing coverage records for each strain. Inferred from species if not specified and --species is specified

### --previous_isotypes (optional if --species is specified)

Location of previous isotype group calls file. Inferred from species if not specified and --species is specified. This file should be tab-separated with at least 3 columns: strain, isotype, and isotype_ref_strain.

### --cutoff (optional if --species is specified)

The minimum concordance required for isotype clustering. Inferred from species if not specified and --species is specified

## -output-dir (optional)

__default__ = WI-{today's date} where the date is formatted as `YYYYMMDD` 

A directory in which to output results


# Output

```
├── gtcheck.tsv
├── isotype_groups.tsv
├── wi_isotype_sample_sheets.txt
└── isotype_comparison.pdf

```

Most files should be obvious. A few are detailed below.

* __gtcheck.txt__ - Contains all of the pairwise genotype comparisons of valid genotypes for every strain.
* __isotype_groups.tsv__ - A tab-separated file with strains, their isotype group, and the isotype reference strain.
* __wi_isotype_sample_sheet.txt__ - A sample sheet for rerunning the wi-gatk pipeline for isotype reference strains only.
* __isotype_comparison.pdf__ - A set of plots showing joins and splits in the new isotype groups compared to the previous isotype groups.

# Relevant Docker Images

* `andersenlab/numpy` ([link](https://hub.docker.com/r/andersenlab/numpy)): Docker image is created within this pipeline using GitHub actions. Whenever a change is made to `env/Dockerfile` or `.github/workflows/build_docker.yml` GitHub actions will create a new docker image and push if successful



