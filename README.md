# PSASS-workflow

PSASS-workflow is a pipeline to search for sex signal in pooled-sequencing data. It is meant to be a more reproducible, portable, and maintainable update to [PSASS-process](https://github.com/RomainFeron/PSASS-process). This pipeline was developed for Linux only, and was specifically designed and tested to work on HPCs (clusters).

PSASS-workflow is implemented with Snakemake and uses Conda to handle software dependencies (although the pipeline can use local version of the software if necessary).

Although the PSASS suite was developed specifically to study sex determination, it can be used on any pooled-sequencing data with some parameter adjustments.

## Setup

### Requirements

PSASS-workflow requires Conda and Snakemake to be installed and setup. If you do not already have them installed, refer to the [Installing Conda and Snakemake](#installing-conda-and-workflowmake) section below.

### Clone the PSASS-workflow repository

Because PSASS-workflow is a snakemake pipeline, you should clone the repository for each new analysis. This way, the pipeline and the data are in the same repository, ensuring that the analyses are entirely reproducible.

**HTTPS:**

```bash
git clone https://github.com/SexGenomicsToolkit/PSASS-workflow.git
```

**SSH (recommended, but requires setup):**

```bash
git clone git@github.com:SexGenomicsToolkit/PSASS-workflow.git
```

### Installing Conda and Snakemake

If this is the first time you are using Conda and/or snakemake, follow [this guide](https://gist.github.com/RomainFeron/da9df092656dd799885b612fedc9eccd).

## Running PSASS-workflow

### General principle

For a snakemake pipeline like PSASS-vis, a directory contains both the pipeline's implementation and the data used for a single run; this setup ensures reproducibility of the analyses.

To run PSASS, the first step is to clone the [GitHub repository](https://github.com/SexGenomicsToolkit/PSASS-workflow) or copy another run's directory. Then, the data required for this run should be placed in the directory. The recommended organization is a `data/` folder at the root of the directory, with a `genome` subdirectory containing the assembly fasta file and a `reads` subdirectory containing all the reads files within the `data` firectory.

Then, the run can be configured by editing the file `config.yaml` with the desired parameter values. To execute the pipeline, you can then simply run `snakemake` from within the main directory.

### Editing the config file

### Running the pipeline

### Special instructions to run PSASS-workflow on HPCs (clusters)

