# PSASS-snake

PSASS-snake is a pipeline to search for sex signal in pooled-sequencing data. It is meant to be a more reproducible, portable, and maintainable update to [PSASS-process](https://github.com/RomainFeron/PSASS-process). This pipeline was developed for Linux only, and was specifically designed and tested to work on HPCs (clusters).

PSASS-snake is implemented with Snakemake and uses Conda to handle software dependencies (although the pipeline can use local version of the software if necessary).

Although the PSASS suite was developed specifically to study sex determination, it can be used on any pooled-sequencing data with some parameter adjustments.

## Setup

### Requirements

PSASS-snake requires Conda and Snakemake to be installed and setup. If you do not already have them installed, refer to the [Installing Conda and Snakemake](#installing-conda-and-snakemake) section below.

### Clone the PSASS-snake repository

Because PSASS-snake is a snakemake pipeline, you should clone the repository for each new analysis. This way, the pipeline and the data are in the same repository, ensuring that the analyses are entirely reproducible.

**HTTPS:**

```bash
git clone https://github.com/RomainFeron/PSASS-snake.git
```

**SSH (recommended, but requires setup):**

```bash
git clone git@github.com:RomainFeron/PSASS-snake.git
```

## Running PSASS-snake

### General principle

For a snakemake pipeline like PSASS-vis, a directory contains both the pipeline's implementation and the data used for a single run; this setup ensures reproducibility of the analyses.

To run PSASS, the first step is to clone the [GitHub repository](https://github.com/RomainFeron/Psass-snake) or copy another run's directory. Then, the data required for this run should be placed in the directory. The recommended organization is a `data/` folder at the root of the directory, with a `genome` subdirectory containing the assembly fasta file and a `reads` subdirectory containing all the reads files within the `data` firectory.

Then, the run can be configured by editing the file `config.yaml` with the desired parameter values. To execute the pipeline, you can then simply run `snakemake` from within the main directory.

### Editing the config file

### Running the pipeline

### Special instructions to run PSASS-snake on HPCs (clusters)


## Installing Conda and Snakemake

If this is the first time you are using Conda and/or snakemake, follow these instructions:

### Install Conda

Full instructions available [here](https://conda.io/projects/conda/en/latest/user-guide/install/linux.html).

**Download the installer script and run it:**

```bash
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

Follow instructions from the prompt.

**Restart your shell:**

```bash
source ~/.bashrc
```

**Note:** If you're not using bash as your shell, source the appropriate file (*e.g.* `~/.zshrc` if you're using zsh)

**Initialize conda for your shell:**

```bash
conda init bash
# Replace bash with whatever shell you are using
```

**Update conda:**

The Conda version from the official installer is not always the latest update. To make sure Conda is up-to-date, run:

```bash
conda update conda
```

### Install Snakemake

The recommended way to install and run Snakemake is to create a conda environment specifically for it:

```bash
# Create a new empty environment called "snakemake"
conda create --name snakemake
# Activate the environment "snakemake"
conda activate snakemake
# Install snakemake from the Bioconda channel (conda-forge contains dependencies)
conda install -c conda-forge -c bioconda snakemake
```

You can now activate the environment snakemake and run snakemake from it. It is advised to keep the environment as clean as possible, *i.e.* only install software related to running snakemake.
