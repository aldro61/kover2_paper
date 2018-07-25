# Data

Use the following steps to reproduce any dataset that was used in the paper.

## Step 1: Choose a dataset

The datasets are organized by species and antibiotic. Each directory contains a metadata file with two columns: 1) the genome identifiers (PATRIC) and 2) the phenotypes.

## Step 2: Download the genomes

To download all the genomes in a dataset, simply feed the corresponding `metadata.tsv` file to the `download_genomes.py` script. For example:

```
python download_genomes 'mycobacterium tuberculosis/kanamycin/metadata.tsv' genomes_dir
```

Now `genomes_dir` contains the genomes and a summary file called `genome_paths.tsv`, containing the path to each genome, has been created in the current directory. You will need this file in the next step.

## Step 3: Create a Kover dataset

Follow this [tutorial](TODO), but replace the label file by the `metadata.tsv` that you just used and the genome file by the corresponding `genome_paths.tsv`.

**To obtain the exact same datasets as in the paper, use 80% of the data for training and create 10 splits with random seeds 0, 1, 2, 3, 4, 5, 6, 7, 8, 9.**
