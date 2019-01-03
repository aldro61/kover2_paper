# Data

## Reproducing the datasets
Use the following steps to reproduce any dataset that was used in the paper.

### Step 1: Choose a dataset

The datasets are organized by species and antibiotic. Each directory contains a metadata file with two columns: 1) the genome identifiers (PATRIC) and 2) the phenotypes.

### Step 2: Download the genomes

To download all the genomes in a dataset, simply feed the corresponding `metadata.tsv` file to the `download_genomes.py` script. For example:

```
python download_genomes 'mycobacterium tuberculosis/kanamycin/metadata.tsv' genomes_dir
```

Now `genomes_dir` contains the genomes and a summary file called `genome_paths.tsv`, containing the path to each genome, has been created in the current directory. You will need this file in the next step.

### Step 3: Create a Kover dataset

Follow this [tutorial](https://aldro61.github.io/kover/doc_tut_data.html), but replace the label file by the `metadata.tsv` that you just used and the genome file by the corresponding `genome_paths.tsv`.

**To obtain the exact same datasets as in the paper, create 10 splits with 80% of the data for training, 10 folds, and random seeds 0, 1, 2, 3, 4, 5, 6, 7, 8, 9.**

## Source

These datasets were extracted from the PATRIC database (https://www.patricbrc.org/). See our paper for the data acquisition protocol. More info on this database can be found in:

Antonopoulos, D. A., Assaf, R., Aziz, R. K., Brettin, T., Bun, C., Conrad, N., ... & Kenyon, R. W. (2017). PATRIC as a unique resource for studying antimicrobial resistance. Briefings in bioinformatics.

Wattam, A. R., Abraham, D., Dalay, O., Disz, T. L., Driscoll, T., Gabbard, J. L., ... & Machi, D. (2013). PATRIC, the bacterial bioinformatics database and analysis resource. Nucleic acids research, 42(D1), D581-D591.
