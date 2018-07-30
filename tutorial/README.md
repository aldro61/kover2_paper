# Tutorial: genotype-to-phenotype models with Kover

This tutorial will walk you through an application of Kover to a set of genomes labelled according to their phenotypes. We use data from the paper (see [data](../data/)), but you could use your own. Specifically, we will learn a model that predicts azithromycin resistance in *Neisseria gonorrhoeae*.

## Outline

* [Getting the example data](#getting-the-example-data)
* [Creating a Kover dataset](#creating-a-kover-dataset)
* [Splitting the data into training and testing sets](#splitting-the-data-into-training-and-testing-sets)
* [Learning models with:](#learning-models)
  * [Set Covering Machines](#set-covering-machines)
  * [Classification and Regression Trees](#classification-and-regression-trees)
* [Interpreting the learned models](#interpreting-the-learned-models)
  * [Annotating k-mers](#annotating-k-mers)
  * [Analyzing equivalent rules](#analyzing-equivalent-rules)


## Getting the example data

First, download the [example data](http://graal.ift.ulaval.ca/adrouin/kover-example-data.zip) (~? Mb), which contains the genome of 392 *Neisseria gonorrhoeae* isolates, along with their susceptibility to azithromycin.
## Creating a Kover dataset

Before using Kover to learn a model, we must package the genomic and phenotypic data into a [Kover dataset](doc_dataset.html#creating-a-dataset), which relies on the HDF5 library to store a compressed representation of the data ([details here](https://github.com/aldro61/kover/wiki/Kover-Dataset-Format)).

To create a dataset, use the following command:

```
kover dataset create from-contigs --genomic-data genome_contigs.tsv --phenotype-description "Azithromycin resistance" --phenotype-metadata metadata.tsv --output example.kover --temp-dir contigs/tmp --progress
```

This produces a dataset file called "example.kover". From now on, you no longer need the original data files.

You can use the [kover dataset info](doc_dataset.html#listing-information-about-a-dataset) command to print information about the dataset. For example, t print the number of genomes and k-mers in the dataset, use:

```
kover dataset info --dataset example.kover --genome-count --kmer-count
```

Your dataset contains **392 genomes vs 4 766 702 k-mers**! This is know as the *fat data* setting, which is very different from the *big data* setting in which the number of examples (genomes) is greater than the number of features (k-mers).

## Splitting the data into training and testing sets

In order to measure the accuracy of the model obtained using Kover, we must split the dataset into a training set and a 
validation set. The training set will be used to learn a model and the validation set will be used to estimate its accuracy.
A Kover dataset can contain multiple splits of the data. The command used for splitting a dataset is [kover dataset split](doc_dataset.html#splitting-a-dataset).

Kover implements a machine learning algorithm and thus has [hyperparameters](doc_learning.html#understanding-the-hyperparameters),
which are free parameters that must be tuned to the data at hand. The most widely used method for setting hyperparameter values
is [k-fold cross-validation](doc_learning.html#k-fold-cross-validation).
In this example, we will use 10-fold cross-validation.

The following command creates a split of the data called "example_split", which uses 80% of the genomes for training and
20% for testing. It also creates 10 cross-validation folds. The data are partitioned randomly, using 2 as the random seed.

```
kover dataset split --dataset example.kover --id example_split --train-size 0.80 --folds 10 --random-seed 2 --progress
```

## Learning models

### Set Covering Machines

### Classification and Regression Trees
** Dont forget to mention multiclass can be done and how **

## Interpreting the learned models

### Annotating k-mers
Show a simple example using nucleotide BLAST

### Analyzing equivalent rules
Show a simple example using UGENE
