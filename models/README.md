# Antimicrobial resistance prediction models

Rule-based models are available for the following algorithms:
* Set Covering Machines (model selection: bound)
* Classification and Regression Trees (model selection: bound)

The models are organized by algorithm, species, and antimicrobial agent.


## Format

All models are provided in FASTA format to facilitate their annotation using tools such as [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome). Equivalent rules are also given in separate FASTA files.


## Visualization

A visual representation of each model is given in PDF format. You can generate this representation for your models using the following script:

```
python plot_model.py model.fasta
```
