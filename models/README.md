# Antimicrobial resistance prediction models

Rule-based models are available for the following algorithms:
* Set Covering Machines (bound selection)
* Classification and Regression Trees (bound selection)

The models are organized by algorithm, species, and antimicrobial agent.


## Visualization

A visual representation of each model is given (see [example](https://github.com/aldro61/kover2_paper/tree/master/models/cart_b/mycobacterium%20tuberculosis/pyrazinamide/repeat_5)). You can generate this representation for your models using the following script:

```
python plot_model.py model.fasta
```


## Format

All models are provided in FASTA format to facilitate their annotation using tools such as [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome). Equivalent rules are also given in separate FASTA files.
