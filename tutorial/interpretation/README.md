# Interpreting Kover models

One particularity of the models learned with Kover is that they are highly interpretable. The models make predictions based on rules that capture the presence/absence of k-mers. Below, we show how simple it is to go from model to biological interpretation with these models.

### Annotating k-mers

After learning a model, the results directory contains a file called `model.fasta`, which contains each k-mer in the model, along with an informative header. Conveniently, these FASTA files can be directly inputted into tools such as [Nucleotide BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastSearch).

In this example, we will use the model 

### Analyzing equivalent rules
Show a simple example using UGENE