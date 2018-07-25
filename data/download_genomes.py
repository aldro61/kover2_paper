"""
Download the genomes in a dataset from the PATRIC database

Note: depends on PATRIC tools (https://github.com/aldro61/patric_tools)

"""
import argparse

from os.path import abspath, join
from patric_tools.genomes import download_genome_contigs


def download_genome(g_id, outdir):
    print "... Genome: {0!s} -> {1!s}/{0!s}.fna".format(g_id, outdir)
    download_genome_contigs(g_id, outdir=outdir)


parser = argparse.ArgumentParser(description="Download genomes for a dataset")
parser.add_argument('metadata', type=str, help="The path to the metadata.tsv file for the dataset.")
parser.add_argument('outdir', type=str, help="The directory in which to output the genomes")
args = parser.parse_args()

genome_ids, labels = zip(*[l.strip().split("\t") for l in open(args.metadata, "r")])

print "Downloading genomes"
genome_paths = []
for g_id in genome_ids:
    download_genome(g_id, outdir=args.outdir)
    genome_path = join(abspath(args.outdir), g_id + ".fna")
    genome_paths.append(g_id + "\t" + genome_path)

summary_path = abspath("genome_paths.tsv")
print "Writing summary TSV file:", summary_path
open(summary_path, "w").write("\n".join(genome_paths))
