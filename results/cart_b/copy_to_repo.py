"""
Export models to the github repository that complements the paper

"""
import json
import numpy as np

from os import listdir, mkdir, system
from os.path import *
from urllib import quote

root = "cart_bound"
destination = "/Users/alexandredrouin/dev/git/kover2_paper/models/cart_b"

for d in listdir(root):
    d = join(root, d)

    if not all([isdir(join(d, "train_0.800_seed_{0:d}_10_folds".format(i))) for i in range(10)]):
        continue

    for sd in listdir(d):
        if "seed" not in sd:
            continue

        sd = join(d, sd)
        results = json.load(open(sd + "/results.json", "r"))
        ds_info = json.load(open(results["data"]["path"].replace("/scratch/adrouin/patric_data", "/Users/alexandredrouin/graham/project/data/patric").replace(
            "/scratch/adrouin", "/Users/alexandredrouin/graham/project").replace("dataset.kover", "dataset_info.json"), "r"))
        species = ds_info["species"][0]
        antibiotic = d.split("/")[-1].split("___")[0]
        seed = int(sd.split("/")[-1].split("_")[3])

        # lines = open(join(sd, "report.txt"), "r").readlines()
        # lines = lines[np.where([l.startswith("Model (") for l in lines])[0][0]:]
        # model = "".join([l for l in lines if len(l.strip()) > 0])
        
        path = join(destination, species)
        if not exists(path):
            mkdir(path)
        
        path = join(path, antibiotic)
        if not exists(path):
            mkdir(path)

        path = join(path, "repeat_{0:d}".format(seed + 1))
        if not exists(path):
            mkdir(path)

        system("cp /Users/alexandredrouin/Desktop/cart_paper_results/{0!s}/*model*.fasta '{1!s}'".format(sd, path))
        
        detailed_results_path = "../../../../../results/cart_b/{}/{}/repeat_{}/".format(quote(species), quote(antibiotic), seed)

        readme = \
            """
# Model

Species: *{}*

Antibiotic: {}

<img src="./model.png" width=500 height=500 />

For details, please refer to the [results directory]({}).

""".format(species[0].upper() + species[1:].lower(), antibiotic.title(), detailed_results_path)
        open("{0!s}/README.md".format(path), "w").write(readme)
