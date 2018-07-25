"""
Runs kover CART on all the datasets in some directory structure

"""
import json
import numpy as np
import os

from functools import partial
from kover.dataset import KoverDataset
from math import ceil
from multiprocessing import Pool
from progressbar import ProgressBar, Percentage, Bar, Timer, ETA


def find_kover_datasets(dataset_storage_path):
    datasets = []
    for root, directories, filenames in os.walk(dataset_storage_path, followlinks=True):
        for filename in filenames:
            if filename.endswith(".kover"):
                datasets.append((os.path.basename(os.path.normpath(root)), os.path.abspath(os.path.join(root,filename))))
    return datasets


def is_already_computed(split, dataset_uuid, output_path):
    """
    Limitations
    -----------
    Doesn't check that the configuration of the run was the same.

    """

    # Verify the results directory structure
    if not all(os.path.exists(os.path.join(output_path, f)) for f in ["config.json", "model.fasta", "report.txt", "results.json"]):
        return False

    # Check if the results are for the same dataset and the same split
    r = json.load(open(os.path.join(output_path, "results.json")))
    if r["data"]["uuid"] == dataset_uuid and r["data"]["split"] == split:
        return True
    else:
        return False


def get_dataset_uncomputed_splits(dataset_path, output_path, splits=None):
    # Open the datasets
    ds = KoverDataset(dataset_path)

    # Find which splits should be computed
    if splits is None:
        splits = [s.name for s in ds.splits]
    elif isinstance(splits, list):
        splits = splits
    else:
        raise Exception("Splits must be a list or None. If None, will run on all splits.")

    # Check which splits have already been computed
    splits_to_compute = []
    for s in splits:
        if not is_already_computed(s, ds.uuid, output_path):
            splits_to_compute.append(s)

    return splits_to_compute


def generate_command(job_info, criterion, max_depth, min_samples_split, class_importances, hp_choice_method, n_cpu, bound_max_genome_size):
    dataset_path = job_info[0]
    split = job_info[1]
    output_path = job_info[2]

    split_results_path = os.path.join(output_path, split)
    if not is_already_computed(split, KoverDataset(dataset_path).uuid, split_results_path):
        cmd = "kover learn tree --dataset {0!s} " \
             "--split {1!s} " \
             "--criterion {2!s} " \
             "--max-depth {3:d} " \
             "--min-samples-split {4!s} " \
             "--class-importance {5!s} " \
             "--hp-choice {6!s} " \
             "--bound-max-genome-size {7:d} " \
             "--n-cpu {8:d} " \
             "--output-dir {9!s}".format(dataset_path,
                                         split,
                                         criterion,
                                         max_depth,
                                         ' '.join([str(x) for x in min_samples_split]),
                                         ' '.join([str(x) for x in class_importances]),
                                         hp_choice_method,
                                         bound_max_genome_size,
                                         n_cpu,
                                         split_results_path)
        return cmd
    return ""


def generate_dispatchers(commands, n_nodes, n_concurrent_per_node, walltime):
    # Shuffle commands to avoid having all long jobs on a node
    commands = np.array(commands)
    np.random.shuffle(commands)

    # Write a dispatcher for each node
    dispatcher_names = []
    n_jobs_per_node = int(ceil(1.0 * len(commands) / n_nodes))
    for i in xrange(n_nodes):
        node_commands = commands[i * n_jobs_per_node : (i + 1) * n_jobs_per_node]

        # Group commands by n_concurrent_per_node
        grouped_commands = []
        n_jobs_per_serial_block = int(ceil(1.0 * len(node_commands) / n_concurrent_per_node))
        for j in xrange(n_concurrent_per_node):
            grouped_commands.append(node_commands[j * n_jobs_per_serial_block : (j + 1) * n_jobs_per_serial_block])

        dispatcher = \
"""#!/bin/bash
#SBATCH --time={0!s}  # dd-hh:mm
#SBATCH --job-name=amr_cart_bound
#SBATCH --output=node_{1:d}.out
#SBATCH --error=node_{1:d}.err
#SBATCH --account=rrg-marchand

#SBATCH --mem=48G
#SBATCH -c 8

#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

{2!s}

wait

""".format(walltime, i, "\n\n".join("(%s)&" % "; ".join(c) for c in grouped_commands if len(c) > 0))

        f_name = "dispatch_node_%d.msub" % i
        open(f_name, "w").write(dispatcher)
        dispatcher_names.append(f_name)

    return dispatcher_names


if __name__ == "__main__":
    method_name = "cart_bound"
    dataset_storage_path = os.path.abspath(os.path.expandvars("/scratch/adrouin/patric_data/single_species_datasets"))
    results_storage_path = os.path.abspath(os.path.expandvars("/scratch/adrouin/patric_experiments/single_species_datasets/{0!s}".format(method_name)))

    kover_hp_choice_method = "bound"
    kover_criterion = "gini"
    kover_max_depth = 20
    kover_min_samples_split = [2]
    kover_class_importances = [0.25, 0.5, 0.75, 1.0]
    kover_bound_max_genome_size = 10461658
    kover_n_concurrent_hp = 2

    # Cluster parameters
    n_nodes = 50
    n_parallel_per_node = 4
    walltime = "24:00:00"

    # Experiment parameters
    splits_to_compute = ["train_only"] + ["train_0.800_seed_{0:d}_10_folds".format(i) for i in xrange(10)]

    if not os.path.exists(results_storage_path):
        os.mkdir(results_storage_path)

    # Find all datasets in this path
    datasets = find_kover_datasets(dataset_storage_path)

    # Find which computations must be performed on each dataset/split
    print "Finding jobs to compute..."
    jobs = []
    for i, (name, dataset_path) in enumerate(datasets):
        output_path = os.path.join(results_storage_path, name)
        if not os.path.exists(output_path):
            os.mkdir(output_path)
        jobs += [(dataset_path, s, output_path) for s in get_dataset_uncomputed_splits(dataset_path, output_path,
                                                                                       splits_to_compute)]
    commands = np.array([generate_command(job_info,
                                          criterion=kover_criterion,
                                          max_depth=kover_max_depth,
                                          class_importances=kover_class_importances,
                                          min_samples_split=kover_min_samples_split,
                                          hp_choice_method=kover_hp_choice_method,
                                          bound_max_genome_size=kover_bound_max_genome_size,
                                          n_cpu=kover_n_concurrent_hp)
                         for job_info in jobs])
    commands = commands[commands != ""]

    if len(commands) > 0:
        if len(commands) <= n_parallel_per_node * n_nodes:
            n_nodes = int(ceil(1.0 * len(commands) / n_parallel_per_node))
        dispatchers = generate_dispatchers(commands, n_nodes, n_parallel_per_node, walltime)
        launch_command = "; sleep 2; ".join("sbatch %s" % d for d in dispatchers)
        print "Launching (%d jobs)...\n%s" % (len(commands), launch_command)
        os.system(launch_command)
#        os.system("rm %s" % " ".join(d for d in dispatchers))
    else:
        print "All jobs are already computed. Nothing to do."
    print "Done."
