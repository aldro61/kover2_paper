"""
Runs the learner on all the datasets in some directory structure

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
    if not all(os.path.exists(os.path.join(output_path, f)) for f in ["config.json", "results.json"]):
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


def generate_command(job_info, regularizer_values, n_features, n_cpu):
    dataset_path = job_info[0]
    split = job_info[1]
    output_path = job_info[2]

    split_results_path = os.path.join(output_path, split)
    if not is_already_computed(split, KoverDataset(dataset_path).uuid, split_results_path):
        cmd = "python ./learner.py --dataset {0!s} " \
             "--split {1!s} " \
             "--regularizer {2!s} " \
             "--features-k {3:d} " \
             "--ncpu {4:d} " \
             "--output {5!s}".format(dataset_path,
                                         split,
                                         ' '.join([str(x) for x in regularizer_values]),
                                         n_features,
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
#SBATCH --job-name=amr_chi2_l1_logistic
#SBATCH --output=node_{1:d}.out
#SBATCH --error=node_{1:d}.err
#SBATCH --account=rpp-corbeilj

#SBATCH --mem=96G
#SBATCH -c 10

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
    method_name = "chi2_l1_logistic"
    dataset_storage_path = os.path.abspath(os.path.expandvars("/scratch/adrouin/patric_data/single_species_datasets"))
    results_storage_path = os.path.abspath(os.path.expandvars("/scratch/adrouin/patric_experiments/single_species_datasets/{0!s}".format(method_name)))

    regularizer_values = np.logspace(-7, 7, 10)
    n_features = int(1e6)
    n_concurrent_hp = 10

    # Cluster parameters
    n_nodes = 50
    n_parallel_per_node = 1
    walltime = "48:00:00"

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
                                          regularizer_values=regularizer_values,
                                          n_features=n_features,
                                          n_cpu=n_concurrent_hp)
                         for job_info in jobs])
    commands = commands[commands != ""]

    if len(commands) > 0:
        if len(commands) <= n_parallel_per_node * n_nodes:
            n_nodes = int(ceil(1.0 * len(commands) / n_parallel_per_node))
        dispatchers = generate_dispatchers(commands, n_nodes, n_parallel_per_node, walltime)
        launch_command = "; ".join("sbatch %s" % d for d in dispatchers)
        print "Launching (%d jobs)...\n%s" % (len(commands), launch_command)
        os.system(launch_command)
#        os.system("rm %s" % " ".join(d for d in dispatchers))
    else:
        print "All jobs are already computed. Nothing to do."
    print "Done."
