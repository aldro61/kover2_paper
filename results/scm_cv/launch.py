"""
Runs kover SCM (cross-validation) on all the datasets in some directory structure

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


def generate_command(job_info, model_types, p_values, max_rules, max_equiv_rules,
					   hp_choice_method, bound_max_genome_size, random_seed, n_cpu):
	dataset_path = job_info[0]
	split = job_info[1]
	output_path = job_info[2]

	split_results_path = os.path.join(output_path, split)
	if not is_already_computed(split, KoverDataset(dataset_path).uuid, split_results_path):
		cmd = "kover learn scm " \
					  "--dataset %s " \
					  "--split %s " \
					  "--model-type %s " \
					  "--p %s " \
					  "--max-rules %d " \
					  "--max-equiv-rules %d " \
					  "--hp-choice %s " \
					  "--bound-max-genome-size %d " \
					  "--random-seed %d " \
					  "--n-cpu %d " \
					  "--output-dir %s" % (dataset_path,
										   split,
										   ' '.join(model_types),
										   ' '.join([str(p) for p in p_values]),
										   max_rules,
										   max_equiv_rules,
										   hp_choice_method,
										   bound_max_genome_size,
										   random_seed,
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
#SBATCH --job-name=amr_scm_cv
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
	dataset_storage_path = os.path.abspath(os.path.expandvars("/scratch/adrouin/patric_data/single_species_datasets"))
	results_storage_path = os.path.abspath(os.path.expandvars("/scratch/adrouin/patric_experiments/single_species_datasets/scm_cv"))

	kover_model_types = ["conjunction", "disjunction"]
	kover_p_values = [1.0, 1.778, 3.162, 5.623, 10.0]
	kover_max_rules = 20
	kover_max_equiv_rules = 100000
	kover_hp_choice_method = "cv"
	kover_bound_max_genome_size = 10461658
	kover_random_seed = 42
	kover_n_concurrent_hp = 2

	# Colosse parameters
	n_nodes = 55
	n_parallel_per_node = 4
	walltime = "24:00:00"

	# Experiment parameters
	splits_to_compute = ["train_only", "train_0.800_seed_0_10_folds", "train_0.800_seed_1_10_folds", "train_0.800_seed_2_10_folds", "train_0.800_seed_3_10_folds", "train_0.800_seed_4_10_folds", "train_0.800_seed_5_10_folds", "train_0.800_seed_6_10_folds", "train_0.800_seed_7_10_folds", "train_0.800_seed_8_10_folds", "train_0.800_seed_9_10_folds"]

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

	print "Generating dispatchers..."
	commands = np.array([generate_command(job_info, model_types=kover_model_types,
													p_values=kover_p_values,
													max_rules=kover_max_rules,
													max_equiv_rules=kover_max_equiv_rules,
													hp_choice_method=kover_hp_choice_method,
													bound_max_genome_size=kover_bound_max_genome_size,
													random_seed=kover_random_seed,
													n_cpu=kover_n_concurrent_hp) for job_info in jobs])
	commands = commands[commands != ""]

	if len(commands) > 0:
		if len(commands) <= n_parallel_per_node * n_nodes:
			n_nodes = int(ceil(1.0 * len(commands) / n_parallel_per_node))
		dispatchers = generate_dispatchers(commands, n_nodes, n_parallel_per_node, walltime)
		launch_command = "; ".join("sbatch %s" % d for d in dispatchers)
		print "Launching (%d jobs)...\n%s" % (len(commands), launch_command)
		os.system(launch_command)
		#os.system("rm %s" % " ".join(d for d in dispatchers))
	else:
		print "All jobs are already computed. Nothing to do."
	print "Done."
