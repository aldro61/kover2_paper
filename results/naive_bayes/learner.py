"""
A naive Bayes learner for kover datasets

"""
import argparse
import gc
import h5py as h
import json
import numpy as np
import os

from functools import partial
from joblib import delayed, Parallel
from kover.dataset.ds import KoverDataset
from kover.learning.experiments.metrics import _get_binary_metrics
from kover.utils import _unpack_binary_bytes_from_ints as unpack
from sklearn.metrics import accuracy_score


def tiebreaker(model_type, feature_idx, thresholds, kind, rule_risks):
    # 0 = greater, 1 = less_equal
    print thresholds, kind
    return 0
    tie_rule_risks = rule_risks[best_utility_idx]
    if model_type == "conjunction":
        result = best_utility_idx[np.isclose(tie_rule_risks, tie_rule_risks.min())]
    else:
        # Use max instead of min, since in the disjunction case the risks = 1.0 - conjunction risks (inverted ys)
        result = best_utility_idx[np.isclose(tie_rule_risks, tie_rule_risks.max())]
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", type=str, help="The path of the Kover dataset.", required=True)
    parser.add_argument("--split", type=str, help="The name of the train/test split to run on.", required=True)
    parser.add_argument("--output", type=str, help="The output directory", required=True)
    args = parser.parse_args()

    print "Training and testing on split with best hyperparameters"
    # Load the data
    ds = KoverDataset(args.dataset)
    split = ds.get_split(args.split)

    train_idx = split.train_genome_idx[...]; test_idx = split.test_genome_idx[...]
    # X = np.vstack((unpack(ds.kmer_matrix[i][feature_idx].reshape([1, -1])) for i in range(ds.kmer_matrix.shape[0])))
    y = ds.phenotype.metadata[...]

    y_train = y[train_idx]
    y_test = y[test_idx]

    # Class priors
    print "Computing class priors"
    log_pos_prior = np.log(1.0 * (y_train == 1).sum() / len(y_train))
    log_neg_prior = np.log(1.0 * (y_train == 0).sum() / len(y_train))

    from kover.learning.common.rules import KmerRuleClassifications
    kmer_classifications = KmerRuleClassifications(dataset=ds.kmer_matrix, n_rows=ds.genome_count)
    train_pos_idx = train_idx[y[train_idx] == 1]
    train_neg_idx = train_idx[y[train_idx] == 0]

    print "Learning positive class probabilities"
    # XXX: Adding constants is like adding two examples, one that has the k-mer and one that doesn't
    p_pos = 1. * (kmer_classifications.sum_rows(train_pos_idx)[: ds.kmer_count].astype(np.uint) + 1) / (len(train_pos_idx) + 2)
    p_pos = np.vstack((1. - p_pos, p_pos))
    log_p_pos = np.log(p_pos)
    print "Learning negative class probabilities"
    p_neg = 1. * (kmer_classifications.sum_rows(train_neg_idx)[: ds.kmer_count].astype(np.uint) + 1) / (len(train_neg_idx) + 2)
    p_neg = np.vstack((1. - p_neg, p_neg))
    log_p_neg = np.log(p_neg)

    if np.any(np.isclose(p_pos, 0.)):
        print "Zero value encountered in positive probas for dataset", args.dataset
        raise ValueError()
    if np.any(np.isclose(p_neg, 0.)):
        print "Zero value encountered in negative probas for dataset", args.dataset
        raise ValueError()

    print "Training completed."

    # Predict on all examples
    example_predictions = None
    for i in xrange(ds.kmer_matrix.shape[0]):

        ex_block_pos_likelihood = np.zeros(64 if ds.kmer_matrix.dtype == np.uint64 else 32)
        ex_block_neg_likelihood = np.zeros(64 if ds.kmer_matrix.dtype == np.uint64 else 32)

        block_size = 1000000
        for block_start_idx in xrange(0, ds.kmer_count, block_size):
            data = unpack(ds.kmer_matrix[i][block_start_idx : block_start_idx + block_size].reshape([1, -1]))

            block_log_p_pos = log_p_pos[:, block_start_idx : block_start_idx + block_size]
            block_log_p_neg = log_p_neg[:, block_start_idx : block_start_idx + block_size]

            tmp_pos = np.sum(block_log_p_pos[data, range(0, block_log_p_pos.shape[1])], axis=1)
            tmp_neg = np.sum(block_log_p_neg[data, range(0, block_log_p_pos.shape[1])], axis=1)

            ex_block_pos_likelihood += tmp_pos
            ex_block_neg_likelihood += tmp_neg

            print "Genomes:", i + 1, "/", ds.kmer_matrix.shape[0], "   Kmer:", block_start_idx + 1, "/", ds.kmer_count

        pos_score = log_pos_prior + ex_block_pos_likelihood
        neg_score = log_neg_prior + ex_block_neg_likelihood
        binary_predictions = (pos_score > neg_score).astype(np.uint8)
        scores = np.vstack((neg_score, pos_score)).T

        if example_predictions is None:
            example_predictions = binary_predictions
            example_class_scores = scores
        else:
            example_predictions = np.hstack((example_predictions, binary_predictions))
            example_class_scores = np.hstack((example_class_scores, scores))

    example_predictions = example_predictions[: ds.genome_count]
    example_class_scores = example_class_scores[: ds.genome_count]

    # Evaluate the final model
    from kover.learning.experiments.metrics import _get_binary_metrics
    train_metrics = _get_binary_metrics(predictions=example_predictions[train_idx], answers=y_train)
    if len(test_idx) > 0:
        test_metrics = _get_binary_metrics(predictions=example_predictions[test_idx], answers=y_test)
    else:
        test_metrics = {}
    all_predictions_binary = example_predictions.tolist()
    all_predictions_proba = (example_class_scores[:, 1] / (example_class_scores[:, 0] + example_class_scores[:, 1])).tolist()  # this is a sketchy way to get a proba

    # Save results
    results = {"metrics": {"train": train_metrics, "test": test_metrics},
               "data": {"path": args.dataset, "split": args.split, "uuid": ds.uuid},
               "predictions": {"binary": all_predictions_binary, "proba": all_predictions_proba},
               "model": {"n_features": ds.kmer_count}}
    if not os.path.exists(args.output):
        os.mkdir(args.output)
    json.dump(results, open(os.path.join(args.output, "results.json"), "w"))

    # Save configuration
    json.dump(vars(args), open(os.path.join(args.output, "config.json"), "w"))
