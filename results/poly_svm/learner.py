"""
A polynomial kernel support vector machine for kover datasets (uses the whole feature space)

"""
import argparse
import h5py as h
import json
import numpy as np
import os

from itertools import product
from joblib import delayed, Parallel
from kover.dataset.ds import KoverDataset
from kover.utils import _unpack_binary_bytes_from_ints as unpack
from sklearn.metrics import accuracy_score
from sklearn.svm import SVC


SVC_MAX_ITER = 1e7


def _normalize_gram_matrix(K):
    normX = np.sqrt(K.diagonal())
    return ((K/normX).T/normX).T


def do_cv(regularizer, degree, normalize_kernel, dataset_path, split_name, kernel_matrix):
    # Load the kover dataset
    ds = KoverDataset(dataset_path)
    split = ds.get_split(split_name)

    # Transform the linear kernel matrix into a polynomial kernel matrix
    kernel_matrix = (kernel_matrix.astype(np.float) + 1)**degree

    if normalize_kernel:
        kernel_matrix = _normalize_gram_matrix(kernel_matrix)

    fold_scores = []
    for fold in split.folds:

        # Load data
        fold_train_idx = fold.train_genome_idx[...]; fold_test_idx = fold.test_genome_idx[...]
        y = ds.phenotype.metadata[...]
        K_train = kernel_matrix[fold_train_idx][:, fold_train_idx]; y_train = y[fold_train_idx]
        K_test = kernel_matrix[fold_test_idx][:, fold_train_idx]; y_test = y[fold_test_idx]

        # Evaluate classifier
        model = SVC(C=regularizer, kernel="precomputed", probability=True, max_iter=SVC_MAX_ITER).fit(K_train, y_train)
        fold_scores.append(1.0 - accuracy_score(y_true=y_test, y_pred=model.predict(K_test)))  # risk

    return np.mean(fold_scores)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", type=str, help="The path of the Kover dataset.")
    parser.add_argument("--split", type=str, help="The name of the train/test split to run on.")
    parser.add_argument("--regularizer", type=float, nargs='+', help="The weight of the L2 regularizer.")
    parser.add_argument("--degree", type=float, nargs='+', help="The degree of the polynomial kernel.")
    parser.add_argument("--normalize-kernel", type=int, nargs='+', help="Whether or not to normalize the kernel. (0 = false, 1 = true)")
    parser.add_argument("--ncpu", type=int, help="The number of CPUs to use")
    parser.add_argument("--output", type=str, help="The output directory")
    args = parser.parse_args()

    args.normalize_kernel = [True if x == 1 else False for x in args.normalize_kernel]

    # Load the kernel matrix
    kernel_matrix = h.File(os.path.join(os.path.abspath(os.path.dirname(args.dataset)), "similarity.h5"), "r")["linear_kernel"][...]

    print "Cross-validation"
    print args.normalize_kernel
    hps = list(product(args.regularizer, args.degree, args.normalize_kernel))
    cv_scores = np.array(Parallel(args.ncpu)(delayed(do_cv)(regularizer=reg, degree=deg, normalize_kernel=norm, dataset_path=args.dataset, split_name=args.split, kernel_matrix=kernel_matrix) for reg, deg, norm in hps))
    print "The cv scores are:", cv_scores
    for h, s in zip(hps, cv_scores):
        print h, s

    # Find the best hyperparameters
    best_hps = hps[cv_scores.argmin()]
    best_regularizer = best_hps[0]
    best_degree = best_hps[1]
    best_normalized = best_hps[2]

    print "Training and testing on split with best hyperparameters"
    # Load the data
    ds = KoverDataset(args.dataset)
    split = ds.get_split(args.split)

    # Transform the linear kernel matrix into a polynomial kernel matrix
    kernel_matrix = (kernel_matrix.astype(np.float) + 1)**best_degree

    if best_normalized:
        kernel_matrix = _normalize_gram_matrix(kernel_matrix)

    train_idx = split.train_genome_idx[...]; test_idx = split.test_genome_idx[...]
    y = ds.phenotype.metadata[...]
    K_train = kernel_matrix[train_idx][:, train_idx]; y_train = y[train_idx]
    K_test = kernel_matrix[test_idx][:, train_idx]; y_test = y[test_idx]

    # Create the model
    model = SVC(C=best_hps[0], kernel="precomputed", probability=True, max_iter=SVC_MAX_ITER).fit(K_train, y_train)

    # Evaluate the final model
    from kover.learning.experiments.metrics import _get_binary_metrics
    train_metrics = _get_binary_metrics(predictions=model.predict(K_train), answers=y_train)
    if len(test_idx) > 0:
        test_metrics = _get_binary_metrics(predictions=model.predict(K_test), answers=y_test)
    else:
        test_metrics = {}
    all_predictions_binary = model.predict(kernel_matrix[:, train_idx]).tolist()
    all_predictions_proba = model.predict_proba(kernel_matrix[:, train_idx])[:, 1].tolist()

    # Save results
    results = {"metrics": {"train": train_metrics, "test": test_metrics},
               "data": {"path": args.dataset, "split": args.split, "uuid": ds.uuid},
               "predictions": {"binary": all_predictions_binary, "proba": all_predictions_proba},
               "model": {"n_features": ds.kmer_count},
               "cv": {"score": cv_scores.min(), "hps": {"normalize_kernel": best_normalized,
                                                        "regularizer": best_regularizer,
                                                        "degree": best_degree}}}
    if not os.path.exists(args.output):
        os.mkdir(args.output)
    json.dump(results, open(os.path.join(args.output, "results.json"), "w"))

    # Save configuration
    json.dump(vars(args), open(os.path.join(args.output, "config.json"), "w"))
