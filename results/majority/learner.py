"""
A linear kernel support vector machine for kover datasets (uses the whole feature space)

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


def do_cv(regularizer, normalize_kernel, dataset_path, split_name, kernel_matrix):
    # Load the kover dataset
    ds = KoverDataset(dataset_path)
    split = ds.get_split(split_name)

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


class DummyMajority(object):
    def fit(self, x, y):
        self.pred_ = int(np.mean(y) > 0.5)
        return self

    def predict(self, X):
        return np.ones(X.shape[0], dtype=np.uint8) * self.pred_


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", type=str, help="The path of the Kover dataset.")
    parser.add_argument("--split", type=str, help="The name of the train/test split to run on.")
    parser.add_argument("--output", type=str, help="The output directory")
    args = parser.parse_args()

    # Load the data
    ds = KoverDataset(args.dataset)
    split = ds.get_split(args.split)

    train_idx = split.train_genome_idx[...]; test_idx = split.test_genome_idx[...]
    X = np.zeros(shape=(ds.genome_count, 1))
    y = ds.phenotype.metadata[...]
    X_train = X[train_idx]; y_train = y[train_idx]
    X_test = X[test_idx]; y_test = y[test_idx]

    # Create the model
    model = DummyMajority().fit(X_train, y_train)

    # Evaluate the final model
    from kover.learning.experiments.metrics import _get_binary_metrics
    train_metrics = _get_binary_metrics(predictions=model.predict(X_train), answers=y_train)
    if len(test_idx) > 0:
        test_metrics = _get_binary_metrics(predictions=model.predict(X_test), answers=y_test)
    else:
        test_metrics = {}
    all_predictions = model.predict(X).tolist()

    # Save results
    results = {"metrics": {"train": train_metrics, "test": test_metrics},
               "data": {"path": args.dataset, "split": args.split, "uuid": ds.uuid},
               "predictions": all_predictions,
               "model": {"n_features": 0, "predicted_value": model.pred_}}
    if not os.path.exists(args.output):
        os.mkdir(args.output)
    json.dump(results, open(os.path.join(args.output, "results.json"), "w"))

    # Save configuration
    json.dump(vars(args), open(os.path.join(args.output, "config.json"), "w"))
