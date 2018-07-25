"""
A L1-regularized logistic regression with chi-square feature selection for kover datasets

"""
import argparse
import gc
import h5py as h
import json
import numpy as np
import os

from joblib import delayed, Parallel
from kover.dataset.ds import KoverDataset
from kover.learning.experiments.metrics import _get_binary_metrics
from kover.utils import _unpack_binary_bytes_from_ints as unpack
from sklearn.metrics import accuracy_score
from sklearn.svm import LinearSVC
from sklearn.linear_model import LogisticRegression

def do_cv(regularizer, dataset_path, split_name):
    # Load the kover dataset
    ds = KoverDataset(dataset_path)
    split = ds.get_split(split_name)

    fold_scores = []
    for fold in split.folds:
        # Select features
        print("... loading HDF5 data")
        p_values = h.File(os.path.join(os.path.abspath(os.path.dirname(args.dataset)),
                                       "chi2_feature_pvalues.h5"),
                          "r")[args.split]["{0!s}_pvalues".format(fold.name)][...]
        print("... extracting top k feature indices")
        feature_idx = p_values.argsort()[: args.features_k]
        feature_idx = np.sort(feature_idx).tolist()  # Requirement of h5py, idx must be a sorted list

        # Load data
        print("... loading the data")
        fold_train_idx = fold.train_genome_idx[...]; fold_test_idx = fold.test_genome_idx[...]
        X = np.vstack((unpack(ds.kmer_matrix[i][feature_idx].reshape([1, -1])) for i in range(ds.kmer_matrix.shape[0])))
        y = ds.phenotype.metadata[...]
        X_train = X[fold_train_idx]; y_train = y[fold_train_idx]
        X_test = X[fold_test_idx]; y_test = y[fold_test_idx]

        # Evaluate classifier
        print("... fitting model")
        model = LogisticRegression(C=regularizer, penalty='l1', dual=False, class_weight='balanced', random_state=42).fit(X_train, y_train)
        fold_scores.append(1.0 - accuracy_score(y_true=y_test, y_pred=model.predict(X_test)))  # risk
        del model; gc.collect()

    return np.mean(fold_scores)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", type=str, help="The path of the Kover dataset.", required=True)
    parser.add_argument("--split", type=str, help="The name of the train/test split to run on.", required=True)
    parser.add_argument("--regularizer", type=float, nargs='+', help="The weight of the regularizer.", required=True)
    parser.add_argument("--features-k", type=int, help="The number of features to consider (k smallest p-values)", required=True)
    parser.add_argument("--ncpu", type=int, help="The number of CPUs to use", required=True)
    parser.add_argument("--output", type=str, help="The output directory", required=True)
    args = parser.parse_args()

    try:
        h.File(os.path.join(os.path.abspath(os.path.dirname(args.dataset)),
                                       "chi2_feature_pvalues.h5")).close()
    except:
        print "Error: the p-value file for dataset", args.dataset, "does not exist. Bye."
        exit()

    print("Cross-validation")
    regularizer_values = np.array(args.regularizer)
    cv_scores = np.array(Parallel(args.ncpu)(delayed(do_cv)(regularizer=reg, dataset_path=args.dataset, split_name=args.split) for reg in regularizer_values))

    # Print the results
    for hps, score in zip(regularizer_values, cv_scores):
        print(hps, score)
    print()

    print("Training and testing on split with best hyperparameters")
    # Load the data
    ds = KoverDataset(args.dataset)
    split = ds.get_split(args.split)

    # Select features
    p_values = h.File(os.path.join(os.path.abspath(os.path.dirname(args.dataset)),
                                   "chi2_feature_pvalues.h5"),
                      "r")[args.split]["main_pvalues"][...]
    feature_idx = p_values.argsort()[: args.features_k]
    feature_idx = np.sort(feature_idx).tolist()  # Requirement of h5py, idx must be a sorted list

    train_idx = split.train_genome_idx[...]; test_idx = split.test_genome_idx[...]
    X = np.vstack((unpack(ds.kmer_matrix[i][feature_idx].reshape([1, -1])) for i in range(ds.kmer_matrix.shape[0])))
    y = ds.phenotype.metadata[...]
    X_train = X[train_idx]; y_train = y[train_idx]
    X_test = X[test_idx]; y_test = y[test_idx]

    # Find the best hyperparameters and instanciate a model
    best_regularizer = regularizer_values[cv_scores.argmin()]
    model = LogisticRegression(C=best_regularizer, penalty='l1', dual=False, class_weight='balanced', random_state=42).fit(X_train, y_train)

    # Evaluate the final model
    from kover.learning.experiments.metrics import _get_binary_metrics
    train_metrics = _get_binary_metrics(predictions=model.predict(X_train), answers=y_train)
    if len(test_idx) > 0:
        test_metrics = _get_binary_metrics(predictions=model.predict(X_test), answers=y_test)
    else:
        test_metrics = {}
    all_predictions_binary = model.predict(X).tolist()
    all_predictions_proba = model.predict_proba(X)[:, 1].tolist()

    # Save results
    results = {"metrics": {"train": train_metrics, "test": test_metrics},
               "data": {"path": args.dataset, "split": args.split, "uuid": ds.uuid},
               "predictions": {"binary": all_predictions_binary, "proba": all_predictions_proba},
               "model": {"n_features": (~np.isclose(model.coef_, 0.)).sum()}}
    if not os.path.exists(args.output):
        os.mkdir(args.output)
    json.dump(results, open(os.path.join(args.output, "results.json"), "w"))

    # Save configuration
    json.dump(vars(args), open(os.path.join(args.output, "config.json"), "w"))
