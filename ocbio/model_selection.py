u"""Utilities for Parallel Model Selection with IPython

Author: Olivier Grisel <olivier@ogrisel.com>
Licensed: Simplified BSD

Modified for use in opencast-bio project by
Gavin Gray. <s0805516@sms.ed.ac.uk>
"""
from __future__ import division
from collections import namedtuple
import os

from IPython.parallel import interactive
from IPython.parallel import TaskAborted
from scipy.stats import sem
import numpy as np

from sklearn.utils import check_random_state
from itertools import izip
try:
    # sklearn 0.14+
    from sklearn.grid_search import ParameterGrid
except ImportError:
    # sklearn 0.13
    from sklearn.grid_search import IterGrid as ParameterGrid

from mmap_utils import warm_mmap_on_cv_splits
from mmap_utils import persist_cv_splits


def is_aborted(task):
    return isinstance(getattr(task, u'_exception', None), TaskAborted)


@interactive
def compute_evaluation(model, cv_split_filename, params=None,
    train_fraction=1.0, mmap_mode=u'r', scoring_method='accuracy'):
    u"""Function executed on a worker to evaluate a model on a given CV split"""
    # All module imports should be executed in the worker namespace
    from time import time
    from sklearn.externals import joblib

    X_train, y_train, X_test, y_test = joblib.load(
        cv_split_filename, mmap_mode=mmap_mode)

    # Slice a subset of the training set for plotting learning curves
    n_samples_train = int(train_fraction * X_train.shape[0])
    X_train = X_train[:n_samples_train]
    y_train = y_train[:n_samples_train]

    # Configure the model
    if model is not None:
        model.set_params(**params)

    # Fit model and measure training time
    t0 = time()
    model.fit(X_train, y_train)
    train_time = time() - t0

    # Compute score on training set
    train_score = model.score(X_train, y_train, scoring=scoring_method)

    # Compute score on test set
    test_score = model.score(X_test, y_test, scoring=scoring_method)

    # Wrap evaluation results in a simple tuple datastructure
    return (test_score, train_score, train_time,
            train_fraction, params)


# Named tuple to collect evaluation results
Evaluation = namedtuple(u'Evaluation', (
    u'validation_score',
    u'train_score',
    u'train_time',
    u'train_fraction',
    u'parameters'))


class RandomizedGridSearch(object):
    u""""Async Randomized Parameter search."""

    def __init__(self, load_balanced_view, random_state=0):
        self.task_groups = []
        self.lb_view = load_balanced_view
        self.random_state = random_state
        self._temp_files = []

    def map_tasks(self, f, skip_aborted=True):
        if skip_aborted:
            return [f(task) for task_group in self.task_groups
                            for task in task_group
                            if not is_aborted(task)]
        else:
            return [f(task) for task_group in self.task_groups
                            for task in task_group]

    def abort(self):
        for task_group in self.task_groups:
            for task in task_group:
                if not task.ready() and not is_aborted(task):
                    try:
                        task.abort()
                    except AssertionError:
                        pass
        return self

    def wait(self):
        self.map_tasks(lambda t: t.wait(), skip_aborted=True)
        return self

    def completed(self):
        return sum(self.map_tasks(lambda t: t.ready(), skip_aborted=True))

    def total(self):
        return sum(self.map_tasks(lambda t: 1, skip_aborted=False))

    def progress(self):
        c = self.completed()
        if c == 0:
            return 0.0
        else:
            return float(c) / self.total()

    def reset(self):
        # Abort any other previously scheduled tasks
        self.abort()

        # Schedule a new batch of evaluation tasks
        self.task_groups, self.all_parameters = [], []

        # Collect temporary files:
        for filename in self._temp_files:
            os.unlink(filename)
        del self._temp_files[:]

    def launch_for_splits(self, model, parameter_grid, cv_split_filenames,
        pre_warm=True, collect_files_on_reset=False,scoring_method='accuracy'):
        u"""Launch a Grid Search on precomputed CV splits."""

        # Abort any existing processing and erase previous state
        self.reset()
        self.parameter_grid = parameter_grid

        # Mark the files for garbage collection
        if collect_files_on_reset:
            self._temp_files.extend(cv_split_filenames)

        # Warm the OS disk cache on each host with sequential reads instead
        # of having concurrent evaluation tasks compete for the the same host
        # disk resources later.
        if pre_warm:
            warm_mmap_on_cv_splits(self.lb_view.client, cv_split_filenames)

        # Randomize the grid order
        random_state = check_random_state(self.random_state)
        self.all_parameters = list(ParameterGrid(parameter_grid))
        random_state.shuffle(self.all_parameters)

        for params in self.all_parameters:
            task_group = []

            for cv_split_filename in cv_split_filenames:
                task = self.lb_view.apply(compute_evaluation,
                    model, cv_split_filename, params=params, scoring_method=scoring_method)
                task_group.append(task)

            self.task_groups.append(task_group)

        # Make it possible to chain method calls
        return self

    def launch_for_arrays(self, model, parameter_grid, X, y, n_cv_iter=5, train_size=None,
                          test_size=0.25, pre_warm=True, folder=u".", name=None,
                          random_state=None):
        cv_split_filenames = persist_cv_splits(
            X, y, n_cv_iter=n_cv_iter, train_size=train_size, test_size=test_size,
            name=name, folder=folder, random_state=random_state)
        return self.launch_for_splits(model, parameter_grid,
            cv_split_filenames, pre_warm=pre_warm, collect_files_on_reset=True)

    def find_bests(self, n_top=5):
        u"""Compute the mean score of the completed tasks"""
        mean_scores = []

        for params, task_group in izip(self.all_parameters, self.task_groups):
            evaluations = [Evaluation(*t.get())
                           for t in task_group
                           if t.ready() and not is_aborted(t)]

            if len(evaluations) == 0:
                continue
            val_scores = [e.validation_score for e in evaluations]
            train_scores = [e.train_score for e in evaluations]
            mean_scores.append((np.mean(val_scores), sem(val_scores),
                                np.mean(train_scores), sem(train_scores),
                                params))

        return sorted(mean_scores, reverse=True)[:n_top]

    def report(self, n_top=5):
        bests = self.find_bests(n_top=n_top)
        output = u"Progress: {0:02d}% ({1:03d}/{2:03d})\n".format(
            int(100 * self.progress()), self.completed(), self.total())
        for i, best in enumerate(bests):
            output += (u"\nRank {0}: validation: {1:.5f} (+/-{2:.5f})"
                       u" train: {3:.5f} (+/-{4:.5f}):\n {5}".format(
                       i + 1, *best))
        return output

    def __repr__(self):
        return self.report()

    def boxplot_parameters(self, display_train=False):
        u"""Plot boxplot for each parameters independently"""
        import pylab as pl
        results = [Evaluation(*task.get())
                   for task_group in self.task_groups
                   for task in task_group
                   if task.ready() and not is_aborted(task)]

        n_rows = len(self.parameter_grid)
        pl.figure()
        for i, (param_name, param_values) in enumerate(self.parameter_grid.items()):
            pl.subplot(n_rows, 1, i + 1)
            val_scores_per_value = []
            train_scores_per_value = []
            for param_value in param_values:
                train_scores = [r.train_score for r in results
                                if r.parameters[param_name] == param_value]
                train_scores_per_value.append(train_scores)

                val_scores = [r.validation_score for r in results
                              if r.parameters[param_name] == param_value]
                val_scores_per_value.append(val_scores)

            widths = 0.25
            positions = np.arange(len(param_values)) + 1
            offset = 0
            if display_train:
                offset = 0.175
                pl.boxplot(train_scores_per_value, widths=widths,
                    positions=positions - offset)

            pl.boxplot(val_scores_per_value, widths=widths,
                positions=positions + offset)

            pl.xticks(np.arange(len(param_values)) + 1, param_values)
            pl.xlabel(param_name)
            pl.ylabel(u"Val. Score")

class LearningCurve(RandomizedGridSearch):
    u"""Async Learning Curve processing and plotting."""
    def launch_for_splits(self, model, cv_split_dict, params=None,
        pre_warm=True, collect_files_on_reset=False,scoring_method='accuracy'):
        u"""Launch a Grid Search on precomputed CV splits."""

        # Abort any existing processing and erase previous state
        self.reset()

        for train_size in cv_split_dict.keys():
            for cv_split_filenames in cv_split_dict[train_size]:
                # Mark the files for garbage collection
                if collect_files_on_reset:
                    self._temp_files.extend(cv_split_filenames)

                # Warm the OS disk cache on each host with sequential reads instead
                # of having concurrent evaluation tasks compete for the the same host
                # disk resources later.
                if pre_warm:
                    warm_mmap_on_cv_splits(self.lb_view.client, cv_split_filenames)

        task_group = []
        self.task_dict = {}
        
        for train_size in cv_split_dict.keys():
            task_dict[train_size] = []
            for cv_split_filenames in cv_split_dict[train_size]:
                for cv_split_filename in cv_split_filenames:
                    task = self.lb_view.apply(compute_evaluation,
                        model, cv_split_filename, params=params, scoring_method=scoring_method)
                    task_group.append(task)
                    self.task_dict[train_size].append(task)

                self.task_groups.append(task_group)

        # Make it possible to chain method calls
        return self

    def launch_for_arrays(self, model, X, y, train_sizes, n_cv_iter=5, params=None,
                          test_ratio=0.25, pre_warm=True, folder=u".", name=None,
                          random_state=None):
        cv_split_dict = {}
        for train_size in train_sizes:
                test_size = int(test_ratio*train_size)
                name = "{0}_{1}".format(name,train_size)
                cv_split_dict[train_size] = persist_cv_splits(
                    X, y, n_cv_iter=n_cv_iter, train_size=int(train_size), test_size=test_size,
                    name=name, folder=folder, random_state=random_state)
                
        return self.launch_for_splits(model, cv_split_dict,
             pre_warm=pre_warm, params=params, collect_files_on_reset=True)

    def report(self):
        output = u"Progress: {0:02d}% ({1:03d}/{2:03d})\n".format(
            int(100 * self.progress()), self.completed(), self.total())
        return output

    def plot_curve(self):
        u"""Plot the resulting learning curve."""
        import pylab as pl
        from scipy.stats import sem
        mean_scores = {}
        for train_size in self.task_dict.keys():
            task_group = self.task_dict[train_size]
            evaluations = [Evaluation(*t.get())
                           for t in task_group
                           if t.ready() and not is_aborted(t)]
            if len(evaluations) == 0:
                continue
            train_scores = [e.train_score for e in evaluations]
            test_scores = [e.test_score for e in evaluations]
            mean_scores[train_size] = (np.mean(train_scores), sem(train_scores),
                                    np.mean(test_scores), sem(test_scores))
        # now the mean_scores dictionary contains everything required to build the plot
        trainsizes = sorted(self.task_dict.keys())
        mean_train = [mean_scores[train_size][0] for train_size in train_sizes]
        mean_test = [mean_scores[train_size][2] for train_size in train_sizes]
        train_confidence = [mean_scores[train_size][1]*2 for train_size in train_sizes]
        test_confidence = [mean_scores[train_size][3]*2 for train_size in train_sizes]

        #plot the training scores
        pl.figure()
        pl.fill_between(trainsizes, mean_train - train_confidence, mean_train + train_confidence,
                    color = 'b', alpha = .2)
        pl.plot(trainsizes, mean_train, 'o-k', c='b', label='Train score')

        #plot the test scores
        pl.figure()
        pl.fill_between(trainsizes, mean_test - test_confidence, mean_test + test_confidence,
                    color = 'b', alpha = .2)
        pl.plot(trainsizes, mean_test, 'o-k', c='b', label='Train score')

        #extra annotation
        plt.xlabel('Training set size')
        plt.ylabel('Score')
        plt.xlim(0, max(trainsizes))
        plt.ylim((None, 1.0))  # The best possible score is 1.0
        plt.legend(loc='best')
        plt.title('Main train and test scores +/- 2 standard errors') 

        #should return the data required to recreate the plot
        #so that it can be pickled
        return mean_scores
