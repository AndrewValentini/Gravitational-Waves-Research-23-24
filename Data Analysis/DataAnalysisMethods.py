# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#  This files contains a library of useful data anaylsis and  #
#  visualization functions. To use any of these functions     #
#  simply call:                                               #
#                                                             #
#  import DataAnalysisMethods as dam                          #
#                                                             #
#  Warning: most methods have not been tested yet.            #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import numpy as np
import pandas as pd
import math


# Calculates the number of bins needed for the given dataframe using Sturge's Rule
def sturges_rule(df : pd.DataFrame) -> int:
    return math.ceil(np.log2(len(df)) + 1)


# Seperates a dataset into neutron star events and black hole events by mass
def seperate_events(df : pd.DataFrame, mass_1_column_name="mass_1_source", mass_2_column_name="mass_2_source", neuron_star_mass_limit=2.16):
    df_neuron_stars = pd.DataFrame(columns=df.columns)
    df_black_holes = pd.DataFrame(columns=df.columns)

    for _, row in df.iterrows():
        if row[mass_1_column_name] < neuron_star_mass_limit or row[mass_2_column_name] < neuron_star_mass_limit:
            df_neuron_stars = pd.concat([df_neuron_stars, pd.DataFrame([row])], ignore_index=True)
        else:
            df_black_holes = pd.concat([df_black_holes, pd.DataFrame([row])], ignore_index=True)
    
    return df_neuron_stars, df_black_holes


# Seperates dataset into an array of bins. If number of bins is not specified Sturge's Rule is used.
# Returns bins along with the respective ranges of those bins
def bin(df : pd.DataFrame, bin_feature_name : str, n_bins : int = None):
    max_feature = max(df[bin_feature_name])
    min_feature = min(df[bin_feature_name])

    total_bins = n_bins
    if n_bins is None:
        total_bins = sturges_rule(df)

    bin_length = (max_feature - min_feature) / total_bins

    bins = []
    bin_ranges = []

    for bin_index in range(total_bins):
        low = min_feature + bin_length * bin_index
        bin = df[df[bin_feature_name] >= low]

        high = low + bin_length
        bin = bin[bin[bin_feature_name] <= high]
        
        if len(bin[bin_feature_name]) != 0:
            bins.append(bin)
            bin_ranges.append((low, high))

    return np.array(bins), np.array(bin_ranges)


# Seperates dataset into an array of folds. Returns folds along with the respective ranges of those folds
def fold(df : pd.DataFrame, fold_feature_name : str, n_fold : int = 5, ascending : bool = False):
    df_ordered = df.sort_values(by=[fold_feature_name], ascending=ascending)

    folds = []
    fold_ranges = []

    for fold_index in range(n_fold):
        start_index = math.floor(fold_index * len(df_ordered[fold_feature_name]) / n_fold)
        end_index = math.floor((fold_index + 1) * len(df_ordered[fold_feature_name]) / n_fold - 1)
        if fold_index == n_fold - 1:
            end_index = len(df_ordered[fold_feature_name]) - 1

        folds.append(df_ordered[:,start_index:end_index])
        fold_ranges.append((min(folds[fold_feature_name, fold_index]), max(folds[fold_feature_name, fold_index])))

    return np.array(folds), np.array(fold_ranges)


_magic_colors = ["indianred", "peru", "orange", "gold", "yellowgreen", "lime", "darkcyan", "lightseagreen", "dodgerblue", "slateblue", "darkviolet", "purple", "magenta", "palevioletred", "pink"]

# Get some evenly spaced number of colors from _magic_colors
def colors(n_colors : int):
    assert len(_magic_colors) >= n_colors

    exclude = len(_magic_colors) - n_colors

    colors = _magic_colors.copy()
    for i in range(exclude):
        index = math.ceil(i / exclude)
        while not _magic_colors[index] in colors:
            index += 1
            if index == len(_magic_colors):
                index = 0
        colors.remove(_magic_colors[index])

    return np.array(colors, dtype=str)


# Fines the slope of the line of best fit with a given y-intercept. Returns that slope.
def fit_slope(x_feature, y_feature : str, y_intercept : int = 0) -> float:
    y_reduced = y_feature - y_intercept
    
    x = x_feature[:,np.newaxis]
    a, _, _, _ = np.linalg.lstsq(x, y_reduced, rcond=None)
    
    return a[0]


# Calculates the expected y-intercept of a set of bins using the average of the set of
# y-intercepts found by np.polyfit. Sets with only one element are ignored.
def y_intercept(x_feature_name : str, y_feature_name : str, *sets : pd.DataFrame, method : str = "mean") -> float:
    assert method is "mean" or method is "median"

    intercepts = list()

    for set in sets:
        if len(set[y_feature_name]) > 1:
            b = 0
            _, b = np.polyfit(set[x_feature_name], set[y_feature_name], 1)
            intercepts.append(b)
    
    if method is "median":
        return np.median(intercepts)
    return np.mean(intercepts)


# Calculates the expected y-intercept using bins, where each bin is a subset of the
# dataset. Acts simular to median_intercept except that the function itself creates
# the bins of the dataframe. Needs to be given the feature to bin by. If n_bins is None
# Sturge's Rule is used to find the number of bins.
def bin_intercept(df : pd.DataFrame, bin_feature_name : str, x_feature_name : str, y_feature_name : str, n_bins : int = None, method : str = "mean") -> float:
    assert method is "mean" or method is "median"

    intercepts = list()
    bins, _ = bin(df, bin_feature_name, n_bins=n_bins)

    for _bin in bins:
        if len(_bin[y_feature_name]) > 1:
            b = 0
            _, b = np.polyfit(_bin[x_feature_name], _bin[y_feature_name], 1)
            intercepts.append(b)
    
    if method is "median":
        return np.median(intercepts)
    return np.mean(intercepts)


# Calculates the expected y-intercept using fold, where each fold is an equal subset
# of the dataset. Needs to be given the feature to fold by.
def fold_intercept(df : pd.DataFrame, fold_feature_name : str, x_feature_name : str, y_feature_name : str, n_folds : int = 5, ascending : bool = False, method : str = "mean") -> float:
    assert method is "mean" or method is "median"

    intercepts = list()
    folds, _ = fold(df, fold_feature_name, n_folds=n_folds, ascending=ascending)

    for _fold in folds:
        b = 0
        _, b = np.polyfit(_fold[x_feature_name], _fold[y_feature_name], 1)
        intercepts.append(b)
    
    if method is "median":
        return np.median(intercepts)
    return np.mean(intercepts)


# Fits a set of subsets with individual linear lines with a constant y-intercept.
def magic_fit(x_feature_name : str, y_feature_name : str, *sets : pd.DataFrame, method : str = "mean"):
    slopes = np.zeros(len(sets), dtype=float)

    _y_intercept = y_intercept(x_feature_name, y_feature_name, sets, method=method)

    for index in range(len(sets)):
        slopes[index] = fit_slope(sets[index], x_feature_name, y_feature_name, y_intercept=_y_intercept)

    return slopes, _y_intercept