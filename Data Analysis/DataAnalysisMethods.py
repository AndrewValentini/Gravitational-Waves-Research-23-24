# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#  This files contains a library of useful data anaylsis and  #
#  visualization functions. To use any of these functions     #
#  simply call:                                               #
#                                                             #
#  import DataAnalysisMethods as dam                          #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import numpy as np
import pandas as pd
import matplotlib as plt
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
def bin(df : pd.DataFrame, bin_feature_name : str, bins : int = None):
    max_feature = max(df[bin_feature_name])
    min_feature = min(df[bin_feature_name])

    total_bins = bins
    if bins is None:
        total_bins = sturges_rule(df)

    bin_length = (max_feature - min_feature) / total_bins

    bins = []
    bin_ranges = []

    for i in range(total_bins):
        low = min_feature + bin_length * i
        bin = df[df[bin_feature_name] >= low]

        high = low + bin_length
        bin = bin[bin[bin_feature_name] <= high]
        
        if len(bin[bin_feature_name]) != 0:
            bins.append(bin)
            bin_ranges.append((low, high))

    return np.array(bins), np.array(bin_ranges)


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
