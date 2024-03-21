import pandas as pd
import numpy as np

# converts name to GW###### format for visualization purposes
def simplify_name(name):
    if str.startswith(name, "GW"):
        return "GW" + name[2:8]
    return "GW" + name[0:6]

# extracts version number from name
def extract_version(name):
    return name[-1:]

# extracts confidence from catalog.shortName
def extract_group(shortname):
    if "2.1" in shortname:
        return shortname[9:]
    return shortname[7:]

# gets confidence number of group
def extract_confidence(group):
    if group == "confident":
        return 3
    if group == "marginal":
        return 2
    if group == "auxiliary":
        return 1
    return 0

O3_all = pd.read_csv("../Data/GWTC-3.csv")
O3_all["simple_name"] = np.array([simplify_name(name) for name in O3_all["name"]])
O3_all["group"] = np.array([3 for _ in range(O3_all["name"].size)])
O3_all["catalog"] = np.array([extract_group(shortname) for shortname in O3_all["catalog.shortName"]])
O3_all["confidence"] = np.array([extract_confidence(group) for group in O3_all["catalog"]])
O3_all["group_name"] = np.array(["O" + str(group) for group in O3_all["group"]])
O3_all = O3_all[["name", "simple_name", "catalog", "confidence", "version", "group_name", "group", "GPS", "far", "network_matched_filter_snr", "mass_1_source", "mass_2_source", "chirp_mass_source", "redshift", "luminosity_distance", "p_astro", "chi_eff"]]
O3_all = O3_all.rename(columns={"p_astro": "prob", "chi_eff": "effective_spin", "network_matched_filter_snr": "snr", "GPS": "gps", "mass_1_source": "mass1", "mass_2_source": "mass2", "chirp_mass_source": "chirp_mass", "combined_far": "far"})

O2_confident = pd.read_csv("../Data/GWTC-2_1-confident.csv")
O2_confident["simple_name"] = np.array([simplify_name(name) for name in O2_confident["name"]])
O2_confident["group"] = np.array([2 for _ in range(O2_confident["name"].size)])
O2_confident["group_name"] = np.array(["O" + str(group) for group in O2_confident["group"]])
O2_confident["catalog"] = np.array([extract_group(shortname) for shortname in O2_confident["catalog.shortName"]])
O2_confident["confidence"] = np.array([extract_confidence(group) for group in O2_confident["catalog"]])
O2_confident = O2_confident[["name", "simple_name", "catalog", "confidence", "version", "group_name", "group", "GPS", "far", "network_matched_filter_snr", "mass_1_source", "mass_2_source", "chirp_mass_source", "redshift", "luminosity_distance", "p_astro", "chi_eff"]]
O2_confident = O2_confident.rename(columns={"p_astro": "prob", "chi_eff": "effective_spin", "network_matched_filter_snr": "snr", "GPS": "gps", "mass_1_source": "mass1", "mass_2_source": "mass2", "chirp_mass_source": "chirp_mass", "combined_far": "far"})

O2_marginal = pd.read_csv("../Data/GWTC-2_1-marginal.csv")
O2_marginal["simple_name"] = np.array([simplify_name(name) for name in O2_marginal["name"]])
O2_marginal["group"] = np.array([2 for _ in range(O2_marginal["name"].size)])
O2_marginal["group_name"] = np.array(["O" + str(group) for group in O2_marginal["group"]])
O2_marginal["catalog"] = np.array([extract_group(shortname) for shortname in O2_marginal["catalog.shortName"]])
O2_marginal["confidence"] = np.array([extract_confidence(group) for group in O2_marginal["catalog"]])
O2_marginal = O2_marginal[["name", "simple_name", "catalog", "confidence", "version", "group_name", "group", "GPS", "far", "network_matched_filter_snr", "mass_1_source", "mass_2_source", "chirp_mass_source", "redshift", "luminosity_distance", "p_astro", "chi_eff"]]
O2_marginal = O2_marginal.rename(columns={"p_astro": "prob", "chi_eff": "effective_spin", "network_matched_filter_snr": "snr", "GPS": "gps", "mass_1_source": "mass1", "mass_2_source": "mass2", "chirp_mass_source": "chirp_mass", "combined_far": "far"})

O2_auxiliary = pd.read_csv("../Data/GWTC-2_1-auxiliary.csv")
O2_auxiliary["simple_name"] = np.array([simplify_name(name) for name in O2_auxiliary["name"]])
O2_auxiliary["group"] = np.array([2 for _ in range(O2_auxiliary["name"].size)])
O2_auxiliary["group_name"] = np.array(["O" + str(group) for group in O2_auxiliary["group"]])
O2_auxiliary["catalog"] = np.array([extract_group(shortname) for shortname in O2_auxiliary["catalog.shortName"]])
O2_auxiliary["confidence"] = np.array([extract_confidence(group) for group in O2_auxiliary["catalog"]])
O2_auxiliary = O2_auxiliary[["name", "simple_name", "catalog", "confidence", "version", "group_name", "group", "GPS", "far", "network_matched_filter_snr", "mass_1_source", "mass_2_source", "chirp_mass_source", "redshift", "luminosity_distance", "p_astro", "chi_eff"]]
O2_auxiliary = O2_auxiliary.rename(columns={"p_astro": "prob", "chi_eff": "effective_spin", "network_matched_filter_snr": "snr", "GPS": "gps", "mass_1_source": "mass1", "mass_2_source": "mass2", "chirp_mass_source": "chirp_mass", "combined_far": "far"})

O1_all = pd.read_csv("../Data/GWTC-1.csv")
O1_all["simple_name"] = np.array([simplify_name(name) for name in O1_all["name"]])
O1_all["group"] = np.array([1 for _ in range(O1_all["name"].size)])
O1_all["catalog"] = np.array([extract_group(shortname) for shortname in O1_all["catalog.shortName"]])
O1_all["confidence"] = np.array([extract_confidence(group) for group in O1_all["catalog"]])
O1_all["group_name"] = np.array(["O" + str(group) for group in O1_all["group"]])
O1_all = O1_all[["name", "simple_name", "catalog", "confidence", "version", "group_name", "group", "GPS", "far", "network_matched_filter_snr", "mass_1_source", "mass_2_source", "chirp_mass_source", "redshift", "luminosity_distance", "p_astro", "chi_eff"]]
O1_all = O1_all.rename(columns={"p_astro": "prob", "chi_eff": "effective_spin", "network_matched_filter_snr": "snr", "GPS": "gps", "mass_1_source": "mass1", "mass_2_source": "mass2", "chirp_mass_source": "chirp_mass", "combined_far": "far"})

observations = pd.concat([O3_all, O2_confident, O2_marginal, O2_auxiliary, O1_all])

observations["total_mass"] = observations["mass1"] + observations["mass2"]
observations["mass_ratio"] = observations["mass1"] / observations["mass2"]
observations["mass_dos"] = abs(observations["mass2"] - observations["mass1"]) / observations["total_mass"]

observations["is_O1"] = observations["group"] == 1
observations["is_O2"] = observations["group"] == 2
observations["is_O3"] = observations["group"] == 3

confident = observations[observations['confidence'] == 3]

O3 = confident[confident['group'] == 3]
O2 = confident[confident['group'] == 2]
O1 = confident[confident['group'] == 1]

O4_all = pd.read_csv("../Data/real_events_O4_ALL.csv")
O4_all = O4_all[["eventid", "chirp_mass", "combined_far", "mass1", "mass2", "snr", "spin1z", "spin2z", "template_duration", "likelihood"]]
O4_all["group"] = np.array([4 for _ in range(O4_all["eventid"].size)])
O4_all["group_name"] = np.array(["O" + str(group) for group in O4_all["group"]])
O4_all = O4_all.rename(columns={"combined_far": "far"})

O3_all = pd.read_csv("../Data/real_events_O3_ALL.csv")
O3_all = O3_all[["eventid", "chirp_mass", "combined_far", "mass1", "mass2", "snr", "spin1z", "spin2z", "template_duration", "chisq", "likelihood"]]
O3_all["group"] = np.array([3 for _ in range(O3_all["eventid"].size)])
O3_all["group_name"] = np.array(["O" + str(group) for group in O3_all["group"]])
O3_all = O3_all.rename(columns={"combined_far": "far"})

events = pd.concat([O4_all, O3_all])

events["total_mass"] = events["mass1"] + events["mass2"]
events["mass_ratio"] = events["mass1"] / events["mass2"]
events["mass_dos"] = abs(events["mass2"] - events["mass1"]) / events["total_mass"]

events_BBH = events[events["mass2"] > 2.16]
events_BNS = events[events["mass2"] <= 2.16]

O4_events = events[events['group'] == 4]
O4 = O4_events
O3_events = events[events['group'] == 3]

O1_BBH = O1[O1["mass2"] > 2.16]
O2_BBH = O2[O2["mass2"] > 2.16]
O3_BBH = O3[O3["mass2"] > 2.16]
O3_BNS = O3[O3["mass2"] <= 2.16]
O3_events_BBH = O3_events[O3_events["mass2"] > 2.16]
O3_events_BNS = O3_events[O3_events["mass2"] <= 2.16]
O4_BBH = O4_events[O4_events["mass2"] > 2.16]
O4_events_BBH = O4_BBH
O4_BNS = O4_events[O4_events["mass2"] <= 2.16]
O4_events_BNS = O4_BNS

observations["mass_1_source"] = observations["mass1"]
observations["mass_2_source"] = observations["mass2"]
observations["network_matched_filter_snr"] = observations["snr"]
observations["chi_eff"] = observations["effective_spin"]
observations["chirp_mass_source"] = observations["chirp_mass"]
observations["combined_far"] = observations["far"]
observations["p_astro"] = observations["prob"]
observations["network_matched_filter_snr"] = observations["snr"]
observations["M_tot"] = observations["total_mass"]

confident["mass_1_source"] = confident["mass1"]
confident["mass_2_source"] = confident["mass2"]
confident["network_matched_filter_snr"] = confident["snr"]
confident["chi_eff"] = confident["effective_spin"]
confident["chirp_mass_source"] = confident["chirp_mass"]
confident["combined_far"] = confident["far"]
confident["p_astro"] = confident["prob"]
confident["network_matched_filter_snr"] = confident["snr"]
observations["M_tot"] = observations["total_mass"]

O3["mass_1_source"] = O3["mass1"]
O3["mass_2_source"] = O3["mass2"]
O3["network_matched_filter_snr"] = O3["snr"]
O3["chi_eff"] = O3["effective_spin"]
O3["chirp_mass_source"] = O3["chirp_mass"]
O3["combined_far"] = O3["far"]
O3["p_astro"] = O3["prob"]
O3["network_matched_filter_snr"] = O3["snr"]
O3["M_tot"] = O3["total_mass"]
O3_all = O3

O3_BBH["mass_1_source"] = O3_BBH["mass1"]
O3_BBH["mass_2_source"] = O3_BBH["mass2"]
O3_BBH["network_matched_filter_snr"] = O3_BBH["snr"]
O3_BBH["chi_eff"] = O3_BBH["effective_spin"]
O3_BBH["chirp_mass_source"] = O3_BBH["chirp_mass"]
O3_BBH["combined_far"] = O3_BBH["far"]
O3_BBH["p_astro"] = O3_BBH["prob"]
O3_BBH["network_matched_filter_snr"] = O3_BBH["snr"]
O3_BBH["M_tot"] = O3_BBH["total_mass"]
O3_BBH_all = O3_BBH
O3_clean = O3_BBH

O3_BNS["mass_1_source"] = O3_BNS["mass1"]
O3_BNS["mass_2_source"] = O3_BNS["mass2"]
O3_BNS["network_matched_filter_snr"] = O3_BNS["snr"]
O3_BNS["chi_eff"] = O3_BNS["effective_spin"]
O3_BNS["chirp_mass_source"] = O3_BNS["chirp_mass"]
O3_BNS["combined_far"] = O3_BNS["far"]
O3_BNS["p_astro"] = O3_BNS["prob"]
O3_BNS["network_matched_filter_snr"] = O3_BNS["snr"]
O3_BNS["M_tot"] = O3_BNS["total_mass"]
O3_BNS_all = O3_BNS

O2["mass_1_source"] = O2["mass1"]
O2["mass_2_source"] = O2["mass2"]
O2["network_matched_filter_snr"] = O2["snr"]
O2["chi_eff"] = O2["effective_spin"]
O2["chirp_mass_source"] = O2["chirp_mass"]
O2["combined_far"] = O2["far"]
O2["p_astro"] = O2["prob"]
O2["network_matched_filter_snr"] = O2["snr"]
O2["M_tot"] = O2["total_mass"]
O2_all = O2
O2_clean = O2[O2["mass2"] > 2.16]

O1["mass_1_source"] = O1["mass1"]
O1["mass_2_source"] = O1["mass2"]
O1["network_matched_filter_snr"] = O1["snr"]
O1["chi_eff"] = O1["effective_spin"]
O1["chirp_mass_source"] = O1["chirp_mass"]
O1["combined_far"] = O1["far"]
O1["p_astro"] = O1["prob"]
O1["network_matched_filter_snr"] = O1["snr"]
O1["M_tot"] = O1["total_mass"]
O1_all = O1
O1_clean = O1[O1["mass2"] > 2.16]

events["network_matched_filter_snr"] = events["snr"]
events["combined_far"] = events["far"]
events["M_tot"] = events["total_mass"]

events_BBH["network_matched_filter_snr"] = events_BBH["snr"]
events_BBH["combined_far"] = events_BBH["far"]
events_BBH["M_tot"] = events_BBH["total_mass"]

events_BNS["network_matched_filter_snr"] = events_BNS["snr"]
events_BNS["combined_far"] = events_BNS["far"]
events_BNS["M_tot"] = events_BNS["total_mass"]

O3_events["network_matched_filter_snr"] = O3_events["snr"]
O3_events["combined_far"] = O3_events["far"]
O3_events["M_tot"] = O3_events["total_mass"]

O3_events_BBH["network_matched_filter_snr"] = O3_events_BBH["snr"]
O3_events_BBH["combined_far"] = O3_events_BBH["far"]
O3_events_BBH["M_tot"] = O3_events_BBH["total_mass"]

O3_events_BNS["network_matched_filter_snr"] = O3_events_BNS["snr"]
O3_events_BNS["combined_far"] = O3_events_BNS["far"]
O3_events_BNS["M_tot"] = O3_events_BNS["total_mass"]

O4_events["network_matched_filter_snr"] = O4_events["snr"]
O4_events["combined_far"] = O4_events["far"]
O4_events["M_tot"] = O4_events["total_mass"]
O4 = O4_events

O4_events_BBH["network_matched_filter_snr"] = O4_events_BBH["snr"]
O4_events_BBH["combined_far"] = O4_events_BBH["far"]
O4_events_BBH["M_tot"] = O4_events_BBH["total_mass"]
O4_BBH = O4_events_BBH
O4_clean = O4_events_BBH

O4_events_BNS["network_matched_filter_snr"] = O4_events_BNS["snr"]
O4_events_BNS["combined_far"] = O4_events_BNS["far"]
O4_events_BNS["M_tot"] = O4_events_BNS["total_mass"]
O4_BNS = O4_events_BNS

O3_all_predicted = pd.read_csv("../Data Analysis/PredictedData/O3_ALL_predicted.csv")
O3_mock = pd.read_csv("../Data Analysis/PredictedData/O3_mock.csv")