import numpy as np
import pandas as pd
import os
import json
import sys
from scipy.stats import skewnorm

from dose_response_curve_fit import plot_curves
from dose_response_curve_fit import data_processing
from dose_response_curve_fit import curve_fit


# Getting data
def test_get_dose_response_values():
    dose = [
        -9.008419543,
        -8.658961368,
        -8.309449538,
        -7.959793372,
        -7.610479534,
        -7.260981754,
        -6.91150953,
        -6.562090964,
        -6.212539525,
        -5.862962545,
        -5.513569521,
        -5.1640559,
        -4.814457845,
        -4.465085896,
        -4.115601174,
        -3.985577476,
        -3.872808783,
        -3.783367932,
        -3.709238461,
        -3.645936402,
        -3.590698217,
        -3.541699647,
        -3.497673043,
        -3.457701687,
        -3.421101131,
    ]
    response = [
        0.007647275,
        0.010845427,
        0.015381487,
        0.021817056,
        0.030932996,
        0.043862619,
        0.062188355,
        0.088152324,
        0.124945828,
        0.177046441,
        0.250717137,
        0.354867077,
        0.521887267,
        0.667132923,
        0.77282863,
        0.832953559,
        0.875054172,
        0.911847676,
        0.937811645,
        0.956137381,
        0.969067004,
        0.978182944,
        0.984618513,
        0.989154573,
        0.992352725,
    ]
    return dose, response


def read_file(data_file):
    """_summary_

    Args:
        data_file (_type_): _description_

    Returns:
        _type_: _description_
    """
    extension = os.path.splitext(data_file)[1]
    if extension == ".csv":
        data_df = pd.read_csv(data_file, sep=",", index_col=False)
    elif extension == ".tsv":
        data_df = pd.read_csv(data_file, sep="\t", index_col=False)
    elif extension == ".xlsx":
        data_df = pd.read_excel(data_file)
    else:
        data_df = pd.DataFrame()
    return data_df


def read_tox21_data(data_file, curve_class_limit=None):
    """_summary_

    Args:
        data_file (_type_): _description_

    Returns:
        _type_: _description_
    """
    dose_response = read_file(data_file)

    response_cols = [
        "DATA0",
        "DATA1",
        "DATA2",
        "DATA3",
        "DATA4",
        "DATA5",
        "DATA6",
        "DATA7",
        "DATA8",
        "DATA9",
        "DATA10",
        "DATA11",
        "DATA12",
        "DATA13",
        "DATA14",
    ]
    dose_cols = [
        "CONC0",
        "CONC1",
        "CONC2",
        "CONC3",
        "CONC4",
        "CONC5",
        "CONC6",
        "CONC7",
        "CONC8",
        "CONC9",
        "CONC10",
        "CONC11",
        "CONC12",
        "CONC13",
        "CONC14",
    ]
    other_cols = [
        "SAMPLE_ID",
        "SAMPLE_DATA_ID",
        "AC50",
        "SAMPLE_NAME",
        "SMILES",
        "HILL_COEF",
        "CURVE_CLASS2",
    ]
    dr_filter = response_cols + dose_cols + other_cols
    dr_filtered = dose_response[dr_filter]

    # Making "other_cols" values accessable through a dictionary so that other
    # column names can be used.
    other_col_keys = [
        "id",
        "data_id",
        "ac50",
        "name",
        "smiles",
        "hill_coef",
        "curve_class",
    ]
    other_col_dict = dict(zip(other_col_keys, other_cols))

    if curve_class_limit:
        dr_filtered = dr_filtered[
            np.abs(dr_filtered[other_col_dict["curve_class"]]) < curve_class_limit
        ]

    return dr_filtered, dose_cols, response_cols, other_col_dict


def read_pubchem_data(data_file, structure_file, curve_class_limit=None):
    # Each row contains multiple repeated assys for the same molecular endpoint
    # This needs to be split up  so that each replicant has its own row.
    def split_pubchem_data(dose_response_data):
        data_list = []
        data_columns = dose_response_data.columns

        # Sokit big data set into a shared and replicates data sets
        shared_columns = [item for item in data_columns if len(item.split("-")) == 1]
        shared_df = dose_response_data[shared_columns]

        replicate1_columns = [
            item
            for item in data_columns
            if (len(item.split("-")) == 2) and item.split("-")[1] == "Replicate_1"
        ]
        replicate1_df = dose_response_data[replicate1_columns]
        # Split off the "Replicate" ending
        replicate1_df.rename(
            lambda column: column.split("-")[0], axis="columns", inplace=True
        )

        replicate2_columns = [
            item
            for item in data_columns
            if (len(item.split("-")) == 2) and item.split("-")[1] == "Replicate_2"
        ]
        replicate2_df = dose_response_data[replicate2_columns]
        # Split off the "Replicate" ending
        replicate2_df.rename(
            lambda column: column.split("-")[0], axis="columns", inplace=True
        )

        replicate3_columns = [
            item
            for item in data_columns
            if (len(item.split("-")) == 2) and item.split("-")[1] == "Replicate_3"
        ]
        replicate3_df = dose_response_data[replicate3_columns]
        # Split off the "Replicate" ending
        replicate3_df.rename(
            lambda column: column.split("-")[0], axis="columns", inplace=True
        )

        return shared_df, [replicate1_df, replicate2_df, replicate3_df]

    dose_response_data = read_file(data_file)
    dose_response_structure = read_file(structure_file)
    shared_df, replicate_list = split_pubchem_data(dose_response_data)

    # Find dose levels from replicate file.
    # This assumes all replicate files have the same number and value for dose
    # Column name looks like "'Activity at 0.0000114526 uM'"
    # Assumes uM units
    conc = [
        (float(item[12:-3])) * 10 ** (-6)
        for item in replicate_list[0].columns
        if item[:11] == "Activity at"
    ]

    conc_columns = [
        item for item in replicate_list[0].columns if item[:11] == "Activity at"
    ]

    response_cols = [f"DATA{i}" for i in range(len(conc))]
    dose_cols = [f"CONC{i}" for i in range(len(conc))]

    column_replace_dict = dict(zip(conc_columns, response_cols))
    for i in range(3):
        # This replaces the "Activity" column names with "DATA{i}"
        replicate_list[i].rename(columns=column_replace_dict, inplace=True)
        replicate_list[i]["Replicate_number"] = i + 1

        # Add shared columns from shared_df to each of the replicate dataframes
        replicate_list[i] = pd.concat([shared_df, replicate_list[i]], axis=1)

        # Join the structure dataframe with each replicate to get the Smiles data
        # for each endpoint
        replicate_list[i] = pd.merge(
            replicate_list[i],
            dose_response_structure,
            left_on="PUBCHEM_SID",
            right_on="PUBCHEM_SUBSTANCE_ID",
        )
        pass

    dr_filtered = pd.concat(replicate_list, ignore_index=True, sort=False)

    # Create columns containing the concentrations extracted previously
    for i in range(len(conc)):
        dr_filtered[f"CONC{i}"] = conc[i]

    # Convert LogAC50 to AC50 (in uM)
    dr_filtered["AC50"] = 10 ** (dr_filtered["Fit_LogAC50"]) * 10**6

    other_cols = [
        "PUBCHEM_SID",
        "Replicate_number",
        "AC50",
        "Smiles",
        "Fit_HillSlope",
        "Fit_CurveClass",
    ]

    dr_filter = response_cols + dose_cols + other_cols
    dr_filtered = dr_filtered[dr_filter]

    dr_filtered["PUBCHEM_SID_UNIQUE"] = (
        dr_filtered.PUBCHEM_SID.astype(str)
        + "_"
        + dr_filtered.Replicate_number.astype(str)
    )

    # Making "other_cols" values accessable through a dictionary so that common
    # column names can be used.
    # This is necessary to make the code that accesss the returned dataframe work for
    # two different sources of ata
    generic_other_col_keys = [
        "id",
        "data_id",
        "ac50",
        "name",
        "smiles",
        "hill_coef",
        "curve_class",
    ]
    generic_other_cols = [
        "PUBCHEM_SID",
        "PUBCHEM_SID_UNIQUE",
        "AC50",
        "",
        "Smiles",
        "Fit_HillSlope",
        "Fit_CurveClass",
    ]

    other_col_dict = dict(zip(generic_other_col_keys, generic_other_cols))

    if curve_class_limit:
        dr_filtered = dr_filtered[
            np.abs(dr_filtered[other_col_dict["curve_class"]]) < curve_class_limit
        ]

    pass
    return dr_filtered, dose_cols, response_cols, other_col_dict


if __name__ == "__main__":
    #    Test the curve fit with made up data
    # # Get test data
    # dose, response = test_get_dose_response_values()
    # response = -1*np.array(response)
    # a, loc, scale, ier, mesg = analyze_dose_response(dose, response)

    # # #####################

    data_source = "pubchem"

    if data_source == "tox21":
        #  Tox21 data
        data_file = "/Users/shsnyder/Documents/projects/dose_response_data/tox21-ache-p4_curve_fit_lt_3_hill_log/tox21-ache-p4.tsv"

        dr_filtered, dose_cols, response_cols, other_col_dict = read_tox21_data(
            data_file, curve_class_limit=3
        )
    elif data_source == "pubchem":
        data_file = "/Users/shsnyder/Documents/projects/dose_response_data/HERG_AID_1671200/HERG_AID_1671200_datatable.csv"
        structure_file = "/Users/shsnyder/Documents/projects/dose_response_data/HERG_AID_1671200/PubChem_substanceID_structures_HERG.tsv"
        (
            dr_filtered,
            dose_cols_all,
            response_cols_all,
            other_col_dict,
        ) = read_pubchem_data(data_file, structure_file, curve_class_limit=3)
    else:
        sys.exit(1)

    if data_source == "tox21":
        logdose = True
    elif data_source == "pubchem":
        logdose = True

    for index, row in dr_filtered.iterrows():
        # limit run for test purposes
        # max_index = 100
        # if index > max_index:
        #     break

        # Remove columns where the response is NA
        # FIlter columns
        filtered_pairs = [
            dr_pair
            for dr_pair in zip(dose_cols_all, response_cols_all)
            if not pd.isna(row[dr_pair[1]])
        ]
        # Unzip and restore column definition
        filtered_cols = list(zip(*filtered_pairs))
        dose_cols = list(filtered_cols[0])
        response_cols = list(filtered_cols[1])

        (
            response,
            dose,
            hill_ac50_in_uM,
            sample_name,
            smiles,
            hill_coef,
        ) = data_processing.retrieve_data_metrics(
            row, response_cols, dose_cols, other_col_dict, logdose=logdose
        )

        # Early abort of run if no Hill parameters
        # If data is good it may be possible to do a skewnormal fit
        # but will assume not for the moment
        if not hill_ac50_in_uM or pd.isna(hill_ac50_in_uM):
            continue

        a, loc, scale, ier, mesg = curve_fit.analyze_dose_response(dose, response)

        #  Save metrics in dictionary which will be saved as a JSON file3
        data_stats = {
            "data": {"data_summary": {}, "skewnormal_params": {}, "hill_params": {}}
        }

        data_summary = data_stats["data"]["data_summary"]

        data_summary["sample_id"] = row[other_col_dict["id"]]
        data_summary["sample_data_id"] = row[other_col_dict["data_id"]]
        data_summary["sample_name"] = sample_name
        data_summary["smiles"] = smiles
        data_summary["curve_class"] = row[other_col_dict["curve_class"]]
        data_summary["size_of_dataset"] = len(response)

        hill_params = data_stats["data"]["hill_params"]
        hill_params["hill_coef"] = hill_coef

        # in micromolar
        if hill_ac50_in_uM:
            hill_params["AC50"] = hill_ac50_in_uM * 1e-6
        else:
            hill_params["AC50"] = None

        skewnormal_params = data_stats["data"]["skewnormal_params"]
        skewnormal_params["a"] = a
        skewnormal_params["loc"] = loc
        skewnormal_params["scale"] = scale
        skewnormal_params["ier"] = ier
        skewnormal_params["mesg"] = mesg
        # ################################

        # If fit successful, plot response and fit
        if curve_fit.successful_skewnorm_curvefit(a, loc, scale, ier):
            if hill_ac50_in_uM:
                plot_curves.plot_skewnorm_hill(
                    dose,
                    response,
                    a,
                    loc,
                    scale,
                    hill_params["hill_coef"],
                    hill_params["AC50"],
                    #    plot_file=None)
                    plot_file=f'{data_file}_{row[other_col_dict["data_id"]]}_both.jpg',
                )
            else:
                plot_curves.plot_skewnorm_curve_fit(
                    dose,
                    response,
                    a,
                    loc,
                    scale,
                    plot_file=f'{data_file}_{row[other_col_dict["data_id"]]}_both.jpg',
                )
            mean_sn, var = skewnorm.stats(a, loc, scale, moments="mv")
            median_sn = skewnorm.median(a, loc, scale)
        else:
            mean_sn = var = median_sn = None

        if median_sn:
            skewnormal_params["mean"] = mean_sn.item()
            skewnormal_params["var"] = var.item()
            skewnormal_params["median_sn"] = 10**median_sn

        # Debug
        print(data_stats)
        # End Debug

        with open(
            f'{data_file}_{row[other_col_dict["data_id"]]}_both.json',
            "w",
            encoding="utf-8",
        ) as f:
            json.dump(data_stats, f, ensure_ascii=False, indent=4)
