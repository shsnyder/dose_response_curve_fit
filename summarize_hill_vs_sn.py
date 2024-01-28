import json
import os
import pandas as pd


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


def create_summary(data_folder):
    file_data = []
    for path, directories, files in os.walk(data_folder):
        for file in files:
            name, extension = os.path.splitext(file)
            if extension in [".json"]:
                with open(os.path.join(path, file), "r") as f:
                    fit_data = json.load(f)
            else:
                continue
            median_sn = fit_data["data"]["skewnormal_params"].get("median_sn", None)
            if median_sn is None:
                continue
            file_data.append(
                (
                    name,
                    fit_data["data"]["data_summary"]["sample_id"],
                    fit_data["data"]["data_summary"]["sample_data_id"],
                    fit_data["data"]["data_summary"]["curve_class"],
                    fit_data["data"]["data_summary"]["smiles"],
                    fit_data["data"]["skewnormal_params"].get("median_sn", None),
                    fit_data["data"]["hill_params"]["AC50"],
                    fit_data["data"]["skewnormal_params"]["median_sn"]
                    - fit_data["data"]["hill_params"]["AC50"],
                    100
                    * (
                        fit_data["data"]["skewnormal_params"]["median_sn"]
                        - fit_data["data"]["hill_params"]["AC50"]
                    )
                    / fit_data["data"]["hill_params"]["AC50"],
                )
            )
    return pd.DataFrame(
        file_data,
        columns=(
            "DATA_SET_NAME",
            "SAMPLE_ID",
            "SAMPLE_DATA_ID",
            "CURVE_CLASS",
            "SMILES",
            "AC50_SN",
            "AC50_HILL",
            "AC50_DIFF_SN-H",
            "AC50_PCT_DIFF_SN-H",
        ),
    )


if __name__ == "__main__":
    data_folder = "/Users/shsnyder/Documents/projects/dose_response_data/tox21-ache-p4_curve_fit_lt_3_hill_log"
    summary_output_file = (
        f"{data_folder}/summary_tox21-ache-p4_curve_fit_lt_3_hill_log.csv"
    )
    # output_file = "/Users/shsnyder/Documents/projects/dose_response_data/tox21-ache-p5_curve_fit_all_hill_log/tox21-ache-p5_curve_fit_all_hill_log.csv"

    summary_df = create_summary(data_folder)
    summary_df = summary_df.sort_values(by=["AC50_PCT_DIFF_SN-H"], ascending=True)
    summary_df.to_csv(summary_output_file, sep=",", encoding="utf-8")

    # ############ Update the original data file with skewnormal AC50
    original_data_file_no_path = "tox21-ache-p4.tsv"
    original_data_file = f"{data_folder}/{original_data_file_no_path}"
    updated_data_file = f"{data_folder}/updated_{original_data_file_no_path}.csv"

    original_data_df = read_file(original_data_file)
    filtered_summary = summary_df[["SAMPLE_DATA_ID", "AC50_SN", "AC50_PCT_DIFF_SN-H"]]
    # JSON file has all units in Moles, original datafile in uM
    # For merge with original data convert AC50_SN to uM for consistancy

    filtered_summary["AC50_SN"] = filtered_summary["AC50_SN"].multiply(1e6)

    updated_data_df = original_data_df.merge(filtered_summary, on=["SAMPLE_DATA_ID"])
    updated_data_df.to_csv(updated_data_file, sep=",", encoding="utf-8")
