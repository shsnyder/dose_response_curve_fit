import json
import os
import pandas as pd


def create_summary(data_folder, output_file):
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
            "AC50_SN",
            "AC50_HILL",
            "AC50_DIFF_SN-H",
            "AC50_PCT_DIFF_SN-H",
        ),
    )


if __name__ == "__main__":
    data_folder = "/Users/shsnyder/Documents/projects/dose_response_data/tox21-ache-p5_curve_fit_all_hill_log"
    output_file = "/Users/shsnyder/Documents/projects/dose_response_data/tox21-ache-p5_curve_fit_all_hill_log/tox21-ache-p5_curve_fit_all_hill_log.csv"

    summary_df = create_summary(data_folder, output_file)
    summary_df = summary_df.sort_values(by=["AC50_PCT_DIFF_SN-H"], ascending=True)
    summary_df.to_csv(output_file, sep=",", encoding="utf-8")
