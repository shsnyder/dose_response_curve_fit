import numpy as np


def retrieve_data_metrics(row, response_cols, dose_cols, logdose=True):
    """_summary_

    Args:
        row (_type_): _description_
        response_cols (_type_): _description_
        dose_cols (_type_): _description_
        logdose (bool, optional): _description_. Defaults to True.

    Returns:
        _type_: _description_
    """
    response = row[response_cols]
    dose = list(row[dose_cols])
    if logdose:
        dose = list(np.log10(list(dose)))
    median = row.get("AC50", None)
    sample_name = row["SAMPLE_NAME"]
    smiles = row["SMILES"]
    hill_coef = row.get("HILL_COEF", None)
    return response, dose, median, sample_name, smiles, hill_coef


# Scale response
def scale_response(response):
    """_summary_

    Args:
        response (_type_): _description_

    Returns:
        _type_: _description_
    """
    start_value = np.mean(response[:3])
    end_value = np.mean(response[-3:])
    range = end_value - start_value
    scaled_response = (response - start_value) / range
    scaled_response[scaled_response < 0] = 0
    return scaled_response, start_value, range


def restore_response(scaled_response, start_value, scale):
    """_summary_

    Args:
        scaled_response (_type_): _description_
        start_value (_type_): _description_
        scale (_type_): _description_

    Returns:
        _type_: _description_
    """
    restored_response = scale * scaled_response + start_value
    return restored_response


def hill_eqn(x, hill_coef, ac50):
    """_summary_

    Args:
        x (_type_): _description_
        hill_coef (_type_): _description_
        ac50 (_type_): _description_

    Returns:
        _type_: _description_
    """
    return 1 / (1 + ((ac50 / x) ** hill_coef))


def inv_hill_eqn(f, hill_coef, ac50):
    """_summary_

    Args:
        f (_type_): _description_
        hill_coef (_type_): _description_
        ac50 (_type_): _description_

    Returns:
        _type_: _description_
    """
    return ac50 * 1 / (1 / f - 1) ** (1 / hill_coef)
