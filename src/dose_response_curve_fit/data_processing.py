import numpy as np

def retrieve_data_metrics(row, response_cols, dose_cols, logdose=True):

    response = row[response_cols]
    dose = list(row[dose_cols])
    if logdose:
        dose        = list(np.log10(list(dose)))
    median      = row.get('AC50', None)
    sample_name = row['SAMPLE_NAME']
    smiles      = row['SMILES']
    hill_coef   = row.get('HILL_COEF', None)
    return response, dose, median, sample_name, smiles, hill_coef

# Scale response
def scale_response(response):
    start_value = np.mean(response[:3])
    end_value   = np.mean(response[-3:])
    range = end_value - start_value
    scaled_response = (response - start_value)/range
    scaled_response[scaled_response<0] = 0
    return scaled_response, start_value, range

def restore_response(scaled_response, start_value, scale):
    restored_response = scale * scaled_response + start_value
    return restored_response

def hill_eqn(x, hill_coef, ac50):
    return 1/(1 + ((ac50/x)**hill_coef))

def inv_hill_eqn(f, hill_coef, ac50):
    return ac50 * 1/(1/f -1)**(1/hill_coef)