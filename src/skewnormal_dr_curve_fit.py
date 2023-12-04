import pandas as pd
import os
import json
from scipy.stats import skewnorm

from dose_response_curve_fit import plot_curves
from dose_response_curve_fit import data_processing 
from dose_response_curve_fit import curve_fit 

# Getting data
def test_get_dose_response_values():
    dose       = [-9.008419543, -8.658961368, -8.309449538, -7.959793372, -7.610479534,\
                -7.260981754, -6.91150953, -6.562090964, -6.212539525, -5.862962545,\
                -5.513569521, -5.1640559, -4.814457845, -4.465085896, -4.115601174,\
                -3.985577476, -3.872808783, -3.783367932, -3.709238461, -3.645936402,\
                -3.590698217, -3.541699647, -3.497673043, -3.457701687, -3.421101131 ]
    response   = [ 0.007647275, 0.010845427, 0.015381487, 0.021817056, 0.030932996,\
                    0.043862619, 0.062188355, 0.088152324, 0.124945828, 0.177046441, \
                    0.250717137, 0.354867077, 0.521887267, 0.667132923, 0.77282863, \
                    0.832953559, 0.875054172, 0.911847676, 0.937811645, 0.956137381, \
                    0.969067004, 0.978182944, 0.984618513, 0.989154573, 0.992352725 ]
    return dose, response

def read_file(data_file):
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

def read_tox21_data(data_file):
    dose_response = read_file(data_file)

    response_cols = ['DATA0','DATA1', 'DATA2', \
                    'DATA3', 'DATA4', 'DATA5', 'DATA6', \
                    'DATA7', 'DATA8','DATA9', 'DATA10', \
                    'DATA11', 'DATA12', 'DATA13', 'DATA14']
    dose_cols     = ['CONC0', 'CONC1', 'CONC2', 'CONC3',\
                    'CONC4', 'CONC5', 'CONC6', 'CONC7',
                    'CONC8', 'CONC9', 'CONC10', 'CONC11', 
                    'CONC12', 'CONC13', 'CONC14']
    dr_filter     = response_cols + dose_cols + ['SAMPLE_ID','AC50', 'SAMPLE_NAME', 'SMILES', 'HILL_COEF']
    dr_filtered   = dose_response[ dr_filter]  

    return dr_filtered, dose_cols, response_cols  

if __name__ == '__main__':
    #    Test the curve fit with made up data
    # # Get test data 
    # dose, response = test_get_dose_response_values()
    # response = -1*np.array(response)
    # a, loc, scale, ier, mesg = analyze_dose_response(dose, response)

    # # #####################

    #  Tox21 data
    data_file = '/Users/shsnyder/Documents/projects/seizure_project/Tox21_ache_testf/tox21-ache-p5.tsv'

    dr_filtered, dose_cols, response_cols = read_tox21_data(data_file)

    for index, row in dr_filtered.iterrows():

        # limit run for test purposes
        max_index = 100
        if index > max_index: break

        response, dose, hill_ac50_in_uM, sample_name, smiles, hill_coef = data_processing.retrieve_data_metrics(row,
                                                                    response_cols,
                                                                    dose_cols,
                                                                    logdose=True)

        # Early abort of run if no Hill parameters
        # If data is good it may be possible to do a skewnormal fit 
        # but will assume not for the moment
        if not hill_ac50_in_uM or pd.isna(hill_ac50_in_uM):
            continue

        a, loc, scale, ier, mesg = curve_fit.analyze_dose_response(dose, response)

        #  Save metrics in dictionary which will be saved as a JSON file3
        data_stats = {'data':
                            {'data_summary':{},
                            'skewnormal_params':{},
                            'hill_params':{}
                            }
                    }
        
        data_summary = data_stats['data']['data_summary']

        data_summary['sample_id']   = row['SAMPLE_ID']
        data_summary['sample_name'] = sample_name
        data_summary['smiles']      = smiles

        hill_params = data_stats['data']['hill_params']
        hill_params['hill_coef '] = hill_coef

        # in micromolar
        if hill_ac50_in_uM:
            hill_params['AC50']       = hill_ac50_in_uM * 1e-6
        else:
            hill_params['AC50']       = None

        skewnormal_params = data_stats['data']['skewnormal_params']
        skewnormal_params['a']        = a
        skewnormal_params['loc']      = loc
        skewnormal_params['scale']    = scale
        skewnormal_params['ier']      = ier
        skewnormal_params['mesg']     = mesg
        # ################################

        # If fit successful, plot response and fit
        if curve_fit.successful_skewnorm_curvefit(a, loc, scale, ier):
            if hill_ac50_in_uM:
                plot_curves.plot_skewnorm_hill(dose, response,  
                                a, loc, scale,
                                hill_params['hill_coef '], hill_params['AC50'],
                                #    plot_file=None)
                                plot_file=f'{data_file}_{row["SAMPLE_ID"]}_both.jpg' )
            else:
                plot_curves.plot_skewnorm_curve_fit (dose, response,  
                                a, loc, scale,
                                plot_file=f'{data_file}_{row["SAMPLE_ID"]}_both.jpg' )
            mean_sn, var = skewnorm.stats(a, loc, scale, moments='mv')
            median_sn    = skewnorm.median(a, loc, scale)
        else:
            mean_sn = var = median_sn = None

        if median_sn:
            skewnormal_params['mean']        = mean_sn.item()
            skewnormal_params['var']         = var.item()
            skewnormal_params['median_sn']   = 10**median_sn

        # Debug
        print (data_stats)
        # End Debug

        with open(f'{data_file}_{row["SAMPLE_ID"]}_both.json', 'w', encoding='utf-8') as f:
            json.dump(data_stats, f, ensure_ascii=False, indent=4)

