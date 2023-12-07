import numpy as np
from scipy.optimize import curve_fit

from dose_response_curve_fit import data_processing


def func_skewnorm_cdf(x, a, loc, scale):
    """_summary_

    Args:
        x (_type_): _description_
        a (_type_): _description_
        loc (_type_): _description_
        scale (_type_): _description_

    Returns:
        _type_: _description_
    """
    from scipy.stats import skewnorm

    return skewnorm.cdf(x, a, loc, scale)


def find_skewnorm_curve_parameters(func, x, y):
    """_summary_

    Args:
        func (_type_): _description_
        x (_type_): _description_
        y (_type_): _description_

    Returns:
        _type_: _description_
    """
    try:
        init_loc = (x[0] + x[-1]) / 2
        popt, pcov, infodict, mesg, ier = curve_fit(
            func,
            x,
            y,
            method="trf",
            full_output=True,
            p0=[0, init_loc, 1],
            # ((lower_bound0, lower_bound1, ..., lower_boundn), (upper_bound0, upper_bound1, ..., upper_boundn))
            bounds=[(-100.0, x[0], -np.inf), (100.0, x[-1], np.inf)],
            # gtol=None,
            xtol=None,
        )
        a, loc, scale = popt
    except Exception as e:
        print(e)
        a = loc = scale = ier = 0
        mesg = "Exceptioo"
    return a, loc, scale, ier, mesg


def analyze_dose_response(dose, response):
    """_summary_

    Args:
        dose (_type_): _description_
        response (_type_): _description_

    Returns:
        _type_: _description_
    """
    trans_response, start, range = data_processing.scale_response(response)

    # Curve fit
    a, loc, scale, ier, mesg = find_skewnorm_curve_parameters(
        func_skewnorm_cdf, dose, trans_response
    )
    print(f"Fit message: {mesg}")

    return a, loc, scale, ier, mesg


def successful_skewnorm_curvefit(a, loc, scale, ier):
    """_summary_

    Args:
        a (_type_): _description_
        loc (_type_): _description_
        scale (_type_): _description_
        ier (_type_): _description_

    Returns:
        _type_: _description_
    """
    successful_curvefit_values = [1, 2, 3, 4]
    reported_curvefilt_failure = ier not in successful_curvefit_values
    if reported_curvefilt_failure:
        return False
    else:
        # Curve fit failed to fit but returned false return code
        # (a, loc, scale) == (0, 0, 0) is for a caught exception
        # (a, loc, scale) == (1, 1, 1) is what the curve_fit function returns when the
        #                              data transition is to drastic to fit
        failed_curve_fit = ((a, loc, scale) == (0, 0, 0)) or (
            (a, loc, scale) == (1, 1, 1)
        )
        return not failed_curve_fit
