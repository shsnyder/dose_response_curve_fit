import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from scipy.stats import skewnorm

from dose_response_curve_fit import data_processing


# Plotting functions
def plot_skewnorm(x, y, a, loc, scale, plot_file=None):
    """_summary_

    Args:
        x (_type_): _description_
        y (_type_): _description_
        a (_type_): _description_
        loc (_type_): _description_
        scale (_type_): _description_
        plot_file (_type_, optional): _description_. Defaults to None.
    """
    # Plot original data

    x1 = np.linspace(
        skewnorm.ppf(0.00001, a, loc, scale), skewnorm.ppf(0.999, a, loc, scale), 100
    )
    sns.lineplot(
        x="dose_pred",
        y="response_pred",
        # data=pd.DataFrame({'dose_pred':x1,
        data=pd.DataFrame(
            {"dose_pred": x1, "response_pred": skewnorm.cdf(x1, a, loc, scale)}
        ),
        color="blue",
        lw=5,
        alpha=0.6,
        label="Skewnormal CDF fit",
    )

    if plot_file:
        plt.savefig(plot_file)
    else:
        plt.show()
    plt.clf()


def plot_skewnorm_curve_fit(x, y, a, loc, scale, plot_file=None):
    """_summary_

    Args:
        x (_type_): _description_
        y (_type_): _description_
        a (_type_): _description_
        loc (_type_): _description_
        scale (_type_): _description_
        plot_file (_type_, optional): _description_. Defaults to None.
    """
    # Plot original data
    dose_response = pd.DataFrame({"Dose": np.power(10, x), "Response": y})
    ax = sns.scatterplot(
        x="Dose", y="Response", color="red", data=dose_response, label="Response data"
    )

    # Plot fit data
    # need scaling metrics to convert skewnormal fit to actual data set
    start = np.mean(y[:3])
    end = np.mean(y[-3:])
    range = end - start

    x1 = np.linspace(
        skewnorm.ppf(0.00001, a, loc, scale), skewnorm.ppf(0.999, a, loc, scale), 100
    )
    sns.lineplot(
        x="dose_pred",
        y="response_pred",
        data=pd.DataFrame(
            {
                "dose_pred": np.power(10, x1),
                "response_pred": data_processing.restore_response(
                    skewnorm.cdf(x1, a, loc, scale), start, range
                ),
            }
        ),
        color="blue",
        lw=5,
        alpha=0.6,
        label="Skewnormal CDF fit",
    )
    ax.text(
        0.05,
        0.4,
        f"a: {a:.3f}\nloc: {loc:.3f}\nscale:{scale:.3f}",
        ha="left",
        va="bottom",
        transform=ax.transAxes,
    )
    ax.set(title=f"Skewnormal curve fit")

    plt.xscale("log")
    if plot_file:
        plt.savefig(plot_file)
    else:
        plt.show()
    plt.clf()


def plot_hill_eqn(hill_coef, ac50):
    """_summary_

    Args:
        hill_coef (_type_): _description_
        ac50 (_type_): _description_
    """
    fit_eqn_limit_lower = 0.01
    fit_eqn_limit_upper = 0.9

    ac50 = ac50

    x1 = np.linspace(
        data_processing.inv_hill_eqn(fit_eqn_limit_lower, hill_coef, ac50),
        data_processing.inv_hill_eqn(fit_eqn_limit_upper, hill_coef, ac50),
        100,
    )

    hill_data = data_processing.hill_eqn(x1, hill_coef, ac50)

    ax = sns.lineplot(
        x="dose_pred",
        y="response_pred",
        data=pd.DataFrame({"dose_pred": x1, "response_pred": hill_data}),
        color="green",
        lw=5,
        alpha=0.6,
        label="Hill equation",
    )

    ax.set(title=f"Hill equation")
    plt.xscale("log")
    plt.show()


def plot_skewnorm_hill(x, y, a, loc, scale, hill_coef, ac50, plot_file=None):
    """_summary_

    Args:
        x (_type_): _description_
        y (_type_): _description_
        a (_type_): _description_
        loc (_type_): _description_
        scale (_type_): _description_
        hill_coef (_type_): _description_
        ac50 (_type_): _description_
        plot_file (_type_, optional): _description_. Defaults to None.
    """
    # Plot original data
    dose_response = pd.DataFrame({"Dose": np.power(10, x), "Response": y})
    ax = sns.scatterplot(
        x="Dose", y="Response", color="red", data=dose_response, label="Response data"
    )

    # Plot fit data
    # need scaling metrics to convert skewnormal fit to actual data set
    start = np.mean(y[:3])
    end = np.mean(y[-3:])
    range = end - start
    # log_ac50 = np.log10(ac50)

    # fit_eqn_limit_lower = 0.00001
    # fit_eqn_limit_upper = 0.99999

    fit_eqn_limit_lower = 0.0000001
    fit_eqn_limit_lower_hill = 0.01
    fit_eqn_limit_upper = 0.9999999
    fit_eqn_limit_upper_hill = 0.999

    x1 = np.linspace(
        skewnorm.ppf(fit_eqn_limit_lower, a, loc, scale),
        skewnorm.ppf(fit_eqn_limit_upper, a, loc, scale),
        100,
    )

    x2 = np.logspace(
        np.log10(
            data_processing.inv_hill_eqn(fit_eqn_limit_lower_hill, hill_coef, ac50)
        ),
        np.log10(
            data_processing.inv_hill_eqn(fit_eqn_limit_upper_hill, hill_coef, ac50)
        ),
        100,
    )

    # x2 = np.linspace(inv_hill_eqn(fit_eqn_limit_lower_hill, hill_coef, ac50),
    #                  inv_hill_eqn(fit_eqn_limit_upper_hill, hill_coef, ac50),
    #                  100)

    # Plot Skewnormal
    ax = sns.lineplot(
        x="dose_pred",
        y="response_pred",
        data=pd.DataFrame(
            {
                "dose_pred": np.power(10, x1),
                "response_pred": data_processing.restore_response(
                    skewnorm.cdf(x1, a, loc, scale), start, range
                ),
            }
        ),
        color="blue",
        lw=5,
        alpha=0.6,
        label="Skewnormal CDF fit",
    )

    hill_data = data_processing.hill_eqn(x2, hill_coef, ac50)
    transformed_hill_data = data_processing.restore_response(hill_data, start, range)
    sns.lineplot(
        x="dose_pred",
        y="response_pred",
        data=pd.DataFrame({"dose_pred": x2, "response_pred": transformed_hill_data}),
        color="green",
        lw=5,
        alpha=0.6,
        label="Hill equation fit",
    )
    ax.set(title=f"Dose response curve fits ")
    plt.xscale("log")
    if plot_file:
        plt.savefig(plot_file)
    else:
        plt.show()
    plt.clf()
