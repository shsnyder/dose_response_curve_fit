import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import skewnorm
from matplotlib import pyplot as plt


def plot_sn(
    x, a_range, loc, scale, plot_file=False, type="p", plot_title="Skewnormal Plot"
):
    data = pd.DataFrame({})
    for a in a_range:
        if type == "p":
            pdf = pd.DataFrame({a: skewnorm.pdf(x, a, loc, scale)})
            y_label = "Density"
        else:
            pdf = pd.DataFrame({a: skewnorm.cdf(x, a, loc, scale)})
            y_label = "Cumulative Probability"
        data = pd.concat([data, pdf], axis=1, ignore_index=True)

    data.columns = [str(x) for x in list(a_range)]

    sns.lineplot(data=data, lw=2).set(title=plot_title, ylabel=y_label)

    plt.legend(title="a value")

    if plot_file:
        plt.savefig(plot_file)
    else:
        plt.show()
    plt.clf()


if __name__ == "__main__":
    x = np.arange(-2, 2, 0.1)
    x = np.float32(x)

    plot_sn(
        x,
        a_range=range(-4, 1, 1),
        loc=0,
        scale=1,
        type="p",
        plot_file="/Users/shsnyder/Documents/projects/dose_response_curve_fit/sn_pd_neg_a.png",
        plot_title="Skewnormal PDF",
    )
    plot_sn(
        x,
        a_range=range(0, 5, 1),
        loc=0,
        scale=1,
        type="p",
        plot_file="/Users/shsnyder/Documents/projects/dose_response_curve_fit/sn_pd_pos_a.png",
        plot_title="Skewnormal PDF",
    )
    plot_sn(
        x,
        a_range=range(-4, 5, 1),
        loc=0,
        scale=1,
        type="p",
        plot_file="/Users/shsnyder/Documents/projects/dose_response_curve_fit/sn_pd_all_a.png",
        plot_title="Skewnormal PDF",
    )
    plot_sn(
        x,
        a_range=range(-4, 1, 1),
        loc=0,
        scale=1,
        type="c",
        plot_file="/Users/shsnyder/Documents/projects/dose_response_curve_fit/sn_cd_neg_a.png",
        plot_title="Skewnormal CDF",
    )
    plot_sn(
        x,
        a_range=range(0, 5, 1),
        loc=0,
        scale=1,
        type="c",
        plot_file="/Users/shsnyder/Documents/projects/dose_response_curve_fit/sn_cd_pos_a.png",
        plot_title="Skewnormal CDF",
    )
    plot_sn(
        x,
        a_range=range(-4, 5, 1),
        loc=0,
        scale=1,
        type="c",
        plot_file="/Users/shsnyder/Documents/projects/dose_response_curve_fit/sn_cd_all_a.png",
        plot_title="Skewnormal CDF",
    )

    plot_sn(
        x,
        a_range=[0],
        loc=0,
        scale=1,
        type="p",
        plot_file="/Users/shsnyder/Documents/projects/dose_response_curve_fit/sn_pd_zero_a.png",
        plot_title="Normal PDF",
    )
    plot_sn(
        x,
        a_range=[0],
        loc=0,
        scale=1,
        type="c",
        plot_file="/Users/shsnyder/Documents/projects/dose_response_curve_fit/sn_cd_zero_a.png",
        plot_title="Normal CDF",
    )
