import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import gridspec
from scipy import stats
import math


def r2(x, y):
    return stats.pearsonr(x, y)[0] ** 2

class TectoOverallCompare(object):
    def __init__(self):
        pass

    def plot(self, df):

        fig = plt.figure()

        gs = gridspec.GridSpec(2,2)
        ax = fig.add_subplot(gs[0])
        ax.set_title("correlation")
        sns.regplot(x="dG_normalized", y="dG_predicted",data=df, ax=ax)
        r_correlation = r2(df["dG_normalized"], df["dG_predicted"])
        ax.annotate('r2 = ' + str(round(r_correlation, 2)), xy=(1, 1),
            xycoords='axes fraction', fontsize=16,
            horizontalalignment='right', verticalalignment='top')

        ax = fig.add_subplot(gs[1])
        sns.distplot(df["dG_normalized"] - df["dG_predicted"], ax=ax)

        ax = fig.add_subplot(gs[2])

        df = df.sort(['dG'])

        sns.barplot(df['dG_normalized'], df["dG_normalized"] - df["dG_predicted"], ax=ax,                              linewidth=0.1)
        ax.xaxis.set_major_formatter(plt.NullFormatter())
        ax.set_xlabel('sorted dGs')
        ax.set_ylabel("error in prediction")

        sns.plt.tight_layout(w_pad=0)


class TectoBasepairStepHistrogram(object):
    def __init__(self):
        pass

    def plot(self, df):

        fig, ax = plt.subplots(figsize=(4,10))

        x_axis = list(df.columns.values)
        sns.heatmap(df[1:], cmap="coolwarm")
        ax.set_xticklabels(x_axis, rotation='vertical')
        #ax.set_yticklabels([round(x,2) for x in df['abs_diff'][::-1]])
        dummy = [i for i in range(len(df))]
        #plt.setp(plt.gca().get_ymajorticklabels(),
        #        size=8, rotation=-360)
        plt.setp(plt.gca().get_xmajorticklabels(),
                size=7)


class TectoBasepairStepCorrelation(object):
    def __init__(self):
        pass

    def plot(self, df):
        