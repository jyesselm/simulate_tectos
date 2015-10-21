import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import gridspec
from scipy import stats
import math


def r2(x, y):
    return stats.pearsonr(x, y)[0] ** 2

class PlotFactory(object):
    def __init__(self):
        pass

    @staticmethod
    def get_plot(pname):
        if   pname == 'TectoOverallCompare':          return TectoOverallCompare()
        elif pname == 'TectoBasepairStepHistrogram':  return TectoBasepairStepHistrogram()
        elif pname == 'TectoBasepairStepCorrelation': return TectoBasepairStepCorrelation()
        elif pname == 'TectoBasepairStepPositionCorrelation':
            return TectoBasepairStepPositionCorrelation()
        elif pname == 'TectoBasepairStepCorrelation': return TectoExhuastiveHelixSummary()
        else:
            raise ValueError("cannot get plot with name "+pname)


class TectoOverallCompare(object):
    def __init__(self):
        self.fig_count = 1
        self.name = "TectoOverallCompare"

    def savefig(self, df, dir=None):
        self.plot(df)
        if dir is not None:
            plt.savefig(dir+"/fig_0.png")
        else:
            plt.savefig("fig_0.png")

    def plot(self, df):

        fig = plt.figure(figsize=(10,3))

        gs = gridspec.GridSpec(1,3)
        ax = fig.add_subplot(gs[0])
        ax.set_title("correlation")
        sns.regplot(x="dG_normalized", y="dG_predicted",data=df, ax=ax)
        r_correlation = r2(df["dG_normalized"], df["dG_predicted"])
        ax.annotate('r2 = ' + str(round(r_correlation, 2)), xy=(1, 1),
            xycoords='axes fraction', fontsize=16,
            horizontalalignment='right', verticalalignment='top')

        ax = fig.add_subplot(gs[1])
        sns.distplot(df["dG_normalized"] - df["dG_predicted"], ax=ax)
        ax.set_xlabel('diff in prediction')
        ax.set_ylabel("fraction")

        ax = fig.add_subplot(gs[2])
        df = df.sort(['dG'])

        sns.barplot(df['dG_normalized'], df["dG_normalized"] - df["dG_predicted"], ax=ax,                              linewidth=0.1)

        ax.xaxis.set_major_formatter(plt.NullFormatter())
        ax.set_xlabel('sorted dGs')
        ax.set_ylabel("error in prediction")

        sns.plt.tight_layout(w_pad=0.5)


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
        self.fig_count = 2
        self.name = "TectoBasepairStepCorrelation"

    def _get_data(self, bp_df):
        new_df = pd.DataFrame(columns=["dataset", "step appearance", "abs(dG - dG_predicted)"])
        pos = 0

        abs_diff = bp_df["abs_diff"]

        r_corrs = []
        datasets = []

        count = 0
        for i, c in bp_df.iteritems():
            if i == "abs_diff":
                continue
            #print i, round(r2(abs_diff, c),2)
            r_corrs.append(r2(abs_diff, c))
            datasets.append(i)

            for j, e in enumerate(c):
                new_df.loc[pos] = [i, e, abs_diff[j]]
                pos += 1

        r_df = pd.DataFrame(columns=("dataset", "r2_correlation"))
        r_df['dataset'] = datasets
        r_df['r2_correlation'] = r_corrs

        return new_df, r_df

    def savefig(self, bp_df, dir=None):
        if dir is not None:
            dir = dir + "/"
        else:
            dir = ""

        new_df, r_df = self._get_data(bp_df)

        sns.set_style("white")

        sns.lmplot(x="step appearance", y="abs(dG - dG_predicted)",
                   col="dataset", hue="dataset", data=new_df,
                   ci=None, palette="muted", size=4, col_wrap=4,
                   scatter_kws={"s": 50, "alpha": 1})

        nfig = 0
        plt.savefig(dir+"fig_"+str(nfig)+".png")
        nfig += 1

        fig = plt.figure()

        ax = sns.barplot(x="dataset", y="r2_correlation",data=r_df,color="#4878CF")
        ax.set_ylabel("r2_correlation")
        for item in ax.get_xticklabels():
            item.set_rotation(45)

        plt.savefig(dir+"fig_"+str(nfig)+".png")
        nfig += 1

    def plot(self, bp_df):

        new_df, r_df = self._get_data(bp_df)

        sns.set_style("white")

        sns.lmplot(x="step appearance", y="abs(dG - dG_predicted)",
                   col="dataset", hue="dataset", data=new_df,
                   ci=None, palette="muted", size=4, col_wrap=4,
                   scatter_kws={"s": 50, "alpha": 1})

        sns.set_style("whitegrid")

        fig = plt.figure()

        ax = sns.barplot(x="dataset", y="r2_correlation",data=r_df,color="#4878CF")
        ax.set_ylabel("r2_correlation")
        for item in ax.get_xticklabels():
            item.set_rotation(45)
        #ax.set_xticklabels(rotation=90)


class TectoBasepairStepPositionCorrelation(object):
    def __init__(self):
        pass

    def plot(self, bp_df):
        #bp_df["bin"] = bp_df.apply (lambda row: self.assign_groups(row, 0.05),axis=1)

        bp_df.sort(['abs_diff'])

        bin_size = float(len(bp_df)) / 5.0
        print bin_size

        #bins = np.linspace(bp_df.abs_diff.min(), bp_df.abs_diff.max(), 10)
        #groups = bp_df.groupby(np.digitize(bp_df.abs_diff, bins))

        df_new = pd.DataFrame(columns=("name", "abs_diff", "percent", "count"))

        pos = "p1,p2,p3,p4,p5,p6,p7".split(",")
        loc = 0
        r_loc = 0

        r_pf = pd.DataFrame(columns=("name", "r2_correlation"))

        for p in pos:
            steps = bp_df[p].unique()
            for i, df in enumerate(np.array_split(bp_df, 5)):
                bin_count = len(df)
                for s in steps:
                    count = len(df[df[p] == s])
                    df_new.loc[loc] = [str(p)+"-"+s, i, float(count)/float(bin_count), bin_count]
                    loc += 1

        for s in steps:
            last_df = df_new[df_new.name == str(p)+"-"+s]
            r2_corr = r2(last_df["abs_diff"], last_df["percent"])
            r_pf.loc[r_loc] = [str(p)+"-"+s, r2_corr]

        df_new.sort(["name"])
        sns.lmplot(x="abs_diff", y="percent",
                  col="name", hue="name", data=df_new,
                  ci=None, palette="muted", size=4, col_wrap=7,
                  scatter_kws={"s": 50, "alpha": 1})


#plots for fixed data
class TectoExhuastiveHelixSummary(object):
    def __init__(self):
        pass

    def plot(self, df):
        ax = sns.distplot(df["dG_predicted"], kde=False)
        ax.set_ylabel("Count")
        ax.set_title("dG Distribution of Helical Sequences")
