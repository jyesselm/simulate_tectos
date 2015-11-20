import argparse
import subprocess
import pandas as pd
import math
import numpy as np
import seaborn as sns
from scipy import stats
from rnamake import base, option

sns.set_style("white")
sns.set_context("talk")

def r2(x, y):
    return stats.pearsonr(x, y)[0] ** 2


def abs_error(x, y):
    diff = np.abs(np.array(x) - np.array(y))
    return np.average(diff)


def parse_args():
    pass


def process_results(f_name, norm_seq=""):
    df = pd.read_csv(f_name)
    value_list = [norm_seq]

    if norm_seq == "":
        max = df['dG'].max()
        df1 = df[df.dG.isin([max])]
        lowest = df1.iloc[0]

    else:
        df1 = df[df.sequence.isin(value_list)]
        lowest = df1.iloc[0]

    predicted_ddGs = []
    actual_ddGs = []
    for i, r in df.iterrows():
        p = 1.9872041e-3*298*math.log(float(lowest['avg_hit_count'])/float(r['avg_hit_count']))
        act = r['dG'] - lowest['dG']
        predicted_ddGs.append(p)
        actual_ddGs.append(act)

    df['predicted_ddG'] = predicted_ddGs
    df['actual_ddG']    = actual_ddGs

    return df


def plot_results(df):
    print r2(df['predicted_ddG'], df['actual_ddG'])
    print abs_error(df['predicted_ddG'], df['actual_ddG'])
    sns.lmplot(x="actual_ddG", y="predicted_ddG", data=df)
    sns.plt.show()


class SimulateTectosWrapper(base.Base):
    def __init__(self):
        self.setup_options_and_constraints()
        self.path = '/Users/josephyesselman/projects/RNAMake/rnamake/lib/RNAMake/cmake/build/simulate_tectos'

    def setup_options_and_constraints(self):
        options = { 'cseq' : "",
                    'css'  : "",
                    'fseq' : "",
                    'fss'  : "",
                    'n'    : 1}

        self.exec_options = { o : 1 for o in "cseq css fseq fss".split()}
        self.default_options = {}
        self.options = option.Options(options)
        for k in self.exec_options.keys():
            self.default_options[k] = self.option(k)

    def run(self, **options):
        self.options.dict_set(options)

        avg = 0
        count = 0
        for i in range(self.option('n')):
            try:
                cmd = self._get_command()
                out = subprocess.check_output(cmd, shell=True)
                hits = int(out.rstrip())
                avg += hits
                count += 1
            except:
                pass

        return avg / count

    def _get_command(self):
        s = self.path
        for opt in self.exec_options.keys():
            val  = self.option(opt)
            dval = self.default_options[opt]
            if val == dval:
                continue
            if opt == "css" or opt == "fss":
                s += " -" + opt + " \'" + val + "\'"
            else:
                s += " -" + opt + " " + val
        return s
