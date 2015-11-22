import argparse
import subprocess
import pandas as pd
import math
import numpy as np
import seaborn as sns
import os
from scipy import stats
from rnamake import base, option, vienna

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
    frac = float(lowest['avg_hit_count']) / float((1000000 - lowest['avg_hit_count']))
    lowest_log = math.log(frac)
    for i, r in df.iterrows():
        p = 1.9872041e-3*298*(lowest_log - math.log(float(r['avg_hit_count'])/float(1000000 - r['avg_hit_count'])))
        #p1 = 1.9872041e-3*298*math.log(float(lowest['avg_hit_count'])/float(r['avg_hit_count']))
        act = r['dG'] - lowest['dG']
        predicted_ddGs.append(p)
        actual_ddGs.append(act)

    df['predicted_ddG'] = predicted_ddGs
    df['actual_ddG']    = actual_ddGs

    return df


def plot_results(df ):
    print "r^2",r2(df['predicted_ddG'], df['actual_ddG'])
    print "abs error", abs_error(df['predicted_ddG'], df['actual_ddG'])
    sns.lmplot(x="actual_ddG", y="predicted_ddG", data=df)


class SimulateTectosWrapper(base.Base):
    def __init__(self):
        self.setup_options_and_constraints()

        if not os.environ.get('RNAMAKE'):
            raise ValueError("you might want to set RNAMAKE in your .bashrc")

        self.path = os.environ.get('RNAMAKE') + \
                    "/rnamake/lib/RNAMake/cmake/build/simulate_tectos"

        if not os.path.isfile(self.path):
            raise  ValueError("simulate_tectos program does not exist please compile c++ code")

    def setup_options_and_constraints(self):
        options = { 'cseq'      : "",
                    'css'       : "",
                    'fseq'      : "",
                    'fss'       : "",
                    'extra_mse' : "",
                    'n'         : 1}

        self.exec_options = { o : 1 for o in "cseq css fseq fss extra_mse".split()}
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


class SimulateTectosWrapperRun(base.Base):
    def __init__(self):
        self.setup_options_and_constraints()

    def setup_options_and_constraints(self):
        options = { 'devel'    : 0,
                    'n'        : 1,
                    'extra_mse': "",
                    'out_file' : 'results.csv'}

        self.options = option.Options(options)

    def run(self, df, **options):
        self.options.dict_set(options)
        v = vienna.Vienna()

        if 'sequence' not in df:
            raise ValueError("df must include sequence as a column")

        avg_hit_counts = []
        for i,r in df.iterrows():
            stw = SimulateTectosWrapper()
            seq = r['sequence']
            if 'secondary_structure' not in df:
                ss = v.fold(seq).structure
            else:
                ss = r['secondary_structure']
            avg_hit_count = stw.run(n=self.option('n'), extra_mse=self.option('extra_mse'),
                                    cseq=seq, css=ss)
            avg_hit_counts.append(avg_hit_count)

        df['avg_hit_count'] = avg_hit_counts
        df.to_csv(self.option('out_file'))
