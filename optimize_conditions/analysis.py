import sys
import pandas as pd
import seaborn as sns
from scipy import stats
import numpy as np
import rnamake.secondary_structure_factory as ssf
import rnamake.basic_io as basic_io

sns.set_style("white")
sns.set_context("talk")

def totuple(a):
    try:
        return tuple(totuple(i) for i in a)
    except TypeError:
        return a

def r2(x, y):
    return stats.pearsonr(x, y)[0] ** 2

df = pd.read_table("summary.txt",
                   names="sequence dG dG_lb dG_ub cutoff_count score d1 d2 r1 r2".split())
df = df[df.apply(lambda x:  not pd.isnull(x['dG']), axis=1)]

d1, d2, r1, r2 = [], [], [], []

for i,r in df.iterrows():
    spl = r['cutoff_count'].split(",")
    df.ix[i, 'score'] = spl[0]
    d1.append(basic_io.str_to_point(spl[1]))
    d2.append(basic_io.str_to_point(spl[2]))
    r1.append(basic_io.str_to_matrix(spl[3]))
    r2.append(basic_io.str_to_matrix(spl[4]))

df.drop('cutoff_count', 1)
df['d1'] = d1
df['d2'] = d2
df['r1'] = r1
df['r2'] = r2

df.to_csv("data_dump_lowest_energy_instead_of_ensemble.csv")

exit()


df = pd.read_table("summary.txt", sep=" ")
df = df[df.apply(lambda x:  not pd.isnull(x['dG']), axis=1)]
df = df[df.apply(lambda x:  x['cutoff_count'] < 23, axis=1)]


print r2(df['dG'], df['cutoff_count'])

sns.lmplot(x='dG', y='cutoff_count',data=df)
sns.plt.ylabel("score")
sns.plt.show()


"""constructs = process.get_constructs_from_file(sys.argv[1])
db = "(((((((..((((((((((((....))))))))))))...)))))))"

for c in constructs:
    ss = ssf.factory.get_structure(c.seq, db)
    if c.diff_abs() > 0.1:
        continue
    if c.dG_lower > 0.3 or c.dG_upper > 0.3:
        continue
    d = {}
    bp_steps = ss.motifs("BP_STEP")
    for bp in bp_steps:
        if bp.sequence() not in d:
            d[bp.sequence()] = 0
        d[bp.sequence()] += 1

    items = d.items()
    items.sort(key=lambda x : x[1])
    print c.diff(), c.dG_lower, c.dG_upper, c.seq"""