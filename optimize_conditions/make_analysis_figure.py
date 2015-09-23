import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
from scipy import stats

import rnamake.secondary_structure_factory as ssf

sns.set_style("white")
sns.set_context("talk")

def r2(x, y):
    return stats.pearsonr(x, y)[0] ** 2

df = pd.read_table("formated_summary.txt")
df = df[df.apply(lambda x:  not pd.isnull(x['dG_normalized']) and x['trial_num'] == 3, axis=1)]

steps =[]
error = []
avg = 0
count = 0

for i, row in df.iterrows():
    seq = row['sequence'][13:19]+"+"+row['sequence'][27:33]
    db  = row['ss'][13:19]+"+"+row['ss'][27:33]
    avg +=  row["dG_normalized"], row["dG_predicted"]
    count += 1

    ss = ssf.factory.get_structure(seq, db)
    d = {}
    bp_steps = ss.motifs("BP_STEP")
    for bp in bp_steps:
        if bp.sequence() not in d:
            d[bp.sequence()] = 0
        d[bp.sequence()] += 1

    items = d.items()
    items.sort(key=lambda x : x[1])
    success = 0
    for v in items:
        if v[1]> 2:
            success = 1
    if success != 1:
        continue
    print row['sequence'], items, row["dG_normalized"], row["dG_predicted"]


fig = plt.figure()
gs = gridspec.GridSpec(2,2)

ax = fig.add_subplot(gs[0])
sns.regplot(x="dG_normalized", y="dG_predicted",data=df, ax=ax)
r_correlation = r2(df["dG_normalized"], df["dG_predicted"])
ax.annotate('r2 = ' + str(round(r_correlation, 2)), xy=(1, 1),
            xycoords='axes fraction', fontsize=16,
            horizontalalignment='right', verticalalignment='top')

ax = fig.add_subplot(gs[1])
sns.distplot(df["dG_normalized"] - df["dG_predicted"], ax=ax)


sns.plt.show()
