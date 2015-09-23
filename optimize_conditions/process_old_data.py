import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
from scipy import stats
import math
import copy
import plot
import datasets


import rnamake.secondary_structure_factory as ssf


sns.set_style("white")
sns.set_context("talk")

def r2(x, y):
    return stats.pearsonr(x, y)[0] ** 2

df = pd.read_table("summary_org.txt", sep=" ",
                   names="sequence dG dG_lb dG_ub cutoff_count".split())

def get_summary_df():

    exp_data = pd.read_table("flowWC.150607_150605.combined.results.sequence.length10.dat",
                             names="id,sequence,dG,eminus,eplus".split(","))
    pred_data = pd.read_table("data/round_two_raw.dat", sep=" ")

    df = pd.merge(exp_data, pred_data, on="sequence", how="outer")
    df = df[df.apply(lambda x:  not pd.isnull(x['dG_predicted']), axis=1)]

    lowest = df[df['sequence'] == 'CTAGGATATGGGGAGGGTTGGGGAACCAACCCTCCCCTAAGTCCTAG']
    lowest_row = None
    for i, r in lowest.iterrows():
        lowest_row = r

    df['dG_normalized'] = df['dG'] - lowest_row['dG']

    return df


#df = sns.load_dataset("anscombe")
#print df

df = get_summary_df()
bp_df = datasets.basepair_step_dependence(df)

plotter = plot.TectoOverallCompare()
plotter.plot(df)

plt.show()

exit()

new_df = pd.DataFrame(columns=["dataset", "step appearance", "abs(dG - dG_predicted)"])
pos = 0

abs_diff = bp_df["abs_diff"]

count = 0
for i, c in bp_df.iteritems():
    if i == "abs_diff":
        continue
    print i, round(r2(abs_diff, c),2)
    for j, e in enumerate(c):
        new_df.loc[pos] = [i, abs_diff[j], e]
        pos += 1

sns.lmplot(x="step appearance", y="abs(dG - dG_predicted)",
           col="dataset", hue="dataset", data=new_df,
           ci=None, palette="muted", size=4, col_wrap=4,
           scatter_kws={"s": 50, "alpha": 1})

plt.show()


exit()



sequences = df['sequence']
unique = { seq : 1 for seq in sequences }

avg_hit_count = []
trial_num = []
seqs = []
dG = []
for seq in unique.keys():
    df_sub = df.loc[df['sequence'] == seq]

    avg_count = 0
    for j, r in df_sub.iterrows():
        avg_count += r['cutoff_count']
    avg_count /= len(df_sub)
    avg_hit_count.append(avg_count)
    trial_num.append(len(df_sub))
    seqs.append(seq)
    dG.append(r['dG'])

df = pd.DataFrame(seqs, columns=['sequence'])
df['dG'] = dG
df['avg_hit_count'] = avg_hit_count
df['trial_num'] = trial_num

lowest = None
for i, row in df.iterrows():
    if lowest is None:
        lowest = row
        continue
    if lowest['dG'] < row['dG']:
        lowest = row

df['dG_normalized'] = df['dG'] - lowest['dG']
dG_predicted = []
for i, row in df.iterrows():
    try:
        prediction = 1.9872041e-3*298*math.log(float(lowest['avg_hit_count'])/float(row['avg_hit_count']))
        dG_predicted.append(prediction)
    except:
        dG_predicted.append('NaN')

df['dG_predicted'] = dG_predicted





exit()




exit()
