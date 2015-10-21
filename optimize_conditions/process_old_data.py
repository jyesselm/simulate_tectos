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

def assign_groups(row, bin_size=0.2):
    group = round(float(row["abs_diff"])/bin_size, 0)*bin_size
    return group

def is_allowed_name(row, allowed_names):
    pass

#df = sns.load_dataset("anscombe")
#print df

df = get_summary_df()
bp_df = datasets.basepair_step_position_dependence(df)

bp_df["bin"] = bp_df.apply (lambda row: assign_groups(row),axis=1)

df_new = pd.DataFrame(columns=("name", "abs_diff", "percent", "count"))

bins = list(bp_df.bin.unique())
pos = "p1,p2,p3,p4,p5,p6,p7".split(",")
loc = 0

allowed_names = ["p1-TG&CA", "p4-AT&AT"]

for p in pos:
    steps = bp_df[p].unique()
    for b in bins:
        df_bin = bp_df[bp_df.bin == b]
        bin_count = len(df_bin)
        for s in steps:
            count = len(df_bin[bp_df[p] == s])
            df_new.loc[loc] = [str(p)+"-"+s, b, float(count)/float(bin_count), bin_count]
            loc += 1

    for s in steps:
        last_df = df_new[df_new.name == str(p)+"-"+s]
        r2_corr = r2(last_df["abs_diff"], last_df["percent"])
        #if r2_corr > 0.7:
        #    allowed_names.append(str(p)+"-"+s)

#df_new2 = df_new
df_new2 = df_new[df_new['name'].isin(allowed_names)]


df_new2.sort(["name"])
sns.lmplot(x="abs_diff", y="percent",
               col="name", hue="name", data=df_new2,
               ci=None, palette="muted", size=4, col_wrap=1,
               scatter_kws={"s": 50, "alpha": 1})

plt.savefig("test.png")
exit()


#while <


exit()

bp_df_sub = bp_df[bp_df.apply(lambda x : x["abs_diff"] > 0 and x["abs_diff"] < .20 , axis=1)]
print bp_df_sub
print len(bp_df_sub), len(bp_df_sub[bp_df_sub.p1 == "GC&GC"])

exit()

plotter = plot.TectoOverallCompare()
plotter.plot(df)

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
