import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
from scipy import stats

import plot
import datasets

sns.set_style("white")
sns.set_context("talk")

df = pd.read_table("summary.txt")

p1 = plot.TectoOverallCompare()
p1.plot(df)

plt.savefig("test_1.pdf", format='pdf')

bp_df = datasets.basepair_step_dependence(df)

p1 = plot.TectoBasepairStepCorrelation()
p1.plot(bp_df)

plt.savefig("test_2.pdf", format='pdf')

exit()


df = pd.read_table("data/old_prediction_data/1bp_flank_mismatch_training.txt")

print df.to_html("test.html")
exit()

sns.lmplot(x="dG_normalized", y="dG_predicted",
           col="name", hue="name", data=df,
           ci=None, palette="muted", size=4, col_wrap=2,
           scatter_kws={"s": 50, "alpha": 1})

for name, df in df.groupby('name'):
    print name, plot.r2(df['dG_normalized'], df['dG_predicted'])

plt.show()