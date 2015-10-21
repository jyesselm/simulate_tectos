import argparse
import seaborn as sns
import pandas as pd

import plot
import datasets
import report

sns.set_style("white")
sns.set_context("talk")

def parse_args():
    parser = argparse.ArgumentParser(description='')

    parser.add_argument('-plots',   help='list of plots to include',
                                    default = None,
                                    required=False)
    parser.add_argument('-pf',      help='file with list of plots',
                                    default = None,
                                    required=False)
    parser.add_argument('-f',       help="file for flags",
                                    default=None,
                                    required=False)
    parser.add_argument('-rn',      help="preset report name",
                                    default=None,
                                    required=False)

    args = parser.parse_args()

    if args.f is not None:
        print "made it"

    return args

args = parse_args()

df = pd.read_table("summary.txt")

bp_df = datasets.basepair_step_position_dependence(df)
plotter = plot.TectoBasepairStepPositionCorrelation()
plotter.plot(bp_df)

sns.plt.savefig("test.png")

exit()

r = report.TectoBasicReport("test/")
r.get_html(df)