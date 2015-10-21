import subprocess
import pandas as pd
import math

prog = "/Users/josephyesselman/projects/RNAMake/rnamake/lib/RNAMake/cmake/build/simulate_tectos"

def reformat_data_file():
    f = open("flowWC.150607_150605.combined.results.sequence.length10.dat")
    lines = f.readlines()
    f.close()

    lines.pop(0)
    data = []
    for l in lines:
        spl = l.split()
        if len(spl) < 5:
            continue
        if float(spl[4]) < 1 and float(spl[3]) < 1:
            data.append([spl[1], spl[2], spl[3], spl[4]])

    data.sort(key= lambda x : x[1], reverse=True)

    f = open("exhustive_helices.dat", "w")

    for d in data:
        f.write(d[0] + " " + d[1] + " " + d[2] + " " + d[3] + "\n")
    f.close()

def run_new():
    f = open("summary.txt", "w")
    #df = pd.read_table("sequences_tiled_length10.txt")
    df = pd.read_table("flowWC.150607_150605.combined.results.sequence.length10.dat",
                          names="id,sequence,dG,eminus,eplus".split(","))
    f.write("sequence dG dG_lb dG_ub cutoff_count\n")
    for i, row in df.iterrows():
        if row['dG'] == 'nan':
            continue
        for j in range(1):
            try:
                output = subprocess.check_output(prog + " -cseq " + row['sequence'], shell=True)
            except:
                continue
            print row['sequence'], row['dG'], row['eminus'], row['eplus'], output,
            f.write(row['sequence'] + "\t" + str(row['dG']) + "\t" + str(row['eminus']) + "\t" + str(row['eplus']) + "\t" + output)
            f.flush()
    f.close()

def get_test_set():
    f  = open("exhustive_helices.dat")
    lines = f.readlines()
    f.close()

    f = open("test_set.dat", "w")

    i = 0
    while i < len(lines) :
        spl = lines[i].split()
        if float(spl[2])  > 0.1 or float(spl[3]) > 0.1:
            i += 1
            continue
        f.write(lines[i])
        i += 19
    f.close()

def fix_old_data():
    df = pd.read_table("sequences_tiled_length10.txt")
    df_new = pd.read_table("summary.txt",sep=" ")

    avg_hit_count = []
    trial_num = []

    for i, row in df.iterrows():
        rows = df_new.loc[df_new['sequence'] == row['sequence']]
        if len(rows) == 0:
            avg_hit_count.append(0)
            trial_num.append(0)
            continue

    avg_count = 0
    for j, r in rows.iterrows():
        avg_count += r['cutoff_count']
        avg_count /= len(rows)
        avg_hit_count.append(avg_count)
        trial_num.append(len(rows))

    df['avg_hit_count'] = avg_hit_count
    df['trial_num'] = trial_num

    lowest = None
    for i, row in df.iterrows():
        if lowest is None:
            lowest = row
            continue
        if lowest['dG'] < row['dG']:
            lowest = row

    #dG_prediction = 1.9872041e-3*298*math.log(lowest[3] / v[3])

    df['dG_normalized'] = df['dG'] - lowest['dG']
    dG_predicted = []
    for i, row in df.iterrows():
        try:
            prediction = 1.9872041e-3*298*math.log(float(lowest['avg_hit_count'])/float(row['avg_hit_count']))
            dG_predicted.append(prediction)
        except:
            dG_predicted.append('nan')

    df['dG_predicted'] = dG_predicted

    df.to_csv("formated_summary.txt", sep="\t")

def get_predicteted_dG(row, lowest_count):
    return 1.9872041e-3*298*math.log(float(lowest_count)/float(row['avg_hit_count']))


run_new()
exit()

"""

df = pd.read_table("test_set.dat", sep=" ")
df['avg_hit_count'] = [0 for x in range(len(df)) ]
df['trials'] = [0 for x in range(len(df)) ]


for i ,r in df.iterrows():
    outputs = []
    for i in range(3):
        try:
            output = subprocess.check_output(prog + " -cseq " + r['sequence'], shell=True)
            outputs.append(output.rstrip())
        except:
            pass
    outputs = [float(x) for x in outputs]
    avg = 0
    for o in outputs:
        avg += o
    avg /= len(outputs)
    df.ix[i, 'avg_hit_count'] = avg
    df.ix[i, 'trials'] = len(outputs)
    print r['sequence'], r['dG'], avg, len(outputs)

df.to_csv("summary.txt", sep="\t")
"""

exp_data = pd.read_table("flowWC.150607_150605.combined.results.sequence.length10.dat",
                          names="id,sequence,dG,eminus,eplus".split(","))
pred_data = pd.read_table("summary_2.txt", sep=" ")
pred_data = pred_data.drop('dG', 1)

df = pd.merge(exp_data, pred_data, on="sequence", how="outer")
df = df[df.apply(lambda x:  not pd.isnull(x['avg_hit_count']), axis=1)]

lowest = df.ix[df['dG'].idxmax()]
df['dG_normalized'] = df['dG'] - lowest['dG']
df['dG_predicted']  = df.apply (lambda row: get_predicteted_dG(row,lowest['avg_hit_count'] ),axis=1)

df.to_csv("summary.txt", sep="\t")














