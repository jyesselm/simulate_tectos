import pandas as pd
import copy
import rnamake.secondary_structure_factory as ssf

def get_steps(df):
    steps = {}
    ref_ss = '(((((((..((((((((((((....))))))))))))...)))))))'
    for i, row in df.iterrows():
        seq = row['sequence'][11:19]+"+"+row['sequence'][27:35]
        db  = ref_ss[11:19]+"+"+ref_ss[27:35]

        ss = ssf.factory.get_structure(seq, db)

        bp_steps = ss.motifs("BP_STEP")
        for bp in bp_steps:
            steps[bp.sequence()] = 0

    return steps


def basepair_step_dependence(df):

    steps = get_steps(df)
    df['abs_diff'] = abs(df["dG_normalized"] - df["dG_predicted"])
    df = df.sort(['abs_diff'])
    ref_ss = '(((((((..((((((((((((....))))))))))))...)))))))'

    sorted_steps = steps.keys()
    sorted_steps.sort()

    new_df = pd.DataFrame()
    new_df['abs_diff'] = []
    for k in sorted_steps:
        new_df[k] = []

    loc = 0

    for i, row in df.iterrows():
        seq = row['sequence'][11:19]+"+"+row['sequence'][27:35]
        db  = ref_ss[11:19]+"+"+ref_ss[27:35]

        current_steps = copy.deepcopy(steps)

        ss = ssf.factory.get_structure(seq, db)
        d = {}
        bp_steps = ss.motifs("BP_STEP")
        for bp in bp_steps:
            current_steps[bp.sequence()] += 1

        sorted_values = [row["abs_diff"]]
        for k in sorted_steps:
            sorted_values.append(current_steps[k])

        new_df.loc[loc] = sorted_values
        loc += 1

    return new_df