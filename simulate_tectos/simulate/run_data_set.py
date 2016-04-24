import sys
import pandas as pd
import simulate_wrapper

if len(sys.argv) < 2:
    raise ValueError("did not specify .csv file")

df = pd.read_csv(sys.argv[1])

stw_run = simulate_wrapper.SimulateTectosWrapperRun()
#stw_run.run(df, extra_mse='/Users/josephyesselman/projects/RNAMake.projects/simulate_tectos/simulate_tectos/build/extra_mse.dat', n=3)
stw_run.run(df, n=3)