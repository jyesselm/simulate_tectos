import subprocess
import os
from rnamake import base, option, vienna

class SimulateTectosWrapper(base.Base):
    def __init__(self):
        self.setup_options_and_constraints()

        if not os.environ.get('RNAMAKE'):
            raise ValueError("you might want to set RNAMAKE in your .bashrc")

        self.path = os.environ.get('RNAMAKE') + \
                    "/rnamake/lib/RNAMake/cmake/build/simulate_tectos_devel"

        if not os.path.isfile(self.path):
            raise  ValueError("simulate_tectos program does not exist please compile c++ code")

    def setup_options_and_constraints(self):
        options = { 'cseq'      : "",
                    'css'       : "",
                    'fseq'      : "",
                    'fss'       : "",
                    'extra_mse' : "",
                    'n'         : 1,
                    's'         : '1000000'}

        self.exec_options = { o : 1 for o in "s cseq css fseq fss extra_mse".split()}
        self.default_options = {}
        self.options = option.Options(options)
        for k in self.exec_options.keys():
            self.default_options[k] = self.option(k)

    def run(self, **options):
        self.options.dict_set(options)
        cmd = self._get_command()
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
                    'fseq'     : "",
                    'fss'      : "",
                    'extra_mse': "",
                    's'        : '1000000',
                    'out_file' : 'results.csv'}

        self.options = option.Options(options)
        self.simulate_options = {x : 1 for x in "devel,n,extra_mse,fss,fseq,s".split(",")}

    def run(self, df, **options):
        self.options.dict_set(options)

        if 'sequence' not in df:
            raise ValueError("df must include sequence as a column")

        stw = SimulateTectosWrapper()
        for k, v in options.iteritems():
            if k in self.simulate_options:
                stw.option(k, v)

        v = vienna.Vienna()

        avg_hit_counts = [-1 for i in range(len(df))]
        for i,r in df.iterrows():
            seq = r['sequence']

            extra_mse = self.option('extra_mse')
            if 'extra_mse' in df:
                extra_mse = r['extra_mse']

            if 'secondary_structure' not in df:
                ss = v.fold(seq).structure
            else:
                ss = r['secondary_structure']

            try:
                avg_hit_count = stw.run(cseq=seq, css=ss, extra_mse=extra_mse)
            except:
                avg_hit_count = -1
            avg_hit_counts[i] = avg_hit_count
            df['avg_hit_count'] = avg_hit_counts

            if i % 10 == 0:
                df.to_csv(self.option('out_file'))


        df['avg_hit_count'] = avg_hit_counts
        df.to_csv(self.option('out_file'))
