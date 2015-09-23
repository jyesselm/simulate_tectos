import sys
import math

class Construct (object):
    def __init__(self, seq, dG, dG_lower, dG_upper, dG_prediction, trials):
        self.seq, self.dG, self.dG_lower, self.dG_upper = seq, dG, dG_lower, dG_upper
        self.dG_prediction, self.trials = dG_prediction, trials

    def diff(self):
        return self.dG - self.dG_prediction

    def diff_abs(self):
        return abs(self.diff())

def get_constructs_from_file(file_name):
    constructs = []
    f = open(file_name)
    lines = f.readlines()
    f.close()

    data = {}
    for l in lines:
        spl = l.split()
        if spl[0] not in data:
            data[spl[0]] = [float(spl[1]), float(spl[2]), float(spl[3]), float(spl[4]), 1]
        else:
            data[spl[0]][3] += float(spl[4])
            data[spl[0]][4] += 1

    lowest = data[spl[0]]
    for v in data.itervalues():
        if v[0] > lowest[0]:
            lowest = v

    for v in data.itervalues():
        if v == lowest:
            continue
        v[0] -= lowest[0]
        v[3] /= v[4]

    lowest[0] = 0
    lowest[3] /= lowest[4]

    for seq, v in data.iteritems():
        dG_prediction = 1.9872041e-3*298*math.log(lowest[3] / v[3])
        c = Construct(seq, v[0], v[1], v[2], dG_prediction, v[4])
        constructs.append(c)

    return constructs

if __name__ == "__main__":

    constructs = get_constructs_from_file(sys.argv[1])

    avg_diff = 0
    avg_unsigned_diff = 0
    count = 0

    for c in constructs:
        avg_diff += c.diff()
        avg_unsigned_diff += c.diff_abs()
        print c.dG, c.dG_prediction
        count += 1

    print avg_diff / count, avg_unsigned_diff / count

