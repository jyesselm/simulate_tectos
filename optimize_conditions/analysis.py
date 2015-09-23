import sys
import process

import rnamake.secondary_structure_factory as ssf

constructs = process.get_constructs_from_file(sys.argv[1])
db = "(((((((..((((((((((((....))))))))))))...)))))))"

for c in constructs:
    ss = ssf.factory.get_structure(c.seq, db)
    if c.diff_abs() > 0.1:
        continue
    if c.dG_lower > 0.3 or c.dG_upper > 0.3:
        continue
    d = {}
    bp_steps = ss.motifs("BP_STEP")
    for bp in bp_steps:
        if bp.sequence() not in d:
            d[bp.sequence()] = 0
        d[bp.sequence()] += 1

    items = d.items()
    items.sort(key=lambda x : x[1])
    print c.diff(), c.dG_lower, c.dG_upper, c.seq