from rnamake import sqlite_library, motif_ensemble, util, motif_factory
from rnamake import secondary_structure_factory as ssf
from rnamake import resource_manager as rm
import os
import pandas as pd


def build_mse_by_topology(m_name, sequence_identity=1):
    """
    Builds an motif ensemble to be use in the tecto simulations using the junction
    topology of a given motif to find all motifs that have the, ignores the flanking
    basepairs

    i.e. if you give a motif name with seq = GAC&GC ss = (.(&)), this code will build an
    ensemble of all A 1x0 bulges in the database.

    .. code-block:: python
        >>> build_mse_by_topology('TWOWAY.1DUQ.7')
        None
        #these files will be generated
        #extra_mse/TWOWAY.1DUQ.7-G408-H422.me
        #extra_mse/TWOWAY.1DUQ.7-G409-H420.me

    :param m_name: The name of the motif that you need an ensemble for
    :type  m_name: string
    :return: None, generates files ensemble files in /extra_mse/

    :param sequence_identity: selects only motifs that have the same junction sequence
        other then the flanking basepairs. default : 1, turned on
    :type  sequence_identity: int

    """

    mlib = sqlite_library.MotifSqliteLibrary("twoway")
    mlib.load_all()

    seeds = mlib.get_multi(name=m_name)

    for s in seeds:
        motifs = []
        ss1 = ssf.factory.secondary_structure_from_motif(s)
        res1 = ss1.residues()
        for m in mlib.all():
            if len(m.residues()) != len(s.residues()):
                continue
            ss2 = ssf.factory.secondary_structure_from_motif(m)
            res2 = ss2.residues()
            fail = 0
            for i in range(len(res1)):
                if res1[i].dot_bracket != res2[i].dot_bracket:
                    fail = 1
                    break
                if res1[i].dot_bracket != ".":
                    continue
                if res1[i].name != res2[i].name and sequence_identity == 1:
                    fail = 1
                    break
            if not fail:
                motifs.append(m)

        energies = [1 for m in motifs]
        me = motif_ensemble.MotifEnsemble()
        me.setup(motifs[0].end_ids[0], motifs, energies)
        str = me.to_str()
        for m in motifs:
            m.to_pdb(m.name + ".pdb")
        #return
        f = open("extra_mse/" + s.name + "-" + s.ends[0].name() + ".me", "w")
        f.write(str)
        f.close()


def build_bp_step_mse_for_topology(m_name):
    """
    Creates an ensemble using a basepair distribution of the first and last
    basepair of a motif. Not recommended to do for anything large then a 1-0 bulge

    .. code-block:: python
        >>> build_bp_step_mse_for_topology('TWOWAY.1DUQ.7')
        None
        #these files will be generated
        #extra_mse/TWOWAY.1DUQ.7-G408-H422.me
        #extra_mse/TWOWAY.1DUQ.7-G409-H420.me

    :param m_name: The name of the motif that you need an ensemble for
    :type  m_name: string
    :return: None, generates files ensemble files in /extra_mse/
    """

    mlib = sqlite_library.MotifSqliteLibrary("twoway")

    seeds = mlib.get_multi(name=m_name)

    for s in seeds:
        end_id  = s.ends[0].res1.name + s.ends[1].res2.name + "_LL_"
        end_id += s.ends[1].res1.name + s.ends[0].res2.name + "_RR"
        me = rm.manager.get_motif_ensemble(end_id)
        str = me.to_str()
        f = open("extra_mse/" + s.name + "-" + s.ends[0].name() + ".me", "w")
        f.write(str)
        f.close()

def build_mse_from_motif_for_topology(m_name, actual_name, mlib=None):
    if mlib is None:
        mlib = sqlite_library.MotifSqliteLibrary("twoway")
    act_m = mlib.get(name=actual_name)
    m = mlib.get(name=m_name)
    for i, end in enumerate(m.ends):
        mi = mlib.get(name=m_name, end_name=end.name())
        mi.to_pdb("m."+str(i)+".pdb")
        me = motif_ensemble.MotifEnsemble()
        all_ms = [mi]
        all_scores = [1]
        mi.name = actual_name + "-" + str(i)
        me.setup(mi.end_ids[0], all_ms, all_scores)
        print mi.end_ids[0], act_m.end_ids[i], act_m.ends[i].name()
        f = open("extra_mse/" + actual_name + "-" + act_m.ends[i].name() + ".me", "w")
        f.write(me.to_str())
        f.close()

def build_mse_from_motif_for_topology_2(org_m, new_m, mlib):

    new_m_key = new_m.name + "-" + new_m.ends[0].name()

    fsum = open("extra_mse/" + org_m.name + "/" + new_m_key + "/extra_mse.dat", "w")

    pos = os.path.abspath(__file__)
    basedir = util.base_dir(pos)

    for i, end in enumerate(org_m.ends):
        mi = mlib.get(name=new_m.name, end_name=new_m.ends[i].name())
        all_ms = [mi]
        all_scores = [1]
        me = motif_ensemble.MotifEnsemble()
        me.setup(org_m.end_ids[i], all_ms, all_scores)

        org_m_key = org_m.name + "-" + end.name()

        path = "extra_mse/" + org_m.name + "/" + new_m_key + "/" + org_m_key + ".me"
        f = open(path, "w")
        f.write(me.to_str())
        f.close()

        fsum.write(org_m_key + " " + basedir+path + "\n")
    fsum.close()


def get_extra_mse_file():
    """
    generates the extra_mse file required for custom ensembles in simulate_tectos
    program 

    """
    pos = os.path.abspath(__file__)
    basedir = util.base_dir(pos)
    f = open("extra_mse.dat", "w")
    for f_name in os.listdir(basedir + "/extra_mse/"):
        f.write(f_name[:-3] + " " + basedir + "extra_mse/"+f_name + "\n")
    f.close()


mlib = sqlite_library.MotifSqliteLibrary("twoway")
mlib.load_all()

df = pd.read_csv("../simulate/data_all.csv")
names = df.m_name.unique()


for m1 in mlib.all():
    if m1.name not in names:
        continue

    try:
        os.mkdir("extra_mse/" + m1.name )
    except:
        continue

    for i, m2 in enumerate(mlib.all()):

        os.mkdir("extra_mse/" + m1.name + "/" + m2.name + "-" + m2.ends[0].name())

        build_mse_from_motif_for_topology_2(m1, m2, mlib)
        if i > 4:
            break

#build_mse_from_motif_for_topology("TWOWAY.1S72.68", "TWOWAY.1DUQ.7")
#build_mse_from_motif_for_topology("TWOWAY.1DUQ.7", "TWOWAY.1DUQ.7")
#build_mse_by_topology("TWOWAY.1DUQ.7", sequence_identity=1)
#m = rm.manager.get_motif(name="GC=CG").to_pdb("step.pdb")
#build_bp_step_mse_for_topology("TWOWAY.1DUQ.7")
#get_extra_mse_file()

