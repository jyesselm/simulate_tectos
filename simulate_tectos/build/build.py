from rnamake import sqlite_library, motif_ensemble, util, motif_factory
from rnamake import secondary_structure_factory as ssf
from rnamake import resource_manager as rm
import os


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
        exit()
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


build_mse_by_topology("TWOWAY.1DUQ.7", sequence_identity=0)
#build_bp_step_mse_for_topology("TWOWAY.1DUQ.7")
#get_extra_mse_file()

