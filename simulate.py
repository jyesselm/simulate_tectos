import argparse
import rnamake
import rnamake.resource_manager as rm
import rnamake.motif_tree as motif_tree
import rnamake.secondary_structure_factory as ssfactory
import rnamake.motif_tree_topology as motif_tree_topology
import rnamake.motif_state_ensemble_tree as motif_state_ensemble_tree
import rnamake.thermo_fluc_sampler as thermo_fluc_sampler
import settings


rm.manager.add_motif(settings.RESOURCE_DIR+"/GGAA_tetraloop")
rm.manager.add_motif(settings.RESOURCE_DIR+"/GAAA_tetraloop")

def parse_args():
    parser = argparse.ArgumentParser(description='')

    parser.add_argument('-s',   help='num of steps',
                                required=False)
    parser.add_argument('-cseq',help='chip sequence',
                                default='CTAGGATATGGAAGATCCTCGGGAACGAGGATCTTCCTAAGTCCTAG',
                                required=False)
    parser.add_argument('-css', help='chip structure',
                                default='(((((((..((((((((((((....))))))))))))...)))))))',
                                required=False)
    parser.add_argument('-fseq',help='flow peice sequence',
                                default='CTAGGAATCTGGAAGTACCGAGGAAACTCGGTACTTCCTGTGTCCTAG',
                                required=False)
    parser.add_argument('-fss', help='flow peice stucture',
                                default='((((((....((((((((((((....))))))))))))....))))))',
                                required=False)
    parser.add_argument('-pdbs',help='dont run simulation just produce pdbs',
                                default=0,
                                required=False)
    args = parser.parse_args()
    return args


def old_build():
    flow_steps = "GC=AU AU=AU AU=GC GC=UA UA=AU AU=CG CG=CG CG=GC GC=AU AU=GC".split()
    chip_steps = "GC=AU AU=AU AU=GC GC=AU AU=UA UA=CG CG=CG CG=UA UA=CG CG=GC".split()

    mt = motif_tree.MotifTree()
    mt.add_motif(rm.manager.get_motif(name="GC=GC"))
    mt.add_motif(rm.manager.get_motif(name="GGAA_tetraloop", end_name="A14-A15"))
    mt.add_motif(rm.manager.get_motif(name=flow_steps[0]), parent_end_name="A7-A22")
    for step in flow_steps[1:]:
        mt.add_motif(rm.manager.get_motif(name=step))
    mt.add_motif(rm.manager.get_motif(name="GAAA_tetraloop", end_name="A149-A154"))
    mt.add_motif(rm.manager.get_motif(name=chip_steps[0]), parent_end_name="A222-A251")
    for step in chip_steps[1:]:
        mt.add_motif(rm.manager.get_motif(name=step))
    mt.write_pdbs("orgs")


def new_build():
    seq1 = 'CTAGGATATGGGGGGUUUUUGGGAACAAAAACCCCCCTAAGTCCTAG'
    seq2 = 'CTAGGAATCTGGAAGTACCGAGGAAACTCGGTACTTCCTGTGTCCTAG'

    db1  = '(((((((..((((((((((((....))))))))))))...)))))))'
    db2  = '((((((....((((((((((((....))))))))))))....))))))'

    ss = ssfactory.factory.get_structure(seq1 +"+" + seq2, db1 + "+" + db2, to_RNA=1)
    m1 = rm.manager.get_motif(name='GAAA_tetraloop')
    m2 = rm.manager.get_motif(name='GGAA_tetraloop')
    ss.add_motif(m1.secondary_structure, m1.name)
    ss.add_motif(m2.secondary_structure, m2.name)
    ss_m = ss.motif('GGAA_tetraloop')
    last_end = None
    for i, end in enumerate(ss_m.ends):
        if ss_m.end_ids[i] == 'GGGAAC_LUUUUR_CCUGUGUC_LLULUULL_GAAUCUGG_RRUURURR':
            last_end = end
            break
    conn = ss.motif_topology_from_end(ss.ends[1], last_end=last_end)
    mtt = motif_tree_topology.MotifTreeTopology(conn)
    mt = motif_tree.motif_tree_from_topology(mtt, sterics=0)
    mt.write_pdbs()


class TectoSimulation(rnamake.base.Base):
    def __init__(self):
        self.setup_options_and_constraints()
        self.sampler = None

    def setup_options_and_constraints(self):
        options = { 'temperature' : 298.15,
                    'steps'       : 1000}

        self.options = rnamake.option.Options(options)
        self.constraints = {}

    def setup(self, cmd_args, **options):
        #trim off basepairs on the end do not need them in the simulation
        ss = ssfactory.factory.get_structure(
            cmd_args.cseq[3:-3] + "+" + cmd_args.fseq[3:-3],
            cmd_args.css[3:-3]  + "+" + cmd_args.fss[3:-3],
            to_RNA=1
        )

        #update secondary structure to include tertiary contact motifs
        m1 = rm.manager.get_motif(name='GAAA_tetraloop')
        m2 = rm.manager.get_motif(name='GGAA_tetraloop')
        ss.add_motif(m1.secondary_structure, m1.name)
        ss.add_motif(m2.secondary_structure, m2.name)

        m2 = rm.manager.get_motif(name='GGAA_tetraloop')
        last_end_id = None
        for i, end in enumerate(m2.ends):
            if end.name() == "A1-A6":
                last_end_id = m2.end_ids[i]
                break
        last_end = ss.motif('GGAA_tetraloop').get_end_by_id(last_end_id)
        conn = ss.motif_topology_from_end(ss.ends[1], last_end=last_end)
        mtt = motif_tree_topology.MotifTreeTopology(conn)
        mt = motif_tree.motif_tree_from_topology(mtt, sterics=0)
        mset = motif_state_ensemble_tree.MotifStateEnsembleTree(mt=mt)

        ni = -1
        ei = -1
        for n in mt:
            if n.data.name == 'GGAA_tetraloop':
                ni = n.index
                bp = n.data.get_basepair(name="A1-A6")[0]
                ei = n.data.ends.index(bp)


        if cmd_args.pdbs:
            mt.write_pdbs()
            mst = mset.to_mst()
            n1 = mst.get_node(ni)
            n2 = mst.get_node(mt.last_node().index-1)
            end_state_1 = n1.data.cur_state.end_states[ei]
            end_state_2 = n2.data.cur_state.end_states[1]
            scorer = thermo_fluc_sampler.FrameScorer()
            score = scorer.score(end_state_1, end_state_2)
            print "distance + rotation diff = ",score
            exit()

        self.simulator = thermo_fluc_sampler.ThermoFlucSimulation()
        self.simulator.setup(mset, ni, mt.last_node().index-1, ei, 1)
        self.simulator.run()





if __name__ == '__main__':
    args = parse_args()


    sim = TectoSimulation()
    sim.setup(args)
