import src.scr3_python_colvar_generator as cg
import src.CoM_colvars as colvar
import src.unbinding as unb
import src.read_pdb as rp
import src.output as out
import argparse
import os
import sys

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--step')
    parser.add_argument('-s', "--stride", default=5, type=int)
    parser.add_argument('-l', "--lig", default="LIG")
    parser.add_argument("--top", default="find")
    parser.add_argument('-c', "--cutoff", default=3.5, type=float, help="Initial cutoff for identifying neighbours in A")
    parser.add_argument('-m', "--maxdist", default=9, type=float, help="distance for exclusion fo pairs in A")
    parser.add_argument('-t', "--trajectory")
    parser.add_argument("--cumulative", action='store_true', default=False,
                        help='Reprocess all the trajectories up to the one specified by "-t". ')
    parser.add_argument('-ns', default=10, type=int)
    parser.add_argument('--writeDCD', action='store_true', default=False, help='write the strided DCD')
    parser.add_argument('-p', "--processonly", action='store_true', default=False,
                        help='Do not write VMD and NAMD input. Other outputs will be written.')
    parser.add_argument("--nosave", action='store_true', default=False,
                        help='Do not save the checkpoint. For debug only.')
    args = parser.parse_args()
    if args.step == "step1":
        Unb = unb.Unbinding()
        if os.path.isfile(Unb.checkpoint):
            Unb = Unb.load()
            Unb.newCycle()
        else:
            Unb.traj_length = args.ns
            Unb.ligresname = args.lig
            Unb.traj_length = args.ns
            Unb.maxdist = args.maxdist / 10
            Unb.set_top(args.top)
            Unb.set_namd_template()
        if args.trajectory is not None:
            Unb.cycle = int(args.trajectory) + 1
        c = rp.Cycle(Unb)
        c.readDCD(Unb.top)
        if args.writeDCD:
            c.saveNewDCD(int(args.stride))
        c.getNeighbour(args.lig, args.cutoff, Unb.clusters)
        c.getClusters(Unb.clusters)
    elif args.step == "step2":
        data = cg.readData("mean_analysis.dat")
        g = cg.createPairs(data)
        cg.printPairs(g)
        cg.printSelection(data)
    elif args.step == "step4":
        c = colvar.Contact()
        c.readClusters("new-distances.dat", COM=True)
        c.writeNAMDcolvar("sum.col", STEP4=True, traj_length=int(args.ns))
    elif args.step == "step5":
        c = colvar.Contact()
        c.readClusters("new-distances.dat", COM=True)
        c.printDistances("string")
    elif args.step == "string":
        Unb = unb.Unbinding()
        if not os.path.isfile(Unb.checkpoint):
            print("Checkpoint file is not found. Exiting")
            sys.exit(0)
        Unb = Unb.load()
        c = rp.Cycle(Unb)
        c.getAllPairs(Unb)
        try:
            c.readDCD(Unb.top)
        except rp.DCDnotReadable:
            Unb.writeOutput("DCD file cannot be read for cycle {:d}".format(Unb.cycle))
            sys.exit(0)
        c.createContact()
        c.contact.prepareString(Unb)
    else:
        Unb = unb.Unbinding()
        if args.cumulative:
            if args.trajectory is None:
                print("With --cumulative option, please specify the last trajectory to be processed with option -t.")
                sys.exit(0)
            Unb.reprocess(args)
        else:
            if os.path.isfile(Unb.checkpoint):
                Unb = Unb.load()
                Unb.newCycle()
            else:
                Unb.readClusters()
                Unb.ligresname = args.lig
                Unb.traj_length = args.ns
                Unb.maxdist = args.maxdist / 10
                Unb.set_top(args.top)
                Unb.set_namd_template()
            if args.trajectory is not None:
                Unb.cycle = int(args.trajectory) + 1
            c = rp.Cycle(Unb)
            try:
                c.readDCD(Unb.top)
            except rp.DCDnotReadable:
                Unb.writeOutput("DCD file cannot be read for cycle {:d}".format(Unb.cycle - 1))
                sys.exit(0)
            if args.writeDCD:
                c.saveNewDCD(int(args.stride))
            c.getNeighbour(args.lig, args.cutoff, Unb.clusters)
            c.getClusters(Unb.clusters)
            Unb.history(c)
            c.createContact()
            out.trackDistances(Unb)
            if not args.processonly:
                out.vmdRep(Unb)
                c.setupCycle()
                c.contact.writeNAMDcolvar(
                    os.path.join(c.wrkdir, "traj_{:d}".format(c.number), "sum_{:d}.col".format(c.number)),
                    traj_length=int(Unb.traj_length))
            if Unb.cycle == 1:
                Unb.writeOutput(out.header())
            Unb.writeOutput(out.cycle(c))
            if not args.nosave:
                Unb.save()
