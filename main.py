import argparse
import os
import sys

import src.cycle
import src.unbinding as unb
import src.read_pdb as rp
import src.output as out
from src.contact import UnboundException


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', "--trajectory")
    parser.add_argument('-l', "--lig", default="LIG")
    parser.add_argument("--top", default="find")
    parser.add_argument('-c', "--cutoff", default=3.5, type=float, help="Initial cutoff for identifying neighbours in A")
    parser.add_argument('-m', "--maxdist", default=9, type=float, help="distance for exclusion fo pairs in A")
    parser.add_argument('-ns', default=10, type=int)
    parser.add_argument("--cumulative", action='store_true', default=False,
                        help='Reprocess all the trajectories up to the one specified by "-t". ')
    parser.add_argument('--writeDCD', action='store_true', default=False, help='write the strided DCD')
    parser.add_argument('-s', "--stride", default=5, type=int)
    parser.add_argument('-p', "--processonly", action='store_true', default=False,
                        help='Do not write VMD and NAMD input. Other outputs will be written.')
    parser.add_argument("--nosave", action='store_true', default=False,
                        help='Do not save the checkpoint. For debug only.')
    parser.add_argument("--report", action='store_true', default=False,
                        help='Report status and exit.')
    parser.add_argument("--auto", action='store_true', default=False,
                        help='Go to the background and restart after namd finished.')
    parser.add_argument("--namd", help="NAMD submission script, taking <name>.inp as an input and writing to"
                                       " <name>.out. Use with --auto.")
    parser.add_argument('--maxiter', default=25, type=int, help="Maximum number of iterations. Use with --auto.")
    parser.add_argument("--string", action='store_true', default=False)
    args = parser.parse_args()
    Unb = unb.Unbinding()
    if args.report:
        if os.path.isfile(Unb.checkpoint):
            Unb = Unb.load()
            Unb.report()
        else:
            print("The is no checkpoint file to report of.")
    elif args.string:
        if os.path.isfile(Unb.checkpoint):
            Unb = Unb.load()
            c = src.cycle.Cycle(Unb)
            c.getAllPairs(Unb)
            c.createContact()
            c.contact.prepareString(Unb)
            pass
        else:
            print("The is no checkpoint file, please reprocess first to create a string.")
    elif args.cumulative:
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
            Unb.cutoff = args.cutoff
        if args.auto:
            if args.namd is None:
                print("With --auto option, please specify a script which submits namd as <script> namd.inp,"
                      " producing namd.out as an output.")
                sys.exit(0)
            if os.fork():
                sys.exit(0)
            Unb.auto(args)
        if args.trajectory is not None:
            Unb.cycle = int(args.trajectory) + 1
        c = src.cycle.Cycle(Unb)
        try:
            c.readDCD(Unb.top)
        except rp.DCDnotReadable:
            Unb.writeOutput("DCD file cannot be read for cycle {:d}".format(Unb.cycle - 1))
            sys.exit(0)
        if args.writeDCD:
            c.saveNewDCD(int(args.stride))
        c.getNeighbour(Unb.ligresname, Unb.cutoff, Unb.clusters)
        c.getClusters(Unb.clusters)
        Unb.history(c)
        c.createContact()
        out.trackDistances(Unb)
        if not args.processonly:
            out.vmdRep(Unb)
            try:
                c.setupCycle()
                c.contact.writeNAMDcolvar(
                    os.path.join(c.wrkdir, "traj_{:d}".format(c.number), "sum_{:d}.col".format(c.number)),
                    traj_length=int(Unb.traj_length))
            except UnboundException:
                Unb.writeOutput("No contact remained in the colvar, ligand unbound!")
                sys.exit(0)
        if Unb.cycle == 1:
            Unb.writeOutput(out.header())
        Unb.writeOutput(out.cycle(c))
        if not args.nosave:
            Unb.save()
