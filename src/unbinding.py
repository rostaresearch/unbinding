import os
import sys
import pickle
from read_ligand import readLigandClusters
import read_pdb as rp
import output as out
import numpy as np
import mdtraj as md
# import matplotlib.pyplot as plt


class Unbinding:
    def __init__(self):
        self.wrkdir = os.getcwd()
        self.checkpoint = os.path.join(self.wrkdir, ".checkpoint")
        self.output = os.path.join(self.wrkdir, "unbinding.out")
        # print header
        prmtop = os.path.join(self.wrkdir, "toppar", "complex.prmtop")
        psf = os.path.join(self.wrkdir, "toppar", "complex.psf")
        pdb = os.path.join(self.wrkdir, "toppar", "complex.pdb")
        if os.path.isfile(prmtop):
            self.top = prmtop
        elif os.path.isfile(psf):
            self.top = psf
        elif not os.path.isfile(pdb):
            self.top = pdb
        else:
            self.writeOutput("No topology file is found in {:s}".format(os.path.join(self.wrkdir, "toppar")))
            self.writeOutput("Topology can be: complex.prmtop, complex.psf, complex.pdb")
            sys.exit(0)
        self.cycle = 1
        self.clusters = None
        self.traj_lenght = 10
        self.pairs = []
        self.maxdist = None     # discarding distance larger then this (in nm)
        self.meanfluct = 1      # discarding distance with larger mean fluctuation (in nm)
        self.N_pairs = 0
        self.ligresname = ""

    def history(self, cycle):
        if self.cycle == 1:
            self.pairs.append(cycle.pairs)
            for i in range(len(self.pairs[0])):
                self.pairs[-1][i].ID = i + 1
                self.N_pairs += 1
        else:
            discard = []
            for oldpair in self.pairs[self.cycle - 2]:
                lig = [oldpair.ligand_atom["index"]]
                for a in oldpair.ligandClusterAtoms:
                    lig.append(a["index"])
                lig = list(set(lig))
                pro = [oldpair.atom["index"]]
                for a in oldpair.proteinClusterAtoms:
                    pro.append(a["index"])
                pro = list(set(pro))
                c1 = md.compute_center_of_mass(cycle.prevtraj.atom_slice(lig))
                c2 = md.compute_center_of_mass(cycle.prevtraj.atom_slice(pro))
                dist = np.linalg.norm(c1 - c2, axis=1)
                smdist = rolling_mean(dist, 50)
                # fluc = np.abs(dist - smdist)
                # print(np.mean(fluc), np.std(fluc))
                # p0 = plt.plot(range(5000), dist, range(5000), smdist)
                # p = plt.plot(fluc)
                # plt.show()
                if np.mean(np.abs(dist - smdist)) > self.meanfluct:
                    discard.append(self.pairs[self.cycle - 2].index(oldpair))
                    self.writeOutput("Distance with id {0:02d} is excluded since the fluctuation is {1:.2f}".format(oldpair.ID, np.mean(np.abs(dist - smdist))))
                elif np.mean(dist[-int(len(dist)/4):]) > self.maxdist:
                    discard.append(self.pairs[self.cycle - 2].index(oldpair))
                    self.writeOutput(
                        """Distance with id {0:02d} is excluded since the the avarage distance
                        in the last quarter ot the simulation is {1:.2f} A""".format(
                            oldpair.ID, 10 * np.mean(dist[-int(len(dist) / 4):])))
            for oldpair in self.pairs[self.cycle - 2]:
                if self.pairs[self.cycle - 2].index(oldpair) in discard:
                    continue
                for newpair in cycle.pairs:
                    if newpair.hasAtom(oldpair.atom["index"]) and newpair.hasAtom(oldpair.ligand_atom["index"]):
                        discard.append(self.pairs[self.cycle - 2].index(oldpair))
            for oldpair in self.pairs[self.cycle - 2]:
                if self.pairs[self.cycle - 2].index(oldpair) in discard:
                    continue
                cycle.pairs.append(oldpair)
            self.pairs.append(cycle.pairs)
            for newpair in self.pairs[self.cycle - 1]:
                ASSIGNED = False
                for c in range(len(self.pairs) - 1):
                    for oldpair in self.pairs[c]:
                        if not ASSIGNED:
                            if newpair.hasAtom(oldpair.atom["index"]) and newpair.hasAtom(oldpair.ligand_atom["index"]):
                                newpair.ID = oldpair.ID
                                ASSIGNED = True
                if not ASSIGNED:
                    newpair.ID = self.N_pairs + 1
                    self.N_pairs += 1
        return

    def readClusters(self):
        try:
            self.clusters = readLigandClusters(os.path.join(self.wrkdir, "toppar", "LIG_clusters.dat"))
        except IOError:
            self.writeOutput("No cluster file is found in {:s}".format(os.path.join(self.wrkdir, "toppar", "LIG_clusters.dat")))
            open(os.path.join(self.wrkdir, "toppar", "LIG_clusters.dat"), 'a').close()
            self.writeOutput("Empty ligand cluster file is created. Check to manual to ddefine clusters.")
            self.clusters = readLigandClusters(os.path.join(self.wrkdir, "toppar", "LIG_clusters.dat"))
        return

    def newCycle(self):
        self.cycle += 1
        return

    def refreshWorkdir(self, wd):
        old = self.wrkdir
        self.wrkdir = wd
        self.top = self.top.replace(old, wd)
        self.output = self.output.replace(old, wd)
        return

    def writeOutput(self, string):
        with open(self.output, 'a') as f:
            f.write("{:s}\n".format(string))
        return

    def save(self):
        cf = open(self.checkpoint, 'wb')
        pickle.dump(self, cf)
        cf.close()
        return

    def load(self):
        cf = open(self.checkpoint, 'rb')
        loaded = pickle.load(cf)
        cf.close()
        return loaded

    def reprocess(self, args):
        self.readClusters()
        self.ligresname = args.lig
        self.traj_length = int(args.ns)
        self.maxdist = args.maxdist / 10
        self.writeOutput(out.header())

        for i in range(int(args.trajectory) + 1):
            if i != 0:
                self.newCycle()
            c = rp.Cycle(self)
            try:
                c.readDCD(self.top)
            except rp.DCDnotReadable:
                self.writeOutput("DCD file cannot be read for cycle {:d}".format(i))
                sys.exit(0)
            if args.writeDCD:
                c.saveNewDCD(int(args.stride))
            c.getNeighbour(args.lig, args.cutoff, self.clusters)
            c.getClusters(self.clusters)
            self.history(c)
            c.createContact()
            out.trackDistances(self)
            if not args.processonly and i == int(args.trajectory):
                out.vmdRep(self)
                c.setupCycle()
                c.contact.writeNAMDcolvar(
                    os.path.join(c.wrkdir, "traj_{:d}".format(c.number), "sum_{:d}.col".format(c.number)),
                    traj_length=int(self.traj_length))
            self.writeOutput(out.cycle(c))
        if not args.nosave:
            self.save()
        return


def rolling_mean(x, N):
    newx = np.empty(len(x))
    for i in range(len(x)):
        if i-N < 0:
            newx[i] = np.mean(x[0:i + N + 1])
        else:
            newx[i] = np.mean(x[i - N:i + N + 1])
    return newx
