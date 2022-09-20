# from builtins import Exception
import os
import sys
import mdtraj as md
import numpy as np
import scr3_python_colvar_generator as cg
import CoM_colvars as cv
sys.path.append(os.path.dirname(os.path.abspath(__file__)))


class PDBnotReadable(Exception):
    pass


class DCDnotReadable(Exception):
    pass


class PDB:
    def __init__(self):
        self.structure = None
        return

    def read(self, file):
        try:
            self.structure = md.load(file)
        except:
            raise PDBnotReadable

    def readDCD(self, dcd, pdb, indices=None):
        try:
            self.structure = md.load_dcd(dcd, top=pdb, atom_indices=indices)
            self.structure = self.structure.image_molecules(inplace=True)
        except:
            raise PDBnotReadable


class Cycle:
    def __init__(self, unbinding):
        self.number = unbinding.cycle
        self.wrkdir = unbinding.wrkdir
        self.traj_length = unbinding.traj_length
        self.prevtraj = None
        self.pairs = []
        self.contact = None

    def readDCD(self, top):
        try:
            self.prevtraj = md.load_dcd(os.path.join("traj_{0:d}/traj_{0:d}.dcd".format(self.number - 1)), top=top)
            self.prevtraj = self.prevtraj.image_molecules(inplace=True)
        except:
            raise DCDnotReadable
        self.prevtraj = self.prevtraj.image_molecules(inplace=True)

    def saveNewDCD(self, stride):
        self.prevtraj[range(0, self.prevtraj.n_frames, stride)].save_dcd(os.path.join("traj_{0:d}/traj_{0:d}-wrapped.dcd".format(self.number - 1)))
        return

    def getNeighbour(self, ligres, cutoff, ligandClusters):
        lig = self.prevtraj.top.select("resname {:s} and not type H".format(ligres))
        protein = self.prevtraj.top.select("protein and not type H")
        neighbours = []
        for l in lig:
            neighbours.append(md.compute_neighbors(self.prevtraj, cutoff=0.1*cutoff, haystack_indices=protein, query_indices=[l]))
        data = []
        for i in range(self.prevtraj.n_frames):
            frame = []
            for j in range(len(lig)):
                for pa in neighbours[j][i]:
                    frame.append([lig[j],
                                 self.prevtraj.top._atoms[lig[j]].residue.resSeq,
                                 self.prevtraj.top._atoms[lig[j]].name,
                                 pa,
                                 self.prevtraj.top._atoms[pa].residue.resSeq,
                                 self.prevtraj.top._atoms[pa].name,
                                 10 * np.linalg.norm(self.prevtraj._xyz[i, lig[j]] -
                                                     self.prevtraj._xyz[i, pa]),
                                 self.prevtraj.top._atoms[pa].residue.name,
                                 ligres])
            data.append(frame)
        self.pairs = cg.createPairs(data)
        for i in range(len(self.pairs)):
            self.pairs[i].getProteinClusterAtoms(self.prevtraj.top)
            self.pairs[i].getLigandClusterAtoms(self.prevtraj.top, ligandClusters)
        self.removeDuplicates()
        return

    def getClusters(self, ligandClusters):
        """Obsolete as hes to be done earlier"""
        for i in range(len(self.pairs)):
            self.pairs[i].getProteinClusterAtoms(self.prevtraj.top)
            self.pairs[i].getLigandClusterAtoms(self.prevtraj.top, ligandClusters)
        self.removeDuplicates()
        return

    def removeDuplicates(self):
        remove = []
        for i in range(len(self.pairs)):
            for j in range(i):
                if self.pairs[j].hasAtom(self.pairs[i].atom["index"]) and self.pairs[j].hasAtom(self.pairs[i].ligand_atom["index"]):
                    union = sorted(list(set(self.pairs[i].mask) | set(self.pairs[j].mask)))
                    value = []
                    for k in union:
                        try:
                            iv = self.pairs[i].value[self.pairs[i].mask.index(k)]
                        except ValueError:
                            value.append(self.pairs[j].value[self.pairs[j].mask.index(k)])
                            continue
                        try:
                            jv = self.pairs[j].value[self.pairs[j].mask.index(k)]
                        except ValueError:
                            value.append(self.pairs[i].value[self.pairs[i].mask.index(k)])
                            continue
                        value.append(min([iv, jv]))
                    self.pairs[j].value = value
                    self.pairs[j].mask = union
                    self.pairs[j].count = len(self.pairs[j].mask)
                    remove.append(i)
        uniquepairs = []
        for i in range(len(self.pairs)):
            if i not in remove and self.pairs[i].count > (0.5 * self.prevtraj.n_frames):
                uniquepairs.append(self.pairs[i])
        self.pairs = uniquepairs
        return

    def getAllPairs(self, Unb):
        IDfound = []
        for c in Unb.pairs:
            for pair in c:
                if pair.ID not in IDfound:
                    self.pairs.append(pair)
                    IDfound.append(pair.ID)
        return

    def createContact(self, COM = True):
        self.contact = cv.Contact()
        self.contact.pdb.structure = self.prevtraj[-1]
        if COM:
            for pair in self.pairs:
                if len(pair.ligandClusterAtoms) == 0:
                    l = self.contact.getLigandCluster(pair.ligand_atom["index"])
                else:
                    l = self.contact.getLigandCluster(pair.ligandClusterAtoms[0]["index"])
                    for a in pair.ligandClusterAtoms[1:]:
                        if not self.contact.ligandClusters[l].hasAtom(a["index"]):
                            self.contact.ligandClusters[l].addAtom(a["index"])
                if len(pair.proteinClusterAtoms) == 0:
                    p = self.contact.getProteinCluster(pair.atom["index"])
                else:
                    p = self.contact.getProteinCluster(pair.proteinClusterAtoms[0]["index"])
                    for a in pair.proteinClusterAtoms[1:]:
                        if not self.contact.proteinClusters[p].hasAtom(a["index"]):
                            self.contact.proteinClusters[p].addAtom(a["index"])
                self.contact.ligandClusters[l].contacts += 1
                self.contact.proteinClusters[p].contacts += 1
                self.contact.associations.append((l, p))
                self.contact.weight.append(1.0)
        else:
            for pair in self.pairs:
                if len(pair.ligandClusterAtoms) == 0:
                    l = self.contact.getLigandCluster(pair.ligand_atom["index"])
                else:
                    for a in pair.ligandClusterAtoms:
                        l = self.contact.getLigandCluster(a["index"])
                if len(pair.proteinClusterAtoms) == 0:
                    p = self.contact.getProteinCluster(pair.atom["index"])
                else:
                    for a in pair.proteinClusterAtoms[1:]:
                        p = self.contact.getProteinCluster(a["index"])
            for lc in self.contact.ligandClusters:
                for pair in self.pairs:
                    for l in pair.ligandClusterAtoms:
                        if l["index"] == lc.atoms[0]:
                            if len(pair.proteinClusterAtoms) == 0:
                                lc.contacts += 1
                                self.contact.associations.append((self.contact.ligandClusters.index(lc), self.contact.getProteinCluster(pair.atom["index"])))
                            else:
                                lc.contacts += len(pair.proteinClusterAtoms)
                                for a in pair.proteinClusterAtoms:
                                    self.contact.associations.append((self.contact.ligandClusters.index(lc),
                                                                      self.contact.getProteinCluster(a["index"])))
            for pc in self.contact.proteinClusters: # count for protein clusters doesn't work
                for pair in self.pairs:
                    for p in pair.proteinClusterAtoms:
                        if p["index"] == pc.atoms[0]:
                            if len(pair.ligandClusterAtoms) == 0:
                                pc.contacts += 1
                            else:
                                pc.contacts += len(pair.ligandClusterAtoms)
        return

    def setupCycle(self):
        try:
            os.mkdir(os.path.join(self.wrkdir, "traj_{:d}".format(self.number)))
        except OSError:
            print("WARNING: traj_{:d} already exists. Files might be overwritten in the folder.".format(self.number))
            print("WARNING: Rename the folder to back it up, or remove.")
            sys.exit(0)
        self.writeNamdInput()

    def writeNamdInput(self):
        input = """############
# Force-Field Parameters
amber              yes
parmfile           ../toppar/complex.prmtop
coordinates        ../toppar/complex.pdb

set temp           298;
outputName         traj_{0:d} # base name for output from this run

set inputname      traj_{1:d}.restart;
binCoordinates     ../traj_{1:d}/$inputname.coor;
binVelocities      ../traj_{1:d}/$inputname.vel;
extendedSystem     ../traj_{1:d}/$inputname.xsc;
restartfreq        500;
dcdfreq            1000;
dcdUnitCell        yes;
xstFreq            500;
outputEnergies     125;
outputTiming       500;

# These are specified by AMBER
exclude             scaled1-4
1-4scaling          0.833333
switching           on
vdwForceSwitching   yes;

# You have some freedom choosing the cutoff
cutoff              12.0;
switchdist          10.0;
pairlistdist        16.0;
stepspercycle       20;
pairlistsPerCycle    2;

# Integrator Parameters
timestep            2.0;
rigidBonds          all;
nonbondedFreq       1;
fullElectFrequency  1;

wrapWater           on;
wrapAll             on;
wrapNearest        off;

# PME (for full-system periodic electrostatics)
source ../toppar/checkfft.str

PME                yes;
PMEInterpOrder       6;
PMEGridSizeX     $fftx;
PMEGridSizeY     $ffty;
PMEGridSizeZ     $fftz;

# Constant Pressure Control (variable volume)
useGroupPressure       yes;
useFlexibleCell         no;
useConstantRatio        no;
langevinPiston          on;
langevinPistonTarget  1.01325;
langevinPistonPeriod  50.0;
langevinPistonDecay   25.0;
langevinPistonTemp   $temp;

# Constant Temperature Control
langevin                on;
langevinDamping        1.0;
langevinTemp         $temp;
langevinHydrogen       off;

colvars on
colvarsConfig sum_{0:d}.col

# run
run                {2:d};

""".format(self.number, self.number - 1, int(self.traj_length * 5E5))
        with open(os.path.join(self.wrkdir, "traj_{:d}".format(self.number), "traj_{:d}.inp".format(self.number)), "w") as f:
            f.write(input)
        return

def main():
    example = PDB()
    example.read("examples/4fku.pdb")
    print("stop")


if __name__ == "__main__": main()