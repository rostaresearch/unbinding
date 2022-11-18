import glob
import logging
import os
import sys
import argparse
from subprocess import Popen
from time import sleep


def tail(f, BLOCK_SIZE = 256):
    f.seek(0, 2)
    size = f.tell()
    if size < BLOCK_SIZE:
        BLOCK_SIZE = size
    f.seek(-BLOCK_SIZE, 2)
    a = f.readlines()
    rv = ""
    for b in a:
        rv += str(b)
    return rv


class NoOutputException(Exception):
    pass


class StillRunningException(Exception):
    pass


class g09calcHandler():
    def __init__(self, files, args):
        self.inputs = files
        self.all = len(files)
        self.submittedInputs = []
        self.foundOutputs = []
        self.termOutputs = []
        self.wrongOutputs = []
        self.toBeRemoved = []
        self.workdir = os.getcwd()
        self.command = args.command
        self.options = args.options
        self.NCalc = args.ncalc

    def submit(self,inp):
        os.chdir(self.workdir)
        Popen([str(self.command), inp, self.options])
        #logging.info("PID" + str(os.getpid()) + ": " + inp + " was submitted")
        self.submittedInputs.append(inp)
        self.inputs.remove(inp)

    def seekOut(self,sinp):
        filename = sinp.split('.')[0]
        if len(glob.glob(os.path.join(self.workdir,filename) + ".log")) == 1:
            self.foundOutputs.append(glob.glob(os.path.join(self.workdir, filename) + "*log")[0])
            self.toBeRemoved.append(sinp)
            return
        elif len(glob.glob(os.path.join(self.workdir, filename) + ".out")) == 1:
            self.foundOutputs.append(glob.glob(os.path.join(self.workdir, filename) + "*out")[0])
            self.toBeRemoved.append(sinp)
            return
        else:
            raise NoOutputException

    def checkOut(self,fout):
        with open(fout, 'rb') as outfile:
            if "termination" in tail(outfile):
                if "Normal" in tail(outfile):
                    self.termOutputs.append(fout)
                    self.toBeRemoved.append(fout)
                    return
                if "Error" in tail(outfile):
                    self.wrongOutputs.append(fout)
                    return
                else:
                    logging.info("PID" + str(os.getpid()) + ": "+fout+" stopped unexpectedly.")
                    self.wrongOutputs.append(fout)
                    return
            else:
                raise StillRunningException

    def watch(self):
        for inp in self.inputs:
            try:
                self.seekOut(inp)
            except NoOutputException:
                pass
        for rem in self.toBeRemoved:
            self.inputs.remove(rem)
        self.toBeRemoved = []
        for fout in self.foundOutputs:
            try:
                self.checkOut(fout)
            except StillRunningException:
                pass
        for rem in self.toBeRemoved:
            self.foundOutputs.remove(rem)
        self.toBeRemoved = []
        logging.info("PID" + str(os.getpid()) + ": Watch started. The current state of calculations:")
        logging.info("PID" + str(os.getpid()) + ": total calculations: " + str(self.all))
        logging.info("PID" + str(os.getpid()) + ": submitted: " + str(len(self.submittedInputs)))
        logging.info("PID" + str(os.getpid()) + ": running: " + str(len(self.foundOutputs)))
        logging.info("PID" + str(os.getpid()) + ": normally terminated: " + str(len(self.termOutputs)))
        logging.info("PID" + str(os.getpid()) + ": terminated with error: " + str(len(self.wrongOutputs)))
        logging.info("PID" + str(os.getpid()) + ": remaining: " + str(len(self.inputs)))
        done = len(self.wrongOutputs) + len(self.termOutputs)
        cycles = 0
        while done != self.all:
            running = len(self.submittedInputs) + len(self.foundOutputs)
            while running < self.NCalc:
                if len(self.inputs) == 0:
                    break
                else:
                    self.submit(self.inputs[0])
                    running = len(self.submittedInputs) + len(self.foundOutputs)
            sleep(15)
            for sinp in self.submittedInputs:
                try:
                    self.seekOut(sinp)
                except NoOutputException:
                    pass
            for rem in self.toBeRemoved:
                self.submittedInputs.remove(rem)
            self.toBeRemoved = []
            sleep(15)
            for fout in self.foundOutputs:
                try:
                    self.checkOut(fout)
                except StillRunningException:
                    pass
            for rem in self.toBeRemoved:
                self.foundOutputs.remove(rem)
            self.toBeRemoved = []
            cycles += 1
            if cycles == 40:
                logging.info("PID" + str(os.getpid()) + ": total calculations: " + str(self.all))
                logging.info("PID" + str(os.getpid()) + ": submitted: " + str(len(self.submittedInputs)))
                logging.info("PID" + str(os.getpid()) + ": running: " + str(len(self.foundOutputs)))
                logging.info("PID" + str(os.getpid()) + ": normally terminated: " + str(len(self.termOutputs)))
                logging.info("PID" + str(os.getpid()) + ": terminated with error: " + str(len(self.wrongOutputs)))
                logging.info("PID" + str(os.getpid()) + ": remaining: " + str(len(self.inputs)))
                cycles = 0
            done = len(self.wrongOutputs) + len(self.termOutputs)
        logging.info("\nPID" + str(os.getpid()) + ": The shift is over, watchman may rest.")
        logging.info("PID" + str(os.getpid()) + ": total calculations: " + str(self.all))
        logging.info("PID" + str(os.getpid()) + ": normally terminated: " + str(len(self.termOutputs)))
        logging.info("PID" + str(os.getpid()) + ": terminated with error: " + str(len(self.wrongOutputs))+"\n")
        if len(self.wrongOutputs) != 0:
            logging.info("Wrong outputs are:")
            for wout in self.wrongOutputs:
                logging.info(wout)
        sys.exit(0)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', "--ncalc")
    parser.add_argument('-c', "--command")
    parser.add_argument('-o', "--options")
    parser.add_argument('-e', "--extension")
    parser.add_argument('-l', "--list")
    args = parser.parse_args()
    # set up defaults
    if args.ncalc == None: args.ncalc = 99
    else: args.ncalc = int(args.ncalc)
    if args.command == None: args.command = "rung09_quiet"
    if args.options == None: args.options = ""
    # logfile
    logging.basicConfig(filename="watchlog", level=logging.INFO, format='%(asctime)s %(message)s')
    # set up files
    if args.list == None and args.extension == None:
        print("Try to submit all files. Continue? Y/n")
        while True:
            answer = input("\n")
            if answer == 'Y':
                break
            elif answer == 'n':
                sys.exit(0)
        files = [f for f in os.listdir('.') if os.path.isfile(f)]
        files.remove("watchlog")
    elif args.list != None:
        files = []
        try:
            with open(args.list,'r') as f:
                for line in f:
                    for s in line.split():
                        if os.path.isfile(s):
                            files.append(s)
        except:
            print("Cannot parse list.", args.list)
            sys.exit(0)
    elif args.extension != None:
        files = [f for f in os.listdir('.') if os.path.isfile(f) and f.split('.')[-1] == args.extension]
    if len(files) == 0:
        print("There is no file to submit.")
        sys.exit(0)
    handler = g09calcHandler(files, args)
    if os.fork():
        sys.exit(0)
    handler.watch()


if __name__ == "__main__":
    main()
