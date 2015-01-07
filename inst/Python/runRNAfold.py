#!/bin/env python

# runRNAfold.py
# Run RNAfold for the structural motif analysis
#
#  qsub -b y -pe smp 8 -j y /lmb/home/yoichiro/Python-2.7.1/python ~/test/runRNAfold.py -i ~/test -o /netscr/yoichiro/test6 -d /netscr/yoichiro/STAU1_hiCLIP
# Created by Yoichiro Sugimoto on 02/01/2014.
# Copyright 2014 __MyCompanyName__. All rights reserved.

import os
import sys, subprocess

class Unbuffered:
   def __init__(self, stream):
       self.stream = stream
   def write(self, data):
       self.stream.write(data)
       self.stream.flush()
   def __getattr__(self, attr):
       return getattr(self.stream, attr)

sys.stdout=Unbuffered(sys.stdout)

def listdir_fullpath(d):
   return [os.path.join(d, f) for f in os.listdir(d)]

def listdir_fullpath_extension(directory, extension):
   FileLists = listdir_fullpath(directory)
   FileLists = [each for each in FileLists if each.endswith(extension)]
   return FileLists

def listdir_fullpath_directories(directory):
    FileLists = [os.path.join(directory, f) for f in os.listdir(directory)]
    FileLists = [each for each in FileLists if (each.find('.') == -1)]
    return FileLists

def runUNIXCommands(commands):
    print commands
    try:
        retcode = subprocess.call(commands, shell=True)
        if retcode < 0:
            print >>sys.stderr, "Child was terminated by signal", -retcode
        else:
            print >>sys.stderr, "Child returned", retcode
    except OSError, e:
        print >>sys.stderr, "Execution failed:", e
        
def runRNAFold(inputfile, outputdir, hiCLIPdir):
    os.chdir(outputdir)
    RNAfoldCommand = hiCLIPdir + "/inst/bin/ViennaRNA-2.1.3/Progs/RNAfold -p -d2 --noPS < " + inputfile + " > " + outputdir + "/result.db.txt"
    runUNIXCommands(RNAfoldCommand)

if __name__ == "__main__":
    from optparse import OptionParser

    parser = OptionParser()

    parser.add_option("-i", "--input", dest="inputFile",
                  help="Input directory which contains fastq files, the fastq files shoule be anarchived and the extension should be .fq", metavar="INPUT")
    parser.add_option("-o", "--output",
                  dest="outputDirectry",
                  help="Directry for output, have to spply")
    parser.add_option("-d", "--hiCLIP_project_dir", dest="hiCLIPDir",
                  help="Path to hiCLIP_project")

    (options, args) = parser.parse_args()

    inputfile = options.inputFile
    outputdir = options.outputDirectry
    hiCLIPdir = options.hiCLIPDir
    
    runUNIXCommands("mkdir -p " + outputdir)

    runRNAFold(inputfile, outputdir, hiCLIPdir)
