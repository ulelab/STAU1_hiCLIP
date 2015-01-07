'''
Created on May 20, 2013

@author: yoichiro
'''
#!/bin/env python

# i1_map_HTSeq.py
#
# Created by Yoichiro Sugimoto on 14/08/2012.
# Copyright 2012 __MyCompanyName__. All rights reserved.
# Command example: qsub -j y -b y -pe smp 8 ~/Python-2.7.1/python ~/i1_map_HTSeq.py -d /path/to/folder -p /path/to/python
# qsub -j y -b y -pe smp 8 ~/Python-2.7.1/python /netscr/yoichiro/hiCLIP_project/first_submission/script/individual_anlaysis/i1_map_HTSeq.py -d /netscr/yoichiro/hiCLIP_project -p ~/Python-2.7.1 -r ~/R-2.15.1

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

def listdir_fullpath_extension(directory, extension):
    FileLists = [os.path.join(directory, f) for f in os.listdir(directory)]
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

def run_i1(hiCLIPdir, pythondir, Rdir):

    hiCLIP_mapping_command = pythondir + "/python " + hiCLIPdir + "/inst/Python/runHybridAnalysis_transcripts-2.0.0.py -i " + hiCLIPdir + "/data/unprocessed/fastq/hiCLIP -o " + hiCLIPdir + "/results/mapping/hiCLIP -b " + hiCLIPdir + "/data/unprocessed/barcode/hiCLIP -d " + hiCLIPdir + " -p " + pythondir + " -r " + Rdir
    runUNIXCommands(hiCLIP_mapping_command)

    RP_mapping_command = pythondir + "/python " + hiCLIPdir + "/inst/Python/runRibosomeProfilingAnalysis-4.0.py -i " + hiCLIPdir + "/data/unprocessed/fastq/ribosome_profiling -o " + hiCLIPdir + "/results/mapping/ribosome_profiling -b " + hiCLIPdir + "/data/unprocessed/barcode/ribosome_profiling -d " + hiCLIPdir + " -p " + pythondir + " -r " + Rdir
    runUNIXCommands(RP_mapping_command)

    mRNASeq_mapping_command = pythondir + "/python " + hiCLIPdir + "/inst/Python/runRibosomeProfilingAnalysis-4.0.py -i " + hiCLIPdir + "/data/unprocessed/fastq/mRNASeq -o " + hiCLIPdir + "/results/mapping/mRNASeq -b " + hiCLIPdir + "/data/unprocessed/barcode/mRNASeq -d " + hiCLIPdir + " -p " + pythondir + " -r " + Rdir
    runUNIXCommands(mRNASeq_mapping_command)

    hnRNPc_mapping_command = pythondir + "/python " + hiCLIPdir + "/inst/Python/runRibosomeProfilingAnalysis-4.0.py -i " + hiCLIPdir + "/data/unprocessed/fastq/hnRNPc -o " + hiCLIPdir + "/results/mapping/hnRNPc -b " + hiCLIPdir + "/data/unprocessed/barcode/hnRNPc -d " + hiCLIPdir + " -p " + pythondir + " -r " + Rdir
    runUNIXCommands(hnRNPc_mapping_command)
    
    STAU1_iCLIP_mapping_command = pythondir + "/python " + hiCLIPdir + "/inst/Python/runRibosomeProfilingAnalysis-4.0.py -i " + hiCLIPdir + "/data/unprocessed/fastq/STAU1_iCLIP -o " + hiCLIPdir + "/results/mapping/STAU1_iCLIP -b " + hiCLIPdir + "/data/unprocessed/barcode/STAU1_iCLIP -d " + hiCLIPdir + " -p " + pythondir + " -r " + Rdir
    runUNIXCommands(STAU1_iCLIP_mapping_command)

if __name__ == "__main__":
    from optparse import OptionParser

    parser = OptionParser()

    parser.add_option("-d", "--hiCLIP_project_dir", dest="hiCLIPDir",
                  help="Path to hiCLIP_project")
    parser.add_option("-p", "--python_dir", dest="pythonDir",
                  help="Path to Python")
    parser.add_option("-r", "--R_dir", dest="RDir",
                  help="Path to R")

    (options, args) = parser.parse_args()

    hiCLIPdir = options.hiCLIPDir
    pythondir = options.pythonDir
    Rdir = options.RDir

    run_i1(hiCLIPdir, pythondir, Rdir)
