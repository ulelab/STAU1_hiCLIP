#!/bin/env python

# runHybridAnalysis_transcripts-2.0.0.py
# Strategt was summarised in /Users/yoichiro/Documents/STAU1_project/transcripts_mapping/Method.docx
# The fq file created by this file should be categorised as RUT (Ribosome RNA unmapped and trimmed to 26 nucleotides)
#
#  qsub -b y -pe smp 8 -j y /lmb/home/yoichiro/Python-2.7.1/python /netscr/yoichiro/STAU1_hiCLIP/inst/Python/RNAfold.py -i ~/test -o ~/test -d /netscr/yoichiro/STAU1_hiCLIP
# Created by Yoichiro Sugimoto on 14/08/2012.
# Copyright 2012 __MyCompanyName__. All rights reserved.

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

def runBowtieForRibosomeProfiling(inputdir, outputdir, hiCLIPdir):
    os.chdir(outputdir)
    inputfile = inputdir + "/DB_constraint.faconst"
    commandRNAcofoldHybrid = hiCLIPdir + "/inst/bin/ViennaRNA-2.1.3/Progs/RNAfold  -p -d2 -C < " + inputfile + " > /dev/null"
    runUNIXCommands(commandRNAcofoldHybrid)
    
if __name__ == "__main__":
    from optparse import OptionParser

    parser = OptionParser()

    parser.add_option("-i", "--input", dest="inputDir",
                  help="Input directory which contains fastq files, the fastq files shoule be anarchived and the extension should be .fq", metavar="INPUT")
    parser.add_option("-o", "--output",
                  dest="outputDirectry",
                  help="Directry for output, have to spply")
    parser.add_option("-d", "--hiCLIP_project_dir", dest="hiCLIPDir",
                  help="Path to hiCLIP_project")

    (options, args) = parser.parse_args()

    inputfile = options.inputDir
    outputdir = options.outputDirectry
    hiCLIPdir = options.hiCLIPDir

    runBowtieForRibosomeProfiling(inputfile, outputdir, hiCLIPdir)
