#!/bin/env python

# parse_RNAfold.py
#  qsub -b y -pe smp 8 -j y /lmb/home/yoichiro/Python-2.7.1/python ~/test/parse_RNAfold.py -i /netscr/yoichiro/STAU1_hiCLIP/results/manuscript/specificity/RNAfold/pos/result.db.txt -o /netscr/yoichiro/test -d /netscr/yoichiro/STAU1_hiCLIP

# Created by Yoichiro Sugimoto on 12/12/2013.
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

def runParseRNAfold(inputfile, outputdir, hiCLIPdir):
    final_output = outputdir + "/RNAforgi_result.txt"
    outputdir_temp = outputdir + "/temp"
    outputdir_temp2 = outputdir + "/temp2"
    mkdir_temp = "mkdir -p " + outputdir_temp
    mkdir_temp2 = "mkdir -p " + outputdir_temp2
    runUNIXCommands(mkdir_temp)
    runUNIXCommands(mkdir_temp2)

    os.chdir(outputdir_temp)

    initialize_file = "rm -f " + final_output
    runUNIXCommands(initialize_file)

    lineCount = 0

    inFile = open(inputfile,'r')

    for line in inFile:
        line = line.strip()
        lineCount = lineCount + 1
        if lineCount % 6 == 1:
            read_id = line
            print read_id
        elif lineCount % 6 == 5:
            db = line
            db = db.split(" {")[0]
            print db

            temp_fafile = outputdir + "/temp2/temp.fa"
            tempFaFile = open(temp_fafile,'w')
            tempFaFile.write(db + "\n")
            tempFaFile.close()

            id_command = "echo '" + read_id + "' >> " + final_output
            # RNAforgi_command = "echo '" + db + "' | python /lmb/home/yoichiro/Python/forgi/examples/longest_stem.py - " + " >> " + final_output
            RNAforgi_command = "python /lmb/home/yoichiro/Python/forgi/examples/dotbracket_to_bulge_graph.py " + temp_fafile + " >> " + final_output

            runUNIXCommands(id_command)
            runUNIXCommands(RNAforgi_command)

    inFile.close()
    # file.close()


if __name__ == "__main__":
    from optparse import OptionParser

    parser = OptionParser()

    parser.add_option("-i", "--input", dest="inputFile",
                  help="Input file", metavar="INPUT")
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

    runParseRNAfold(inputfile, outputdir, hiCLIPdir)
