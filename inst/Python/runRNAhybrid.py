## Yoichiro Sugimoto 2012/Aug/17
## input shold be fasta file
## swap barcode from random 3 barcode + 4 experiment barcode + 2 random barcode 
#            to 4 experiment barcode + 5 random barcode

import os
import sys, subprocess

def runUNIXCommands(commands):
    # print commands
    try:
        retcode = subprocess.call(commands, shell=True)
        if retcode < 0:
            print >>sys.stderr, "Child was terminated by signal", -retcode
        else:
            print >>sys.stderr, "Child returned", retcode
    except OSError, e:
        print >>sys.stderr, "Execution failed:", e

def runRNAhybrid(inputfile, outputfile, hiCLIPdir):
    inFile = open(inputfile,'r')
    lineCount = 0

    allHybridOut = ""
    outFile = open(outputfile, "w")
    outFile.write(allHybridOut)
    outFile.close()

    for line in inFile:
        lineCount = lineCount + 1
        if lineCount % 2 == 0:
            read = line.strip()
            seqs = read.split("&")
            commandA = hiCLIPdir + "/inst/bin/RNAhybrid/bin/RNAhybrid -s 3utr_human -m 1000 -n 1000 -c " + seqs[0] + " " + seqs[1] + " >> " + outputfile
            runUNIXCommands(commandA)
    inFile.close()



if __name__ == "__main__":
    from optparse import OptionParser

    parser = OptionParser()

    parser.add_option("-i", "--input", dest="inputFile",
                  help="Input fasta file name, the experimental barcode shoule be located at [3:6] and random barcode at [0:2] and [7:8] on 0-base", metavar="INPUT")
    parser.add_option("-o", "--output",
                  dest="outputFile",
                  help="Output file name, have to spply")
    parser.add_option("-d", "--hiCLIP_project_dir", dest="hiCLIPDir",
                  help="Path to hiCLIP_project")


    (options, args) = parser.parse_args()

    inputfile = options.inputFile
    outputfile = options.outputFile
    hiCLIPdir = options.hiCLIPDir

    runRNAhybrid(inputfile, outputfile, hiCLIPdir)
