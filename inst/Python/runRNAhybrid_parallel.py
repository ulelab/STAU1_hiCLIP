## Yoichiro Sugimoto 2012/Aug/17
## test command: qsub -j y -b y -pe smp 8 ~/Python-2.7.1/python /netscr/yoichiro/STAU1_hiCLIP/inst/Python/runRNAhybrid_parallel.py -i /netscr/yoichiro/STAU1_hiCLIP/results/manuscript/FDR_all/intra.inter.gene.cofold.fa -o /netscr/yoichiro/STAU1_hiCLIP/results/manuscript/FDR_all/intra.inter.gene.cofold.test.RNAhybrid.out -d /netscr/yoichiro/STAU1_hiCLIP
##

import os, sys, subprocess, multiprocessing, commands


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

def executeRNAhybrid(inputseq):
    seqs = inputseq.split("&")
    commandA = hiCLIPdir + "/inst/bin/RNAhybrid/bin/RNAhybrid -s 3utr_human -m 1000 -n 1000 -c " + seqs[0] + " " + seqs[1]
    RNAhybrid_out = commands.getstatusoutput(commandA)[1]
    return RNAhybrid_out

def runRNAhybrid(inputfile, outputfile, hiCLIPdir):
    inFile = open(inputfile,'r')
    lineCount = 0

    inputseq_list = []

    # This block reads input sequences and create the list
    for line in inFile:
        lineCount = lineCount + 1
        if lineCount % 2 == 0:
            read = line.strip()
            inputseq_list.append(read)

    inFile.close()

    cpu_num = int(multiprocessing.cpu_count())
    print "The number of CPU used is: " + str(cpu_num)

    pall_out = Parallel(n_jobs=cpu_num)(delayed(executeRNAhybrid)(i) for i in inputseq_list)
    pall_out = "\n".join(pall_out) + "\n"

    outFile = open(outputfile, "w")
    outFile.write(pall_out)
    outFile.close()

    return 0

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

    joblibPATH = hiCLIPdir + "/inst/bin/joblib-0.8.2-py2.7.egg"
    sys.path.append(joblibPATH)
    from joblib import Parallel, delayed

    runRNAhybrid(inputfile, outputfile, hiCLIPdir)
