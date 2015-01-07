## Yoichiro Sugimoto 2012/Aug/17
## test command: qsub -j y -b y -pe smp 2 ~/Python-2.7.1/python /netscr/yoichiro/STAU1_hiCLIP/inst/Python/runRNAforgi_parallel.py -i /lmb/home/yoichiro/test.db.txt -o /netscr/yoichiro/STAU1_hiCLIP/results/manuscript/RNAfold/ -d /netscr/yoichiro/STAU1_hiCLIP -p ~/Python-2.7.1
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

def executeRNAforgi(forgi_outfile, pythondir):
    commandA = pythondir + "/python " + hiCLIPdir + "/inst/bin/forgi/examples/dotbracket_to_bulge_graph.py " + forgi_outfile
    RNAhybrid_out = commands.getstatusoutput(commandA)[1]
    return ">" + forgi_outfile.split("/")[-1].split(".")[0] + "\n" + RNAhybrid_out

def runRNAhybrid(inputfile, outputfile, hiCLIPdir, pythondir):
    inFile = open(inputfile,'r')
    lineCount = 0

    outfile_list = []

    temp_file = outputfile + "/temp_forgi"
    runUNIXCommands("mkdir -p " + temp_file)

    for line in inFile:
        lineCount = lineCount + 1
        if lineCount % 3 == 1:
            read_id = line.strip()
            outfile_name = temp_file + "/" + read_id.replace(">", "") + ".txt"
            outfile_list.append(outfile_name)
        if lineCount % 3 == 0:
            read = line.strip()
            read = read.split(" ")[0]
            outFile = open(outfile_name, "w")
            outFile.write(read)
            outFile.close()


    inFile.close()

    cpu_num = int(multiprocessing.cpu_count())
    print "The number of CPU used is: " + str(cpu_num)

    pall_out = Parallel(n_jobs=cpu_num)(delayed(executeRNAforgi)(i, j) for i, j in zip(outfile_list, [pythondir] * len(outfile_list)))
    pall_out = "\n".join(pall_out) + "\n"

    forgi_out = outputfile + "/RNAforgi.out.txt"

    outForgiFile = open(forgi_out, "w")
    outForgiFile.write(pall_out)
    outForgiFile.close()

    return 0

if __name__ == "__main__":
    from optparse import OptionParser

    parser = OptionParser()

    parser.add_option("-i", "--input", dest="inputFile",
                  help="Input fasta file name, the experimental barcode shoule be located at [3:6] and random barcode at [0:2] and [7:8] on 0-base", metavar="INPUT")
    parser.add_option("-o", "--output",
                  dest="outputFile",
                  help="Output file name, have to spply")
    parser.add_option("-p", "--python_dir", dest="pythonDir",
                  help="Path to Python")
    parser.add_option("-d", "--hiCLIP_project_dir", dest="hiCLIPDir",
                  help="Path to hiCLIP_project")



    (options, args) = parser.parse_args()

    inputfile = options.inputFile
    outputfile = options.outputFile
    pythondir = options.pythonDir
    hiCLIPdir = options.hiCLIPDir

    joblibPATH = hiCLIPdir + "/inst/bin/joblib-0.8.2-py2.7.egg"
    sys.path.append(joblibPATH)
    from joblib import Parallel, delayed

    runRNAhybrid(inputfile, outputfile, hiCLIPdir, pythondir)
