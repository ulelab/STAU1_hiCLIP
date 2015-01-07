'''
Created on May 20, 2013

@author: yoichiro
'''
#!/bin/env python
#
# 
# qsub -j y -b y -pe smp 2 ~/Python-2.7.1/python /netscr/yoichiro/STAU1_hiCLIP/exec/run_analysis_part.py -d /netscr/yoichiro/STAU1_hiCLIP -p ~/Python-2.7.1 -r ~/R-2.15.1
# qsub -j y -b y -pe smp 2 ~/Python-2.7.1/python /netscr/yoichiro/STAU1_hiCLIP2/exec/run_analysis_part.py -d /netscr/yoichiro/STAU1_hiCLIP2 -p ~/Python-2.7.1 -r ~/R-2.15.1


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
        
def runSweave(sweavefile, hiCLIPdir, Rdir, pythondir):
    Stangle_command = Rdir + "/bin/R -e \"Stangle('" + sweavefile + ".Rnw')\" --args hiCLIPdir=" + hiCLIPdir + " pythondir=" + pythondir + " Rdir=" + Rdir
    sweave_command = Rdir + "/bin/R -e \"Sweave('" + sweavefile + ".Rnw')\" --args hiCLIPdir=" + hiCLIPdir + " pythondir=" + pythondir + " Rdir=" + Rdir
    pdflatex_command = Rdir + "/bin/R CMD pdflatex " + sweavefile + ".tex"
    runUNIXCommands(Stangle_command)
    runUNIXCommands(sweave_command)
    runUNIXCommands(pdflatex_command)
    runUNIXCommands(pdflatex_command)
    runUNIXCommands(pdflatex_command)
    
def run_all(hiCLIPdir, pythondir, Rdir):
    # Reset all, comment out for test run
    # runUNIXCommands("rm -rf " + hiCLIPdir + "/results/*")
    #
    # runUNIXCommands("mkdir -p " + hiCLIPdir + "/results")
    # runUNIXCommands("mkdir -p " + hiCLIPdir + "/results/mapping")
    # runUNIXCommands("mkdir -p " + hiCLIPdir + "/results/mapping/DESeq")
    # runUNIXCommands("mkdir -p " + hiCLIPdir + "/results/mapping/hiCLIP")
    # runUNIXCommands("mkdir -p " + hiCLIPdir + "/results/mapping/mRNASeq")
    # runUNIXCommands("mkdir -p " + hiCLIPdir + "/results/mapping/hnRNPc")
    # runUNIXCommands("mkdir -p " + hiCLIPdir + "/results/mapping/ribosome_profiling")
    # runUNIXCommands("mkdir -p " + hiCLIPdir + "/results/manuscript")
    # runUNIXCommands("mkdir -p " + hiCLIPdir + "/results/manuscript/figures")
    # runUNIXCommands("mkdir -p " + hiCLIPdir + "/results/manuscript/IGV")
    # runUNIXCommands("mkdir -p " + hiCLIPdir + "/results/manuscript/range")
    # runUNIXCommands("mkdir -p " + hiCLIPdir + "/results/manuscript/RNAcoFold")
    # runUNIXCommands("mkdir -p " + hiCLIPdir + "/results/mapping/STAU1_iCLIP")

    os.chdir(hiCLIPdir + "/inst/doc")

    # runSweave("data_preprocessing_1", hiCLIPdir, Rdir, pythondir)
    # runSweave("data_preprocessing_2", hiCLIPdir, Rdir, pythondir)
    # runSweave("data_preprocessing_3", hiCLIPdir, Rdir, pythondir)
    #
    # i1_command = pythondir + "/python " + hiCLIPdir + "/exec/i1_map_HTSeq.py -d " + hiCLIPdir + " -p " + pythondir + " -r " + Rdir
    # runUNIXCommands(i1_command)


    # runSweave("hiCLIP_data_analysis_1", hiCLIPdir, Rdir, pythondir)
    # runSweave("hiCLIP_data_analysis_2", hiCLIPdir, Rdir, pythondir)
    # runSweave("hiCLIP_data_analysis_3", hiCLIPdir, Rdir, pythondir)
    #
    # fa_preparation_command = pythondir + "/python " + hiCLIPdir + "/inst/Python/repeatEnrichmentFa.py -d " + hiCLIPdir
    # runUNIXCommands(fa_preparation_command)
    #
    # runSweave("hiCLIP_data_analysis_4", hiCLIPdir, Rdir, pythondir)
    runSweave("hiCLIP_data_analysis_5", hiCLIPdir, Rdir, pythondir)
    # runSweave("hiCLIP_data_analysis_6", hiCLIPdir, Rdir, pythondir)
    # runSweave("hiCLIP_data_analysis_7", hiCLIPdir, Rdir, pythondir)
    # runSweave("hiCLIP_data_analysis_8", hiCLIPdir, Rdir, pythondir)
    # runSweave("hiCLIP_data_analysis_9", hiCLIPdir, Rdir, pythondir)


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

    run_all(hiCLIPdir, pythondir, Rdir)
