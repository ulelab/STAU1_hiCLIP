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

def runFaPreparation(hiCLIPdir):
    command0 = "mkdir " + hiCLIPdir + "/results/mapping/hiCLIP/nonHybrid/RPlike/Fasta_26_merged"
    command1 = "cat " + hiCLIPdir + "/results/mapping/hiCLIP/nonHybrid/RPlike/Fasta_26/nonHybrid_DOX_LigMinus.26.fa " + hiCLIPdir + "/results/mapping/hiCLIP/nonHybrid/RPlike/Fasta_26/nonHybrid_DOX_LigPlusHighRNase.26.fa " + hiCLIPdir + "/results/mapping/hiCLIP/nonHybrid/RPlike/Fasta_26/nonHybrid_DOX_LigPlus.26.fa > " + hiCLIPdir + "/results/mapping/hiCLIP/nonHybrid/RPlike/Fasta_26_merged/hiCLIP_nonHybrid_merged.26.fa"
    command2 = "mkdir /netscr/yoichiro/STAU1_hiCLIP/results/mapping/mRNASeq/merged_fa"
    command3 = "cat /netscr/yoichiro/STAU1_hiCLIP/results/mapping/mRNASeq/RUTT/id1_wt_1_mRNASeq_RUTT.fa /netscr/yoichiro/STAU1_hiCLIP/results/mapping/mRNASeq/RUTT/id2_wt_2_mRNASeq_RUTT.fa > /netscr/yoichiro/STAU1_hiCLIP/results/mapping/mRNASeq/merged_fa/mRNASeq_merged.26.fa"
    command4 = "cat /netscr/yoichiro/STAU1_hiCLIP/results/mapping/STAU1_iCLIP/RUTT/id1_STAU1iCLIP_RUTT.fa /netscr/yoichiro/STAU1_hiCLIP/results/mapping/STAU1_iCLIP/RUTT/id2_STAU1iCLIP_RUTT.fa /netscr/yoichiro/STAU1_hiCLIP/results/mapping/STAU1_iCLIP/RUTT/id3_STAU1iCLIP_RUTT.fa /netscr/yoichiro/STAU1_hiCLIP/results/mapping/STAU1_iCLIP/RUTT/id4_STAU1iCLIP_RUTT.fa /netscr/yoichiro/STAU1_hiCLIP/results/mapping/STAU1_iCLIP/RUTT/id5_STAU1iCLIP_RUTT.fa /netscr/yoichiro/STAU1_hiCLIP/results/mapping/STAU1_iCLIP/RUTT/id6_STAU1iCLIP_RUTT.fa > /netscr/yoichiro/STAU1_hiCLIP/results/mapping/STAU1_iCLIP/RUTT/maerged_STAU1_iCLIP.fa"
    command5 = "gzip /netscr/yoichiro/STAU1_hiCLIP/results/mapping/hiCLIP/nonHybrid/RPlike/Fasta_26_merged/hiCLIP_nonHybrid_merged.26.fa"
    command6 = "gzip /netscr/yoichiro/STAU1_hiCLIP/results/mapping/mRNASeq/merged_fa/mRNASeq_merged.26.fa"
    command7 = "gzip /netscr/yoichiro/STAU1_hiCLIP/results/mapping/hnRNPc/RUTT/id1_hnRNPc_RUTT.fa"
    command8 = "gzip /netscr/yoichiro/STAU1_hiCLIP/results/mapping/STAU1_iCLIP/RUTT/maerged_STAU1_iCLIP.fa"

    runUNIXCommands(command0)
    runUNIXCommands(command1)
    runUNIXCommands(command2)
    runUNIXCommands(command3)
    runUNIXCommands(command4)
    runUNIXCommands(command5)
    runUNIXCommands(command6)
    runUNIXCommands(command7)
    runUNIXCommands(command8)


if __name__ == "__main__":
    from optparse import OptionParser

    parser = OptionParser()

    parser.add_option("-d", "--hiCLIP_project_dir", dest="hiCLIPDir",
                  help="Path to hiCLIP_project")


    (options, args) = parser.parse_args()

    hiCLIPdir = options.hiCLIPDir

    runFaPreparation(hiCLIPdir)
