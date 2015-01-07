# Yoichiro Sugimoto

import sys, os, subprocess

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

def formatForRNAcofold(input, output):
    file = open(output, "w")

    lineCount = 0

    inFile = open(input,'r')
    for line in inFile:
        line = line.strip()
        lineCount = lineCount + 1
        if lineCount % 4 == 1:
            L_id = line.strip(">")
        elif lineCount % 4 == 2:
            L_sread = line
        elif lineCount % 4 == 3:
            R_id = line.strip(">")
        elif lineCount % 4 == 0:
            R_sread = line
            if L_id.strip("_L") != R_id.strip("_R"):
                print L_id + "\t" + R_id
                raise NameError('The fasta file seems not to be sorted!!')

            readId = ">" + L_id.strip("_L")
            sread = L_sread + "&" + R_sread
            hybridOut = readId + "\n" + sread + "\n"
            file.write(hybridOut)
    inFile.close()
    file.close()


def RNAcofoldAnalysis(hybridInput, nonHybridInput, outputdir, hiCLIPdir, Rdir):

    # 1. Format FASTA file to RNAcofold format (cofoldFa)
    if outputdir[-1] == "/":
        outputdir = outputdir[:-1]

    runUNIXCommands("mkdir -p " + outputdir + "/RNAcofold/cofoldFa")

    hybridFaOutput = outputdir + "/RNAcofold/cofoldFa/Hybrid_" + hybridInput.split("/")[-1].split(".fa")[0] + ".cofoldFa"
    nonHybridFaOutput = outputdir + "/RNAcofold/cofoldFa/nonHybrid_" + nonHybridInput.split("/")[-1].split(".fa")[0] + ".cofoldFa"

    formatForRNAcofold(hybridInput, hybridFaOutput)
    formatForRNAcofold(nonHybridInput, nonHybridFaOutput)

    # 2. RNAcofold analysis
    hybridRNAcofoldOut = outputdir + "/RNAcofold/cofoldOut/Hybrid_" + hybridInput.split("/")[-1].split(".fa")[0] + ".cofoldOut"
    nonHybridRNAcofoldOut = outputdir + "/RNAcofold/cofoldOut/nonHybrid_" + nonHybridInput.split("/")[-1].split(".fa")[0] + ".cofoldOut"


    runUNIXCommands("mkdir -p " + outputdir + "/RNAcofold/cofoldOut")
    commandRNAcofoldHybrid = hiCLIPdir + "/inst/bin/ViennaRNA/Progs/RNAcofold -noPS < " + hybridFaOutput + " > " + hybridRNAcofoldOut
    commandRNAcofoldNonHybrid = hiCLIPdir + "/inst/bin/ViennaRNA/Progs/RNAcofold -noPS < " + nonHybridFaOutput + " > " + nonHybridRNAcofoldOut
    runUNIXCommands(commandRNAcofoldHybrid)
    runUNIXCommands(commandRNAcofoldNonHybrid)

    # 3. Density plot and two-sample Wilcoxon test with R
#    plotFileName = hybridInput.split("/")[-1].split(".fa")[0] + "vs" + nonHybridInput.split("/")[-1].split(".fa")[0] + ".pdf"
#    Rcommand = Rdir + "/bin/Rscript " + hiCLIPdir + "/R/RNAcofold/plotRNAcofoldOut.R " + outputdir + "/RNAcofold/" + " "  + hybridRNAcofoldOut + " " + nonHybridRNAcofoldOut + " " + plotFileName + " " + Rdir
#    runUNIXCommands(Rcommand)


if __name__ == "__main__":
    from optparse import OptionParser

    parser = OptionParser()

    parser.add_option("-i", "--hybrid", dest="hybridInput",
                  help="Hybrid sequence fasta file should be separated by Link and sorted by read name", metavar="INPUT")
    parser.add_option("-n", "--nonhybrid", dest="nonHybridInput",
                  help="Non-Hybrid sequence fasta file should be separated by simulated size and sorted by read name", metavar="INPUT")
    parser.add_option("-o", "--output",
                  dest="outputDirectry",
                  help="Directry for output, have to spply")
    
    parser.add_option("-d", "--hiCLIP_project_dir", dest="hiCLIPDir",
                  help="Path to hiCLIP_project")
    parser.add_option("-r", "--R_dir", dest="RDir",
                  help="Path to R")

    (options, args) = parser.parse_args()

    hybridInput = options.hybridInput
    nonHybridInput = options.nonHybridInput
    outputdir = options.outputDirectry
    hiCLIPdir = options.hiCLIPDir
    Rdir = options.RDir

    RNAcofoldAnalysis(hybridInput, nonHybridInput, outputdir, hiCLIPdir, Rdir)