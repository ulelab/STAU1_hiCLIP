## Yoichiro: remove reads with (2 or more Link) or (LinkL3) in the sequence
## Modified from a code in "http://www.biostars.org/post/show/15446/demultiplex-illumina-with-barcodes-on-identifier-line/"
## Since original Levenshtein function was too slow, I introduced cutadapt packages, which is fater but with a bit issues.

import sys
import subprocess

def runUNIXCommands(commands):
    try:
        retcode = subprocess.call(commands, shell=True)
        if retcode < 0:
            print >>sys.stderr, "Child was terminated by signal", -retcode
        else:
            print >>sys.stderr, "Child returned", retcode
    except OSError, e:
        print >>sys.stderr, "Execution failed:", e

def calcHamming(seq, adaptor):
    return [Levenshtein.hamming(m, adaptor) for m in [seq[n:(n+len(adaptor))] for n in range(len(seq) - (len(adaptor) - 1))]]

def minHamming(seq, adaptor):
    return min([Levenshtein.hamming(m, adaptor) for m in [seq[n:(n+len(adaptor))] for n in range(len(seq) - (len(adaptor) - 1))]])

def trimLinkL3andL3(inputfile, outputName, minLength, LinkNum):
    # This script can be used to select hybrid candidate (minLength = 53 (= 17 + 19 + 17)  & LinkNum = 1) or select reads for non-hybrid (minLength = 53 & LinkNum = 0)

    minLength = int(minLength)
    LinkNum = int(LinkNum)
    
    inFile = open(inputfile,'r')
    inputName = inputfile.split("/")[-1]

    header = ''
    data = ''

    outFile = {}

    Link = "CTGTAGGCACCATACAATG"
    
    lineCount = 0
    hybridCount = 0
    
    outFile = open(outputName, "w")

    for line in inFile:
        AdaptorFlag = False

        lineCount = lineCount + 1

        if lineCount % 2 == 1:
            rname = line.strip()
        else:
            seq = line.strip()
                
            # require to have L3 in reads and allow 2 mismatch for LinkL3 and 1 mismatc for L3

            if len(seq) >= minLength:
                if (sum(i <= 1 for i in calcHamming(seq, Link)) == LinkNum):
                    hybridCount = hybridCount + 1
                    # After adaptor removal, the reads with multiple Link is likely to be low quality reads...
                    outFile.write(rname + "\n" + seq + "\n")
                
    print "The number of total reads processed..."
    print str(lineCount/2)

    inFile.close()
    outFile.close()


def selectMappableHybrid(inputfile, outputName, space):

    # space indicate the mergin of the reads from the end to Link adaptor.
    # longer space help mapping, but less likely to occur

    space = int(space)

    inFile = open(inputfile,'r')

    Link = "CTGTAGGCACCATACAATG"
    
    outFile = open(outputName, "w")
    mappableHybridCount = 0
    lineCount = 0

    for line in inFile:
        lineCount = lineCount + 1

        if lineCount % 2 == 1:
            rname = line.strip()
        else:
            seq = line.strip()

            # Allow two mismatch
            LinkPos = [i for i,x in enumerate(calcHamming(seq, Link)) if x == minHamming(seq, Link)][0]

            if (LinkPos >= space) and (LinkPos <= len(seq) - space - len(Link)):
                mappableHybridCount = mappableHybridCount + 1
                outFile.write(rname + "_L" + "\n" + seq[:LinkPos] + "\n" + rname + "_R" + "\n" + seq[LinkPos + len(Link):] + "\n")

    inFile.close()
    outFile.close()
    print "total reads processed: " + str(lineCount/2)
    print "The number of mappable hybrid reads: " + str(mappableHybridCount)

if __name__ == "__main__":
    import sys
    from optparse import OptionParser

    parser = OptionParser()
    
    parser.add_option("-d", "--hiCLIP_project_dir", dest="hiCLIPDir",
                  help="Path to hiCLIP_project")
    parser.add_option("-i", "--input", dest="inputFile",
                  help="Input fasta file name", metavar="INPUT")
    parser.add_option("-o", "--output",
                  dest="outputDirectry",
                  help="Directry for output, have to spply")

    (options, args) = parser.parse_args()

    hiCLIPdir = options.hiCLIPDir
    inputfile = options.inputFile
    outputdir = options.outputDirectry

    LevPATH = hiCLIPdir + "/inst/bin/python_Levenshtein-0.10.2-py2.7-linux-x86_64.egg"
    sys.path.append(LevPATH)
    import Levenshtein
    
    if outputdir[-1] == "/":
        outputdir = outputdir[:-1]
    inputName = inputfile.split("/")[-1]
    
    # Creating directories
    runUNIXCommands("mkdir -p " + outputdir + "/hybrid/")
    runUNIXCommands("mkdir -p " + outputdir + "/hybrid/hybridTemp/")
    runUNIXCommands("mkdir -p " + outputdir + "/hybrid/hybridTemp/I3T")
    runUNIXCommands("mkdir -p " + outputdir + "/hybrid/hybridTemp/nI3T")
    runUNIXCommands("mkdir -p " + outputdir + "/hybrid/hybridTemp/L3T")
    runUNIXCommands("mkdir -p " + outputdir + "/hybrid/hybridTemp/merged")
    
    runUNIXCommands("mkdir -p " + outputdir + "/hybrid/hybridTemp/AT")
    
    runUNIXCommands("mkdir -p " + outputdir + "/hybrid/Fasta/")
    
    runUNIXCommands("mkdir -p " + outputdir + "/nonHybrid/")
    runUNIXCommands("mkdir -p " + outputdir + "/nonHybrid/Fasta")
    
    # For both
    #    1. trim LinkL3 linker which occures ~18% first and create output I3T/--_I3T.fa in which LinkL3 was trimmed.
    #    2. Then, for the reads which does not have LinkL3 (nI3T/--_nI3T.fa), trim L3 reads (output L3T/--_L3T.fa) (~ 72% of total reads)
    #    3. merge I3T and L3T files as full length read cDNA set (merged/_merged.fa).
    

    # 1
    I3TOut =  outputdir + "/hybrid/hybridTemp/I3T/" + inputfile.split("/")[-1].replace(".fa", "_I3T.fa")
    nI3TOut =  outputdir + "/hybrid/hybridTemp/nI3T/" + inputfile.split("/")[-1].replace(".fa", "_nI3T.fa")
    
    I3Tcommand = hiCLIPdir + "/inst/bin/cutadapt-1.2.1/bin/cutadapt -a CTGTAGGCACCATACAATGAGATCGGAAGAGCGGTTCAG -e 0.06 -n 10 -m 186 --too-short-output=" + I3TOut + " " + inputfile + " > " + nI3TOut 
    print(I3Tcommand)
    runUNIXCommands(I3Tcommand)
    
    # 2
    L3TOut = outputdir + "/hybrid/hybridTemp/L3T/" + inputfile.split("/")[-1].replace(".fa", "_L3T.fa")
    L3Tcommand = hiCLIPdir + "/inst/bin/cutadapt-1.2.1/bin/cutadapt -a AGATCGGAAGAGCGGTTCAG -e 0.06 -n 10 -O 10 -m 186 --too-short-output=" + L3TOut + " " + nI3TOut + " > /dev/null"
    print(L3Tcommand)
    runUNIXCommands(L3Tcommand)
    
    # 3
    mergedOut = outputdir + "/hybrid/hybridTemp/merged/" + inputfile.split("/")[-1].replace(".fa", "_merged.fa")
    mergedCommand = "cat " + L3TOut + " " + I3TOut + " > " + mergedOut
    print(mergedCommand)
    runUNIXCommands(mergedCommand)
      
    # For hybrid only
    hybridTempOut = outputdir + "/hybrid/hybridTemp/AT/AT_" + inputName
    trimLinkL3andL3(mergedOut, hybridTempOut, minLength = 53, LinkNum = 1)
    
    hybridOut = outputdir + "/hybrid/Fasta/hybrid_" + inputName
    selectMappableHybrid(hybridTempOut, hybridOut, space = 17)
    
    # For non hybrid only
    nonHybridOut = outputdir + "/nonHybrid/Fasta/nonHybrid_" + inputName
    trimLinkL3andL3(mergedOut, nonHybridOut, minLength = 53, LinkNum = 0)
    
