'''
Created on Aug 28, 2013

@author: yoichiro
'''
#!/bin/env python

# runBowtieForRibosomeProfiling.py
# The fq file created by this file should be categorised as RUT (Ribosome RNA unmapped and trimmed to 26 nucleotides)
#
# Created by Yoichiro Sugimoto on 14/08/2012.
# Copyright 2012 __MyCompanyName__. All rights reserved.
# Command example: qsub -j y -b y -pe smp 8 python ~/temp_run/runRibosomeProfilingAnalysis-4.0.py -i /netscr/yoichiro/RibosomeProfiling/rowReads/ -o /netscr/yoichiro/RibosomeProfiling/ -b /netscr/yoichiro/RibosomeProfiling/barcode/
# Command example: qsub -j y -b y -pe smp 8 python ~/temp_run/runRibosomeProfilingAnalysis-4.0.py -i /netscr/yoichiro/mRNASeq/rowReads/ -o /netscr/yoichiro/mRNASeq/ -b /netscr/yoichiro/mRNASeq/barcode/


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

def runBowtieForRibosomeProfiling(inputfile, outputdir, barcodedir, hiCLIPdir, pythondir, Rdir):

    if outputdir[-1] == "/":
        outputdir = outputdir[:-1]

    runUNIXCommands("mkdir -p " + outputdir + "/F")
    runUNIXCommands("mkdir -p " + outputdir + "/swappedBarcode")

    runUNIXCommands("mkdir -p " + outputdir + "/indexed")

    runUNIXCommands("mkdir -p " + outputdir + "/RU")
    runUNIXCommands("mkdir -p " + outputdir + "/RUU")
    runUNIXCommands("mkdir -p " + outputdir + "/RUT")

    runUNIXCommands("mkdir -p " + outputdir + "/RUTT")
    runUNIXCommands("mkdir -p " + outputdir + "/randomBarcode")

    runUNIXCommands("mkdir -p " + outputdir + "/TophatOut")
    runUNIXCommands("mkdir -p " + outputdir + "/sam")

    runUNIXCommands("mkdir -p " + outputdir + "/HTSeqCount")
    runUNIXCommands("mkdir -p " + outputdir + "/bed")

# Step 1, FASTQ to FASTA

    fileLists = listdir_fullpath_extension(inputfile, '.fq')

    for fq in fileLists:
        outName = outputdir + "/F/" + fq.split("/")[-1].split(".")[0] + "_F.fa"
        fqToFaCommand = hiCLIPdir + "/inst/bin/fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64/fastq_to_fasta -Q 33 -v -i " + fq + " | " + hiCLIPdir + "/inst/bin/fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64/fastx_renamer -n COUNT -o " + outName
        runUNIXCommands(fqToFaCommand)

# Step 2, swap barcode from random 3 barcode + 4 experiment barcode + 2 random barcode
#            to 4 experiment barcode + 5 random barcode

    faFileLists = listdir_fullpath_extension(outputdir + "/F", '.fa')

    for fa in faFileLists:
        outName = outputdir + "/swappedBarcode/" + fa.split("/")[-1].replace("_F", "_swappedBarcode")
        swapBarcodeCommand = pythondir + "/python " + hiCLIPdir + "/inst/Python/swapBarcodes.py -i " + fa + " -o " + outName
        runUNIXCommands(swapBarcodeCommand)

# Step 3, index library

    faFileLists = listdir_fullpath_extension(outputdir + "/swappedBarcode", 'fa')

    barcodeLists = listdir_fullpath_extension(barcodedir, '.barcode')

    for fa in faFileLists:
        indexedOutDir = outputdir + "/indexed/"
        for barcodes in barcodeLists:
            if barcodes.split("/")[-1].find(fa.split("/")[-1].replace("_swappedBarcode.fa", "")) != -1:
                matchedBarcodes = barcodes
            else:
                continue

        indexCommand = 'cat ' + fa +  ' | ' + hiCLIPdir + '/inst/bin/fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64/fastx_barcode_splitter.pl --bcfile ' + matchedBarcodes  + ' --bol --mismatches 1 --prefix ' + indexedOutDir + ' --suffix ".fa"'
        runUNIXCommands(indexCommand)

## From second run, start from here
# Step 4-Phase 1, map hybrid reads to rRNAs and tRNAs, and collect unmapped reads
    faFileLists = listdir_fullpath_extension(outputdir + "/indexed", '.fa')
    
    for fa in faFileLists:
        outName = outputdir + "/RUU/" + fa.split("/")[-1].replace(".fa", "_RUU.fa")
        bowtieCommand = hiCLIPdir + "/inst/bin/bowtie-0.12.7/bowtie -f -p 8 -v 2 --trim5 9 --trim3 15 --un " + outName + " " + hiCLIPdir + "/data/processed/bowtie_index/hs_rRNAs_and_tRNAs/hs_rRNAs_and_tRNAs " + fa + " > /dev/null"
        runUNIXCommands(bowtieCommand)  
        
# Step 4-Phase 2, map hybrid reads to mtDNA and pre-rRNAs
    faFileLists = listdir_fullpath_extension(outputdir + "/RUU", ".fa")
    
    for fa in faFileLists:
        outName = outputdir + "/RU/" + fa.split("/")[-1].replace("_RUU.fa", "_RU.fa")
        bowtieCommand = hiCLIPdir + "/inst/bin/bowtie-0.12.7/bowtie -f -p 8 -v 2 --trim5 9 --trim3 15 --un " + outName + " " + hiCLIPdir + "/data/processed/bowtie_index/chrM_pre_rRNA/chrM_pre_rRNA " + fa + " > /dev/null"
        runUNIXCommands(bowtieCommand)  


# Step 5, clip L3 linker and discard the RNA length is less than 26 (35 including barcodes)

    trimmedFileLists = listdir_fullpath_extension(outputdir + "/RU", '.fa')

    for fa in trimmedFileLists:
        outName = outputdir + "/RUT/" + fa.split("/")[-1].replace("_RU", "_RUT")
        trimL3Command = hiCLIPdir + "/inst/bin/fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64/fastx_clipper -a AGATCGGAAGAGCGGTTCAG -l 35 -v -n -i " + fa + " -o " + outName
        runUNIXCommands(trimL3Command)

# Step 6, trim barcodes and read length to 26
    RUFileLists = listdir_fullpath_extension(outputdir + "/RUT", '.fa')

    for fa in RUFileLists:
        outName = outputdir + "/RUTT/" + fa.split("/")[-1].replace("_RUT", "_RUTT")
        trim35Command = hiCLIPdir + "/inst/bin/fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64/fastx_trimmer -f 10 -l 35 -Q 33 -i " + fa + " -o " + outName
        runUNIXCommands(trim35Command)

# Step 6-2, get random barcodes
    RUFileLists = listdir_fullpath_extension(outputdir + "/RUT", '.fa')

    for fa in RUFileLists:
        outName = outputdir + "/randomBarcode/" + fa.split("/")[-1].replace("_RUT", "_randomBarcode")
        randomBarcodeCommand = hiCLIPdir + "/inst/bin/fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64/fastx_trimmer -f 5 -l 9 -i " + fa + " -o " + outName
        runUNIXCommands(randomBarcodeCommand)


# Step 7-Phase 3, map to transcripts by bowtie
    RUTTFileLists = listdir_fullpath_extension(outputdir + "/RUTT", '.fa')

    for fa in RUTTFileLists:
        outName = outputdir + "/sam/" + fa.split("/")[-1].replace("_RUTT.fa", ".sam")
        bowtieCommand = hiCLIPdir + "/inst/bin/bowtie-0.12.7/bowtie -f -p 8 -v 2 -m 1 --best --strata --sam " + hiCLIPdir + "/data/processed/bowtie_index/longest_mRNA_and_ncRNA/longest_mRNA_and_ncRNA " + fa + " > " + outName
        runUNIXCommands(bowtieCommand)


# Step 9, create BedGraph file and barcode consideration. If necessary, shift position according to the mapping results.
    samFileLists = listdir_fullpath_extension(outputdir + "/sam", '.sam')
    randomBarcodeLists = listdir_fullpath_extension(outputdir + "/randomBarcode", '.fa')

    for samFile in samFileLists:
        bedOutDir = outputdir + "/bed/"
        for randomBarcodes in randomBarcodeLists:
            if randomBarcodes.split("/")[-1].find(samFile.split("/")[-1].replace(".sam", "")) != -1:
                matchedRandomBarcodes = randomBarcodes
            else:
                continue

        bedCommand = pythondir + '/python ' + hiCLIPdir + '/inst/Python/samToBed-2.0.py ' + samFile +  ' ' + matchedRandomBarcodes  + ' ' + bedOutDir + ' 12'
        runUNIXCommands(bedCommand)



if __name__ == "__main__":
    from optparse import OptionParser

    parser = OptionParser()
    
    parser.add_option("-i", "--input", dest="inputFile",
                  help="Input directory which contains fastq files, the fastq files shoule be anarchived and the extension should be .fq", metavar="INPUT")
    parser.add_option("-o", "--output",
                  dest="outputDirectry",
                  help="Directry for output, have to spply")
    parser.add_option("-b", "--barcode",
                  dest="barcodeDirectory",
                  help="Directry for barcode lists, file have to have the same name as input fq file and .barcode extension")
    parser.add_option("-d", "--hiCLIP_project_dir", dest="hiCLIPDir",
                  help="Path to hiCLIP_project")
    parser.add_option("-p", "--python_dir", dest="pythonDir",
                  help="Path to Python")
    parser.add_option("-r", "--R_dir", dest="RDir",
                  help="Path to R")

    (options, args) = parser.parse_args()

    inputfile = options.inputFile
    outputdir = options.outputDirectry
    barcodedir = options.barcodeDirectory
    hiCLIPdir = options.hiCLIPDir
    pythondir = options.pythonDir
    Rdir = options.RDir
    
    runBowtieForRibosomeProfiling(inputfile, outputdir, barcodedir, hiCLIPdir, pythondir, Rdir)

