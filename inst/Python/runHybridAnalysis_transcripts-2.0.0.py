#!/bin/env python

# runHybridAnalysis_transcripts-2.0.0.py
# Strategt was summarised in /Users/yoichiro/Documents/STAU1_project/transcripts_mapping/Method.docx
# The fq file created by this file should be categorised as RUT (Ribosome RNA unmapped and trimmed to 26 nucleotides)
#
#  qsub -b y -pe smp 8 -j y /lmb/home/yoichiro/Python-2.7.1/python /lmb/home/yoichiro/temp_run/runHybridAnalysis_transcripts-2.0.0.py -i /netscr/yoichiro/hiCLIP/rowReads/ -o /netscr/yoichiro/hiCLIP/ -b /netscr/yoichiro/hiCLIP/barcode/
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

def runBowtieForRibosomeProfiling(inputdir, outputdir, barcodedir, hiCLIPdir, pythondir, Rdir):

    if outputdir[-1] == "/":
       outputdir = outputdir[:-1]

    runUNIXCommands("mkdir -p " + outputdir + "/F")
    runUNIXCommands("mkdir -p " + outputdir + "/swappedBarcode")

    runUNIXCommands("mkdir -p " + outputdir + "/indexed")

    runUNIXCommands("mkdir -p " + outputdir + "/TempT")
    runUNIXCommands("mkdir -p " + outputdir + "/T")
    runUNIXCommands("mkdir -p " + outputdir + "/randomBarcode")
    
    runUNIXCommands("mkdir -p " + outputdir + "/hybrid/")
    runUNIXCommands("mkdir -p " + outputdir + "/hybrid/sam/")
    runUNIXCommands("mkdir -p " + outputdir + "/hybrid/sam_P1_mapped")
    runUNIXCommands("mkdir -p " + outputdir + "/hybrid/sam_P1_mapped_unique")
    runUNIXCommands("mkdir -p " + outputdir + "/hybrid/sam_P3_mapped")
    runUNIXCommands("mkdir -p " + outputdir + "/hybrid/sam_P4_mapped")
    runUNIXCommands("mkdir -p " + outputdir + "/hybrid/fasta_P1_unmapped")
    runUNIXCommands("mkdir -p " + outputdir + "/hybrid/fasta_P2_unmapped")
    runUNIXCommands("mkdir -p " + outputdir + "/hybrid/fasta_P3_unmapped")
    runUNIXCommands("mkdir -p " + outputdir + "/hybrid/bed")
    runUNIXCommands("mkdir -p " + outputdir + "/hybrid/rRNA_unique_bed")
    
    runUNIXCommands("mkdir -p " + outputdir + "/nonHybrid/")
    runUNIXCommands("mkdir -p " + outputdir + "/nonHybrid/sam")
    runUNIXCommands("mkdir -p " + outputdir + "/nonHybrid/sam_P1_mapped")
    runUNIXCommands("mkdir -p " + outputdir + "/nonHybrid/sam_P3_mapped")
    runUNIXCommands("mkdir -p " + outputdir + "/nonHybrid/sam_P4_mapped")
    runUNIXCommands("mkdir -p " + outputdir + "/nonHybrid/fasta_P1_unmapped")
    runUNIXCommands("mkdir -p " + outputdir + "/nonHybrid/fasta_P2_unmapped")
    runUNIXCommands("mkdir -p " + outputdir + "/nonHybrid/fasta_P3_unmapped")
    runUNIXCommands("mkdir -p " + outputdir + "/nonHybrid/bedLike")
    runUNIXCommands("mkdir -p " + outputdir + "/nonHybrid/bed")
    
    runUNIXCommands("mkdir -p " + outputdir + "/nonHybrid/RPlike")
    runUNIXCommands("mkdir -p " + outputdir + "/nonHybrid/RPlike/RD")
    runUNIXCommands("mkdir -p " + outputdir + "/nonHybrid/RPlike/Fasta_26")
    runUNIXCommands("mkdir -p " + outputdir + "/nonHybrid/RPlike/TophatOut")
    runUNIXCommands("mkdir -p " + outputdir + "/nonHybrid/RPlike/sam/")
    runUNIXCommands("mkdir -p " + outputdir + "/nonHybrid/RPlike/bed/")
        

# Step 1, FASTQ to FASTA

    fileLists = listdir_fullpath(inputdir)
    fileLists = [each for each in fileLists if each.endswith('.fq')]

    for fq in fileLists:
        tempOutName = outputdir + "/F/" + fq.split("/")[-1].split(".")[0] + "_F.fa.temp"
        outName = outputdir + "/F/" + fq.split("/")[-1].split(".")[0] + "_F.fa"
        
        fqToFaCommand = hiCLIPdir + "/inst/bin/fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64/fastq_to_fasta -Q 33 -v -n -i " + fq + " -o " + tempOutName
        runUNIXCommands(fqToFaCommand)
        
        faRenameCommand = hiCLIPdir + "/inst/bin/fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64/fastx_renamer -n COUNT -i " + tempOutName + " -o " + outName
        runUNIXCommands(faRenameCommand)
        
        removeTempCommand = "rm " + tempOutName
        runUNIXCommands(removeTempCommand)
        

# Step 2, swap barcode from random 3 barcode + 4 experiment barcode + 2 random barcode 
#            to 4 experiment barcode + 5 random barcode

    faFileLists = listdir_fullpath(outputdir + "/F")
    faFileLists = [each for each in faFileLists if each.endswith('.fa')]

    for fa in faFileLists:
        outName = outputdir + "/swappedBarcode/" + fa.split("/")[-1].replace("_F", "_swappedBarcode")
        swapBarcodeCommand = pythondir + "/python " + hiCLIPdir + "/inst/Python/swapBarcodes.py -i " + fa + " -o " + outName
        runUNIXCommands(swapBarcodeCommand)

# Step 3, index library

    faFileLists = listdir_fullpath(outputdir + "/swappedBarcode")
    faFileLists = [each for each in faFileLists if each.endswith('.fa')]

    barcodeLists = listdir_fullpath(barcodedir)
    barcodeLists = [each for each in barcodeLists if each.endswith('.barcode')]

    for fa in faFileLists:
        indexedOutDir = outputdir + "/indexed/"
        for barcodes in barcodeLists:
            if barcodes.split("/")[-1].find(fa.split("/")[-1].replace("_swappedBarcode.fa", "")) != -1:
                matchedBarcodes = barcodes
            else:
                continue

        indexCommand = 'cat ' + fa +  ' | ' + hiCLIPdir + '/inst/bin/fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64/fastx_barcode_splitter.pl --bcfile ' + matchedBarcodes  + ' --bol --mismatches 1 --prefix ' + indexedOutDir + ' --suffix ".fa"'
        runUNIXCommands(indexCommand)

# Step 4, trim experimental barcodes and extract random barcode
    faFileLists = listdir_fullpath(outputdir + "/indexed")
    faFileLists = [each for each in faFileLists if each.endswith('.fa')]
    
    for fa in faFileLists:
        tempOutName = outputdir + "/TempT/" + fa.split("/")[-1].replace(".fa", "_TempT.fa")
        OutName = outputdir + "/T/" + fa.split("/")[-1].replace(".fa", "_T.fa")
        randomBarcodeOut = outputdir + "/randomBarcode/" + fa.split("/")[-1].replace(".fa", "_randomBarcode.fa")
        
        trimExpBarcodeCommand = hiCLIPdir + "/inst/bin/fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64/fastx_trimmer -f 5 -i " + fa + " -o " + tempOutName
        runUNIXCommands(trimExpBarcodeCommand)
        
        trimRandomBarcodeCommand = hiCLIPdir + "/inst/bin/fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64/fastx_trimmer -f 6 -i " + tempOutName + " -o " + OutName
        runUNIXCommands(trimRandomBarcodeCommand)
        
        extractRandomBarcodeCommand = hiCLIPdir + "/inst/bin/fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64/fastx_trimmer -f 1 -l 5 -i " + tempOutName + " -o " + randomBarcodeOut
        runUNIXCommands(extractRandomBarcodeCommand)
        
# Step 5, select proper hybrid candidate and mappable non-hybrid reads.
    faFileLists = listdir_fullpath_extension(outputdir + "/T", ".fa")
    
    for fa in faFileLists:
        selectReadsCommand = pythondir + "/python " + hiCLIPdir + "/inst/Python/selectHybrid-3.1.0.py -i " + fa + " -o " + outputdir + " -d " + hiCLIPdir
        runUNIXCommands(selectReadsCommand)
        
# Step 6-Phase 1-1, map hybrid reads to rRNAs and tRNAs
    faFileLists = listdir_fullpath_extension(outputdir + "/hybrid/Fasta", ".fa")
    
    for fa in faFileLists:
        outP1MappedName = outputdir + "/hybrid/sam_P1_mapped/P1_" + fa.split("/")[-1].replace(".fa", "").replace("_T", "") + ".sam"
        outP1UnmappedName = outputdir + "/hybrid/fasta_P1_unmapped/" + fa.split("/")[-1].replace("_T", "")
        bowtieCommand = hiCLIPdir + "/inst/bin/bowtie-0.12.7/bowtie -f -p 8 -v 2 -k 1 --best --sam --un " + outP1UnmappedName + " " + hiCLIPdir + "/data/processed/bowtie_index/hs_rRNAs_and_tRNAs/hs_rRNAs_and_tRNAs " + fa + " > " + outP1MappedName
        runUNIXCommands(bowtieCommand)

# Step 6-Phase 1-2 map hybrid reads to rRNAs and tRNAs, only allow unique hit
    faFileLists = listdir_fullpath_extension(outputdir + "/hybrid/Fasta", ".fa")

    for fa in faFileLists:
        outP1UniqueMappedName = outputdir + "/hybrid/sam_P1_mapped_unique/" + fa.split("/")[-1].replace(".fa", "").replace("_T", "") + ".sam"
        bowtieCommand = hiCLIPdir + "/inst/bin/bowtie-0.12.7/bowtie -f -p 8 -m 1 -v 2 --best --sam " + hiCLIPdir + "/data/processed/bowtie_index/hs_rRNAs_and_tRNAs/hs_rRNAs_and_tRNAs " + fa + " > " + outP1UniqueMappedName
        runUNIXCommands(bowtieCommand)
        
# Step 6-Phase 2, map hybrid reads to mtDNA and pre-rRNAs
    faFileLists = listdir_fullpath_extension(outputdir + "/hybrid/fasta_P1_unmapped", ".fa")
    
    for fa in faFileLists:
        outP2UnmappedName = outputdir + "/hybrid/fasta_P2_unmapped/" + fa.split("/")[-1]
        bowtieCommand = hiCLIPdir + "/inst/bin/bowtie-0.12.7/bowtie -f -p 8 -v 2 -k 1 --best --sam --un " + outP2UnmappedName + " " + hiCLIPdir + "/data/processed/bowtie_index/chrM_pre_rRNA/chrM_pre_rRNA " + fa + " > /dev/null"
        runUNIXCommands(bowtieCommand)  

# Step 6-Phase 3, map hybrid reads to transcripts with bowtie
    faFileLists = listdir_fullpath_extension(outputdir + "/hybrid/fasta_P2_unmapped", ".fa")
    
    for fa in faFileLists:
        outP3MappedName = outputdir + "/hybrid/sam_P3_mapped/P3_" + fa.split("/")[-1].replace(".fa", "") + ".sam"
        outP3UnmappedName = outputdir + "/hybrid/fasta_P3_unmapped/" + fa.split("/")[-1]
        bowtieCommand = hiCLIPdir + "/inst/bin/bowtie-0.12.7/bowtie -f -p 8 -v 2 -m 1 --best --strata --sam --un " + outP3UnmappedName + " " + hiCLIPdir + "/data/processed/bowtie_index/longest_mRNA_and_ncRNA/longest_mRNA_and_ncRNA " + fa + " > " + outP3MappedName 
        runUNIXCommands(bowtieCommand)  
        
# Step 6-Phase 4, map hybrid reads to genome with bowtie
    faFileLists = listdir_fullpath_extension(outputdir + "/hybrid/fasta_P3_unmapped", ".fa")
    
    for fa in faFileLists:
        outP4MappedName = outputdir + "/hybrid/sam_P4_mapped/P4_" + fa.split("/")[-1].replace(".fa", "") + ".sam"
        bowtieCommand = hiCLIPdir + "/inst/bin/bowtie-0.12.7/bowtie -f -p 8 -v 2 -m 1 --best --strata --sam " + hiCLIPdir + "/data/processed/bowtie_index/hg19/hg19 " + fa + " > " + outP4MappedName 
        runUNIXCommands(bowtieCommand) 
        
# Step 7. merge three of sam files 
    faFileLists = listdir_fullpath_extension(outputdir + "/hybrid/fasta_P3_unmapped", ".fa")
    
    for fa in faFileLists:
        outP1MappedName = outputdir + "/hybrid/sam_P1_mapped/P1_" + fa.split("/")[-1].replace(".fa", "") + ".sam"
        outP3MappedName = outputdir + "/hybrid/sam_P3_mapped/P3_" + fa.split("/")[-1].replace(".fa", "") + ".sam"
        outP4MappedName = outputdir + "/hybrid/sam_P4_mapped/P4_" + fa.split("/")[-1].replace(".fa", "") + ".sam"
        outMergedSamName = outputdir + "/hybrid/sam/" + fa.split("/")[-1].replace(".fa", "") + ".sam"
        mergeSamCommand = Rdir + "/bin/Rscript " + hiCLIPdir + "/R/mergeSamHybrid.R " + outP1MappedName + " " + outP3MappedName + " " + outP4MappedName + " " +  outMergedSamName
        runUNIXCommands(mergeSamCommand) 
        
# Step 8-1, create bed like file for hybrid reads...
    samFileLists = listdir_fullpath_extension(outputdir + "/hybrid/sam", '.sam')
    randomBarcodeLists = listdir_fullpath_extension(outputdir + "/randomBarcode", '.fa')
    
    for samFile in samFileLists:
        bedOutDir = outputdir + "/hybrid/bed/"
        for randomBarcodes in randomBarcodeLists:
            if samFile.split("/")[-1].replace("hybrid_", "").replace(".sam", "") == randomBarcodes.split("/")[-1].replace("_randomBarcode.fa", ""):
                matchedRandomBarcodes = randomBarcodes
            else:
                continue

        bedCommand = pythondir + '/python ' + hiCLIPdir + '/inst/Python/hiClipBarcode-LeftRight-READ.py ' + samFile +  ' ' + matchedRandomBarcodes  + ' ' + bedOutDir
        runUNIXCommands(bedCommand)

# Step 8-2, create bed like file for hybrid reads uniquly mapped to rRNAs and tRNAs...
    samFileLists = listdir_fullpath_extension(outputdir + "/hybrid/sam_P1_mapped_unique", '.sam')
    randomBarcodeLists = listdir_fullpath_extension(outputdir + "/randomBarcode", '.fa')

    for samFile in samFileLists:
        bedOutDir = outputdir + "/hybrid/rRNA_unique_bed/"
        for randomBarcodes in randomBarcodeLists:
            if samFile.split("/")[-1].replace("hybrid_", "").replace(".sam", "") == randomBarcodes.split("/")[-1].replace("_randomBarcode.fa", ""):
                matchedRandomBarcodes = randomBarcodes
            else:
                continue

        bedCommand = pythondir + '/python ' + hiCLIPdir + '/inst/Python/hiClipBarcode-LeftRight-READ.py ' + samFile +  ' ' + matchedRandomBarcodes  + ' ' + bedOutDir
        runUNIXCommands(bedCommand)


# Step 9-Phase 1, map nonhybrid reads to rRNAs and tRNAs
    faFileLists = listdir_fullpath_extension(outputdir + "/nonHybrid/Fasta", ".fa")
    
    for fa in faFileLists:
        outP1MappedName = outputdir + "/nonHybrid/sam_P1_mapped/P1_" + fa.split("/")[-1].replace(".fa", "").replace("_T", "") + ".sam"
        outP1UnmappedName = outputdir + "/nonHybrid/fasta_P1_unmapped/" + fa.split("/")[-1].replace("_T", "")
        bowtieCommand = hiCLIPdir + "/inst/bin/bowtie-0.12.7/bowtie -f -p 8 -v 2 -k 1 --best --sam --un " + outP1UnmappedName + " " + hiCLIPdir + "/data/processed/bowtie_index/hs_rRNAs_and_tRNAs/hs_rRNAs_and_tRNAs " + fa + " > " + outP1MappedName
        runUNIXCommands(bowtieCommand)  
        
# Step 9-Phase 2, map non-hybrid reads to mtDNA and pre-rRNAs
    faFileLists = listdir_fullpath_extension(outputdir + "/nonHybrid/fasta_P1_unmapped", ".fa")
    
    for fa in faFileLists:
        outP2UnmappedName = outputdir + "/nonHybrid/fasta_P2_unmapped/" + fa.split("/")[-1]
        bowtieCommand = hiCLIPdir + "/inst/bin/bowtie-0.12.7/bowtie -f -p 8 -v 2 -k 1 --best --sam --un " + outP2UnmappedName + " " + hiCLIPdir + "/data/processed/bowtie_index/chrM_pre_rRNA/chrM_pre_rRNA " + fa + " > /dev/null"
        runUNIXCommands(bowtieCommand)  

# Step 9-Phase 3, map non-hybrid reads to transcripts with bowtie
    faFileLists = listdir_fullpath_extension(outputdir + "/nonHybrid/fasta_P2_unmapped", ".fa")
    
    for fa in faFileLists:
        outP3MappedName = outputdir + "/nonHybrid/sam_P3_mapped/P3_" + fa.split("/")[-1].replace(".fa", "") + ".sam"
        outP3UnmappedName = outputdir + "/nonHybrid/fasta_P3_unmapped/" + fa.split("/")[-1]
        bowtieCommand = hiCLIPdir + "/inst/bin/bowtie-0.12.7/bowtie -f -p 8 -v 2 -m 1 --best --strata --sam --un " + outP3UnmappedName + " " + hiCLIPdir + "/data/processed/bowtie_index/longest_mRNA_and_ncRNA/longest_mRNA_and_ncRNA " + fa + " > " + outP3MappedName 
        runUNIXCommands(bowtieCommand)  
        
# Step 9-Phase 4, map non-hybrid reads to genome with bowtie
    faFileLists = listdir_fullpath_extension(outputdir + "/nonHybrid/fasta_P3_unmapped", ".fa")
    
    for fa in faFileLists:
        outP4MappedName = outputdir + "/nonHybrid/sam_P4_mapped/P4_" + fa.split("/")[-1].replace(".fa", "") + ".sam"
        bowtieCommand = hiCLIPdir + "/inst/bin/bowtie-0.12.7/bowtie -f -p 8 -v 2 -m 1 --best --strata --sam " + hiCLIPdir + "/data/processed/bowtie_index/hg19/hg19 " + fa + " > " + outP4MappedName 
        runUNIXCommands(bowtieCommand) 
        
# Step 10. merge three of non-hybrid reads sam files 
    faFileLists = listdir_fullpath_extension(outputdir + "/nonHybrid/fasta_P3_unmapped", ".fa")
    
    for fa in faFileLists:
        outP1MappedName = outputdir + "/nonHybrid/sam_P1_mapped/P1_" + fa.split("/")[-1].replace(".fa", "") + ".sam"
        outP3MappedName = outputdir + "/nonHybrid/sam_P3_mapped/P3_" + fa.split("/")[-1].replace(".fa", "") + ".sam"
        outP4MappedName = outputdir + "/nonHybrid/sam_P4_mapped/P4_" + fa.split("/")[-1].replace(".fa", "") + ".sam"
        outMergedSamName = outputdir + "/nonHybrid/sam/" + fa.split("/")[-1].replace(".fa", "") + ".sam"
        mergeSamCommand = Rdir + "/bin/Rscript " + hiCLIPdir + "/R/mergeSamNonHybrid.R " + outP1MappedName + " " + outP3MappedName + " " + outP4MappedName + " " +  outMergedSamName
        runUNIXCommands(mergeSamCommand) 

        
# Step 11, create bed like file for non hybrid full length reads
    samFileLists = listdir_fullpath_extension(outputdir + "/nonHybrid/sam", '.sam')
    randomBarcodeLists = listdir_fullpath_extension(outputdir + "/randomBarcode", '.fa')
    
    for samFile in samFileLists:
        bedLikeOutDir = outputdir + "/nonHybrid/bedLike/"
        bedOutDir = outputdir + "/nonHybrid/bed/"
        for randomBarcodes in randomBarcodeLists:
            if samFile.split("/")[-1].replace("nonHybrid_", "").replace(".sam", "") == randomBarcodes.split("/")[-1].replace("_randomBarcode.fa", ""):
                matchedRandomBarcodes = randomBarcodes
            else:
                continue

        bedLikeCommand = pythondir + '/python ' + hiCLIPdir + '/inst/Python/hiClipBarcode-one-READ.py ' + samFile +  ' ' + matchedRandomBarcodes  + ' ' + bedLikeOutDir
        bedCommand = pythondir + '/python ' + hiCLIPdir + '/inst/Python/samToBed-2.0.py ' + samFile +  ' ' + matchedRandomBarcodes  + ' ' + bedOutDir + ' 0'
        
        runUNIXCommands(bedLikeCommand)
        runUNIXCommands(bedCommand)          
        
# Step 12-Phase 1, trim nonhybrid reads to 26 nts length
    faFileLists = listdir_fullpath_extension(outputdir + "/nonHybrid/fasta_P2_unmapped", ".fa")
    
    for fa in faFileLists:
        outName = outputdir + "/nonHybrid/RPlike/Fasta_26/" + fa.split("/")[-1].replace(".fa", ".26.fa")
        trim35Command = hiCLIPdir + "/inst/bin/fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64/fastx_trimmer -f 10 -l 35 -Q 33 -i " + fa + " -o " + outName
        runUNIXCommands(trim35Command)
        
        
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
    
    LevenshteinPath = hiCLIPdir + "/inst/bin/python_Levenshtein-0.10.2-py2.7-linux-x86_64.egg"
    sys.path.append(LevenshteinPath)

    runBowtieForRibosomeProfiling(inputfile, outputdir, barcodedir, hiCLIPdir, pythondir, Rdir)
