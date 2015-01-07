## Yoichiro Sugimoto 2012/Aug/17
## input shold be fasta file
## swap barcode from random 3 barcode + 4 experiment barcode + 2 random barcode 
#            to 4 experiment barcode + 5 random barcode

import sys

def swapBarcode(inputfile, outputfile):

    inFile = open(inputfile,'r')
    outFile = open(outputfile, "w")

    lineCount = 0
    for line in inFile:
        lineCount = lineCount + 1
        if lineCount % 2 == 1:
            header = line.strip()
            outFile.write(header + "\n")
        elif lineCount % 2 == 0:
            read = line.strip()
            randomBarcode = read[0:3] + read[7:9]
            experimentBarcode = read[3:7]
            seqRead = read[9:]
            swapped = experimentBarcode + randomBarcode + seqRead
            outFile.write(swapped + "\n")

    inFile.close()
    outFile.close()   

if __name__ == "__main__":
    from optparse import OptionParser

    parser = OptionParser()

    parser.add_option("-i", "--input", dest="inputFile",
                  help="Input fasta file name, the experimental barcode shoule be located at [3:6] and random barcode at [0:2] and [7:8] on 0-base", metavar="INPUT")
    parser.add_option("-o", "--output",
                  dest="outputFile",
                  help="Output file name, have to spply")

    (options, args) = parser.parse_args()

    inputfile = options.inputFile
    outputfile = options.outputFile

    swapBarcode(inputfile, outputfile)
