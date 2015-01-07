import sys, subprocess, os

class Unbuffered:
    def __init__(self, stream):
        self.stream = stream
    def write(self, data):
        self.stream.write(data)
        self.stream.flush()
    def __getattr__(self, attr):
        return getattr(self.stream, attr)

sys.stdout=Unbuffered(sys.stdout)

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
    

def add_color(inputps, constraint, outputfile, wrapperfile, wrapperfile2):
    inFile = open(inputps,'r')
    wrapperFile = open(wrapperfile, "r")
    wrapperFile2 = open(wrapperfile2, "r")
    
    file = open(outputfile, "w")

    ii1 = True
    ii2 = True
    
    line = inFile.readline()
    
    while ii1 == True:
        
        if line.startswith("/pairs") == False:
            bowtieIn = line
            file.write(bowtieIn)
            line = inFile.readline()
        else:
            bowtieIn = line
            file.write(bowtieIn)
            ii1 = False
            line = inFile.readline()
            
            while ii2 == True:
                if line.startswith("init") == False:
                    bowtieIn = line
                    file.write(bowtieIn)
                    line = inFile.readline()
                else:
                    bowtieIn = line
                    file.write(bowtieIn)
                    ii2 = False
                    line = inFile.readline()
    
    for wline in wrapperFile:
        wIn = wline
        file.write(wIn)
        
    for db in constraint:
        if db == ".":
            color_code = "0.5"
        elif db == "1":
            color_code = "0"
        else:
            color_code = "1"
        file.write(color_code + "\n")
        
    for wline2 in wrapperFile2:
        wIn2 = wline2
        file.write(wIn2)

    file.close()

    inFile.close()
    
def add_color_batch(inputfile, psdir, outputdir, wrapperfile, wrapperfile2):
    input_constraint = inputfile 
    inFile = open(input_constraint,'r')
    
    if outputdir[-1] == "/":
        outputdir = outputdir[:-1]
    os.chdir(outputdir)
    
    lineCount = 0
    
    for line in inFile:
        line = line.strip()
        lineCount = lineCount + 1
        if lineCount % 2 == 1:
            gene_id = line.replace("> ", "")
            # print(gene_id)
        else:
            constraint = line
            # print(constraint)
            
            inputps = psdir + "/" + gene_id + "_ss.ps"
            outputfile = outputdir + "/" + gene_id + "_ss_colored.ps"

            if os.path.isfile(inputps):
                add_color(inputps, constraint, outputfile, wrapperfile, wrapperfile2)
            

            

if __name__ == "__main__":
    from optparse import OptionParser

    parser = OptionParser()

    parser.add_option("-i", "--input", dest="inputFile",
                  help="Input constraint file", metavar="INPUT")
    parser.add_option("-p", "--psdirectory",
                  dest="psDirectory",
                  help="Directry for ps files")
    parser.add_option("-o", "--output",
                  dest="outputDirectry",
                  help="Directry for output, have to spply")
    parser.add_option("-w", "--wrappper", dest="wrappperFile",
                  help="wrappperFile")
    parser.add_option("-x", "--wrappper2", dest="wrappperFile2",
                  help="wrappperFile2")

    (options, args) = parser.parse_args()

    inputfile = options.inputFile
    psdir = options.psDirectory
    outputdir = options.outputDirectry
    wrapperfile = options.wrappperFile
    wrapperfile2 = options.wrappperFile2
    #  prefix = options.prefix

    add_color_batch(inputfile, psdir, outputdir, wrapperfile, wrapperfile2)
