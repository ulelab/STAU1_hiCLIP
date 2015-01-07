'''
Created on 20 Sep 2012

@author: Nejc Haberman
'''

'''
input: 
    fin_sam    - path to .SAM file
    fin_random_barcode    - path to random barcode file
    output_PATH   - path for results: total.tab and collapsed.tab
    
barcode_for_hiClip method:
 
The method read .SAM file and write collapsed.TAB and total.TAB file. In both files whe ignore "4" strand. For the output we get IDs from diffrent reads 
for every position. If they have the same position then we have more IDs in a row (id1, id2, id3 for example, then the NUM of reads is 3)

In collpased outout we do the same except for IDs that are on the same position we also check their random barcodes. If they have the same barcodes we put IDs in to ()
For example: (id1, id2, id3), id5, id6  where the NUM_of_reads is 3- these reads are on the same position and first 3 IDs have the same random barcode.
'''


import sys

#method reads random barcodes from fin_name file and returns dictonary of barcodes where the key is read ID
def load_random_barcodes(fin_name):
    fin = open(fin_name, "rt")
    rb = fin.readline()
    rand_barcodes = {}
    while rb:
        id = rb.rstrip('\n').rstrip('\r')
        id = id[1:]
        rb = fin.readline()
        if rb:
            barcode = rb.rstrip('\n').rstrip('\r')
            rand_barcodes.setdefault(id,barcode)
            rb = fin.readline()
    fin.close()
    return rand_barcodes


#from array (with IDs) returns string if IDs
def array_to_str(array):
    string = ""
    for i in range(0, array.__len__() -1):
        string += str(array[i]) + ","
    string += str(array[-1])
    return string



# compares 2 barcodes and ignores N values. example: AATAT is the same as ANTAT 
def compare_barcode (a, b):
    for l, r in map(None, a, b):
        if (l != r):
            if (l != 'N') and (r != 'N'):
                return False
    return True


# from array of IDs and random barcodes, returns string of IDs. IDs with the same barcode are caputred in () and number of uniq IDs. Returns NA,0 for missing barcodes
# method also compares barcodes with missing informacion that includes N
def get_collapsed_ids (ids, random_barcodes):
    collapsed_ids_string = ""
    collapsed_barcodes = {}
    missing_barcodes = {}   #if they include N in barcode. example AANTA
    for i in range(0, ids.__len__()):
        id = ids[i]
        barcode = random_barcodes.get(id)
        if barcode == None:
            print "missing barcode for id: " + str(id)
        else:
            if str(barcode).find('N') > 0:
                missing_barcodes.setdefault(barcode,[]).append(id)
            else:
                collapsed_barcodes.setdefault(barcode, []).append(id)
            if missing_barcodes.__len__() > 0:  #if there are any missing (once that have N) barcodes we add IDs to collapsed one if there are any matching barcodes
                for collapsed_barcode, collapsed_ids in collapsed_barcodes.items():
                    for missing_barcode, missing_ids in missing_barcodes.items():
                        if compare_barcode(collapsed_barcode, missing_barcode):     #if there is a match then we add IDs from missing dicti to matcing one and remove it
                            collapsed_barcodes.setdefault(collapsed_barcode, []).extend(missing_ids)
                            missing_barcodes.pop(missing_barcode)   #remove
    number_of_uniq_ids = 0  
    for barcode, collapsed_ids in collapsed_barcodes.items():   #creating a string of IDs for uniq and collapsed () IDs
        if collapsed_ids.__len__() > 1:
            collapsed_ids_string += "("
            for j in range(0, collapsed_ids.__len__() -1):
                collapsed_ids_string += collapsed_ids[j] + ","
            collapsed_ids_string += collapsed_ids[-1] + "),"
        else:
            collapsed_ids_string += collapsed_ids[-1] +","
        number_of_uniq_ids = collapsed_barcodes.__len__()   #number of ids with uniq barcode
    
    if number_of_uniq_ids > 0:
        return collapsed_ids_string.rstrip(','), number_of_uniq_ids
    else:
        return "NA", 0
        
#returns position which depends on exon/intron junction from read length  
def get_position(position, read_length, strand):
    if strand == '+':
        return position
    if strand == '-':   #read_length ofr example: 6M676N17M1033N3M -> (exon)M(intron)N(exon)M(intron)N(exon)M
        tokens = read_length.split('M')
        exon_intron_NUM = tokens.__len__()
        introns = []
        exons = []
        exons.append(int(tokens[0]))  #first one is always exon
        for i in range(1,exon_intron_NUM -1):   #getting length for all exons and introns
            intron_exon = tokens[i]
            intron, exon = intron_exon.split('N')
            exons.append(int(exon))
            introns.append(int(intron))

        length = sum(exons) + sum(introns)
        position += length -1
        return position
    else:
        print "unknown strand: " + strand
        return None
    
           

def barcode_for_hiClip (sam_fname, barcode_fname, outPATH):
    
    #output files
    tokens = sam_fname.rsplit('/')
    file_name = tokens[-1]
    file_name = file_name.rstrip("sam").rstrip('.')
    fname_out = outPATH + file_name + "-barcode_for_hiClip"
    fout_total = open (fname_out + "-total.tab", 'w')
    fout_collapsed = open(fname_out + "-collapsed.tab", 'w')
    fout_total.write("ids\tchr\t0_based_position\t1_based_position\tstrand\tNUM of reads\n")
    fout_collapsed.write("ids\tchr\t0_based_position\t1_based_position\tstrand\tNUM of reads\n")
    random_barcodes = load_random_barcodes(barcode_fname)
    
    fin_sam = open(sam_fname, "rt")
    line = fin_sam.readline()
    reads = {}
    while line[0] == '@':    #header of .SAM
        line = fin_sam.readline()
        
    while line:
        tokens = line.rstrip('\n').rstrip('\r').split('\t')
        read_id = tokens[0]
        strand = tokens[1]
        chr = tokens[2]
        position = tokens[3]
        read_length = tokens[5]
        if strand != '4':   #we ignore 4 strands
            if strand == '0': strand = '+'
            if strand == '16': strand = '-'
            position = get_position(int(position), read_length, strand)
            reads.setdefault(chr, {}).setdefault(strand,{}).setdefault(str(position),[]).append(read_id)
        line = fin_sam.readline()
        
    for chr, chr_reads in reads.items():
        for strand, strand_reads in chr_reads.items():
            for position, ids in strand_reads.items():
                fout_total.write(array_to_str(ids) + '\t' + chr + '\t' + str(int(position) - 1) + '\t' + position + '\t' + strand + '\t' + str(ids.__len__()) + '\n')
                collapsed_ids_in_string, collapsed_ids_number = get_collapsed_ids(ids, random_barcodes)            
                fout_collapsed.write(collapsed_ids_in_string + '\t' + chr + '\t' + str(int(position) - 1) + '\t' + position + '\t' + strand + '\t' + str(collapsed_ids_number) + '\n')

    fout_total.close()
    fin_sam.close()



if sys.argv.__len__() == 4:
    sam_fname = sys.argv[1]
    barcode_fname = sys.argv[2]
    outPATH = sys.argv[3]
    barcode_for_hiClip (sam_fname, barcode_fname, outPATH)
else:
    print "error:\t3 arguments are needed\n" + '\n' +"example:\t $ python barcode_for_hiClip.py sam_file.sam rand_barcode.fa output_PATH"
    print "\narguments: \nsam_file.sam\t- path to .SAM file\nrand_barcode.fa\t- path to random barcode file\noutput_PATH\t- path for output_files: total.tab and collapsed.tab\ns" 

    
'''
PATH = "/home/nebo/MRC-LMB/MRC-LMB-data/2012.09.20@barcode_for_hiClip-Yoichiro/one-READ/"
barcode_fname = "DOX_LigMinus_randomBarcode.fa"
sam_fname = "nonHybrid_DOX_LigMinus_T.sam"
barcode_for_hiClip (PATH + sam_fname, PATH + barcode_fname, PATH)
'''
