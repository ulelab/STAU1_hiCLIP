#! /usr/bin/env python

'''
Created on 24 Aug 2012

@author: Nejc Haberman
'''        
import sys
import collections

# input: 
#    fin_sam    - path to .SAM file
#    fin_random_barcode    - path to random barcode file
#    output_PATH   - path for results: total.bed and collapsed.bed
#    shift   - integer number for shifting a start position

# sam_file_to_bed method:
# 
# The method read .SAM file and write collapsed.BED/.TAB and total.BED/.TAB file. In both files whe ignore "4" strand. In Collapsed we also ignore duplicate barcodes.
# In .TAB files we also have information about IDs in the first column.


def sam_file_to_bed (fin_sam, fin_random_barcode, output_PATH, shift):
    random_barcodes = load_random_barcodes(fin_random_barcode)
    fin = open(fin_sam, "rt")
    sam = fin.readline()
    reads = {}
    while sam[0] == '@':    #header of .SAM
        sam = fin.readline()
    while sam:
        tokens = sam.rstrip('\n').rstrip('\r').split('\t')
        read_id = tokens[0]
        strand = tokens[1]
        chr = tokens[2]
        position = int(tokens[3])
        read_length = tokens[5]
        
        if strand != '4': #we ignore 4
            if strand == '0':   strand = '+'
            if strand == '16':  strand = '-'
            position = get_position(position, read_length, strand, shift)
            if position == None:
                print "position set problem:\t" + read_id
            else:
                barcode = random_barcodes.get(read_id) #we get random barcode from .fa file
                if barcode == None:
                    print "missing barcode for:\t" + read_id
                else:
                    read = strand, read_id
                    reads.setdefault(chr, {}).setdefault(position, {}).setdefault(strand,[]).append(read_id)
        sam = fin.readline()
    fin.close()
    
    # .BED files
    tokens = fin_sam.rsplit('/')
    file_name = tokens[-1]
    file_name = file_name.rstrip("sam").rstrip('.')
    fname_bed_total = file_name + "_total" + ".bed"
    fname_bed_collapsed = file_name + "_collapsed" + ".bed"
    fname_bed_total_with_ids = file_name + "_IDs_total.tab"
    fname_bed_collapsed_with_ids = file_name + "_IDs_collapsed.tab"
    fout_bed_total = open(output_PATH + fname_bed_total, "w")
    fout_bed_collapsed = open(output_PATH + fname_bed_collapsed, "w")
    fout_tab_total = open(output_PATH + fname_bed_total_with_ids, "w")
    fout_tab_collapsed = open(output_PATH + fname_bed_collapsed_with_ids, "w")
    fout_tab_total.write("IDs\tchr\tposition (0 based)\tposition (1 based)\tcount\n")
    fout_tab_collapsed.write("IDs\tchr\tposition (0 based)\tposition (1 based)\tcount\n")
    fout_bed_total.write(get_bed_header(file_name + "_total"))
    fout_bed_collapsed.write(get_bed_header(file_name + "_collapsed"))
    sum_total = 0
    sum_collapsed = 0
    
    for chr, chr_reads in reads.items():
        sorted_chr_reads = collections.OrderedDict(sorted(chr_reads.items()))   #sort positions in every chromosome
        for position, position_reads in sorted_chr_reads.items():
            for strand, ids in position_reads.items():
                if strand == '-':   strand_char = '-'
                else:   strand_char = ''
                fout_bed_total.write(chr + '\t' + str(position-1) + '\t' + str(position) + '\t' + strand_char + str(ids.__len__()) + '\n')
                fout_tab_total.write(array_to_str(ids) + '\t' + chr + '\t' + str(position-1) + '\t' + str(position) + '\t' + strand_char + str(ids.__len__()) + '\n')
                collapsed_ids_in_string, collapsed_ids_number = get_collapsed_ids(ids, random_barcodes)    
                fout_bed_collapsed.write(chr + '\t' + str(position-1) + '\t' + str(position) + '\t' + strand_char + str(collapsed_ids_number) + '\n')
                fout_tab_collapsed.write(collapsed_ids_in_string + '\t' + chr + '\t' + str(position-1) + '\t' + str(position) + '\t' + strand_char + str(collapsed_ids_number) + '\n')
                sum_total += ids.__len__()
                sum_collapsed += collapsed_ids_number
    fout_bed_total.close()
    fout_bed_collapsed.close()
    fout_tab_total.close()
    fout_tab_collapsed.close()
    print ('\n')
    print (fname_bed_total + "\t successfuly exported")
    print (fname_bed_collapsed + "\t successfuly exported")
    print (fname_bed_total_with_ids + "\t successfuly exported")
    print (fname_bed_collapsed_with_ids + "\t successfuly exported")
    ratio = float(sum_collapsed) / float(sum_total)
    print ("Ratio number (collapsed/total): " + str(ratio) + '\n')


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
    
#returns shifted position
def get_position(position, read_length, strand, shift):
    
    if shift == 0 and strand == '+':
        return position
    #get lenght information for exons and introns
    tokens = read_length.split('M')
    exon_intron_NUM = tokens.__len__()
    introns = []
    exons = []
    exons.append(int(tokens[0]))  #first one is always exon
    for i in range(1,exon_intron_NUM -1):
        intron_exon = tokens[i]
        intron, exon = intron_exon.split('N')
        exons.append(int(exon))
        introns.append(int(intron))
    length = sum(exons) + sum(introns)
    if shift == 0 and strand == '-':
        position += length -1
        return position
    if shift > sum(exons):
        print "shift value is out of range\t" + "position:\t" + str(position) 
        return None
    if strand == '+':
        if shift < exons[0]:
            position += shift
            return position
        else:
            for i in range(1, exons.__len__()):   #sum as many introns as shift overlaps the exons
                if (shift >= sum(exons[0:i])) and (shift < sum(exons[0:i+1])):
                    position += shift + sum(introns[0:i])
                    return position
    elif strand == '-':
        if shift < exons[-1]:
            position += length -1 -shift
            return position
        else:
            len = exons.__len__()
            for i in range(len-1, 0,-1):
                if (shift >= sum(exons[i:len]) and shift < sum(exons[i-1:len])):
                    position += length -1 - sum(introns[i-1:len]) -shift
                    return position 
    else:
        print "unknown strand: " + strand
        return None

    
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
    
    
# compares 2 barcodes and ignores N values. example: AATAT is the same as ANTAT 
def compare_barcode (a, b):
    for l, r in map(None, a, b):
        if (l != r):
            if (l != 'N') and (r != 'N'):
                return False
    return True

#from array (with IDs) returns string if IDs
def array_to_str(array):
    string = ""
    for i in range(0, array.__len__() -1):
        string += str(array[i]) + ","
    string += str(array[-1])
    return string

    
def get_bed_header(fname):
    header = [
        "track type=bedGraph name=" + '"' + fname + '"',
        "description=" + '"' + fname + '"',
        "db=hg19 regions="+'"'+'"',
        "mapped_tos=" + '"' + "hg19" + '"',
        "altColor="+'"'+"200,120,59"+'"',
        "color="+'"'+"120,101,172"+'"',
        "exp_ids="+'"'+'"',
        "maxHeightPixels="+'"'+"100:50:0"+'"',
        "visibility=" + '"' + "full" + '"',
        "priority=" + '"' + "20" +'"',
        "grouping_operation=" + '"' + "sum" + '"',
        "mapped_to=" + '"' + "hg19" + '"',
        "detailed_description=" + '"'+'"',
        "group=" + '"' + '"',
        "annot_vers=" + '"' + "ensembl67" + '"' 
    ]
    return ("%s\n" % " ".join(header))
        
'''
PATH = "/home/nebo/MRC/MRC-LMB-data/2012.08.20@RandomBarcode-Yoichiro/id4-test/"
sam_file_to_bed(PATH + "id4_KD_2_RP_RUTT_accepted_hits.sam",  PATH + "id4_KD_2_RP_randomBarcode.fa", PATH, 0)
'''

if sys.argv.__len__() == 4:
    sam_fname = sys.argv[1]
    fa_fname = sys.argv[2]
    out_PATH = sys.argv[3]
    sam_file_to_bed(sam_fname, fa_fname, out_PATH, 0)
elif sys.argv.__len__() == 5:
    sam_fname = sys.argv[1]
    fa_fname = sys.argv[2]
    out_PATH = sys.argv[3]
    shift_num = sys.argv[4]
    shift = int(shift_num)
    sam_file_to_bed(sam_fname, fa_fname, out_PATH, shift)
else:
    #print str(sys.argv.__len__())
    print "error:\tat least 3 arguments are needed\n" + '\n' +"example:\t $ python samToBed.py sam_file.sam rand_barcode.fa outputPATH shiftNUM"
    print "\narguments: \nsam_file.sam\t- path to .SAM file\nrand_barcode.fa\t- path to random barcode file\noutputPATH\t- path for results: total.bed and collapsed.bed\nshiftNUM\t- integer number for shifting a start position (default is 0)\n" 

    
    
    
    
