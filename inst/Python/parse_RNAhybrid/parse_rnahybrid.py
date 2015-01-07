# parse_rnahybrid.py
# qsub -b y -j y ~/Python-2.7.1/python ~/test/parse_rnahybrid.py -i /netscr/yoichiro/STAU1_hiCLIP/data/processed/fasta/a_concatenated_hybrid.fa -o /netscr/yoichiro/test4 -d /netscr/yoichiro/STAU1_hiCLIP -p hybrid
# qsub -b y -j y ~/Python-2.7.1/python ~/test/parse_rnahybrid.py -i /netscr/yoichiro/STAU1_hiCLIP/data/processed/fasta/a_rp.concatenated_hybrid.fa -o /netscr/yoichiro/test4 -d /netscr/yoichiro/STAU1_hiCLIP -p control


import re, tempfile, shutil, os, errno, subprocess, sys

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
        
        
def pyrimidine_enrichment(seq_a, seq_b):
    Y_a = len(re.findall("[CT]", seq_a))
    Y_b = len(re.findall("[CT]", seq_b))
    float(Y_b) / float(Y_a)
    if float(Y_b) > float(Y_a):
        Y_ratio = float(Y_a) / float(Y_b)
    else:
        Y_ratio = float(Y_b) / float(Y_a)
        
    return Y_ratio

def C_enrichment(seq_a, seq_b, trim_l):
    seq_a = seq_a[:(54 - 2 * trim_l)]
    seq_b = seq_b[:(54 - 2 * trim_l)]
    C_a = len(re.findall("[C]", seq_a))
    G_a = len(re.findall("[G]", seq_a))
    C_b = len(re.findall("[C]", seq_b))
    G_b = len(re.findall("[G]", seq_b))
    
    a_score = C_a - G_a
    b_score = C_b - G_b
    
    if a_score > b_score:
        seq_list = [seq_a, seq_b]
    else:
        seq_list = [seq_b, seq_a]
    return seq_list


def stem_merger(stem_list, m_buldge_list):
    stem_max = max(stem_list)
    for i in range(0, len(stem_list) - 1):
        if m_buldge_list[i] != 1000:
            temp_stem_max = stem_list[i] + stem_list[i + 1]
            stem_max = max(stem_max, temp_stem_max)
            
    return stem_max

def parse_stem_buldge_list(upper_buldge, stem, lower_buldge):
    i = 0
    stem_list = []
    buldge_list = []
    
    max_len = max(len(upper_buldge), len(stem), len(lower_buldge))
    
    while i < max_len - 1:
        s_num = 0
        u_num = 0
        l_num = 0
        
        while (stem[i] != " ") & (i < max_len - 1):
            s_num = s_num + 1
            i = i + 1
            
        while (stem[i] == " ") & (i < max_len - 1):
            if upper_buldge[i] != " ":
                u_num = u_num + 1
            
            if lower_buldge[i] != " ":
                l_num = l_num + 1
                
            i = i + 1
        
        stem_list.append(s_num)
        buldge_list.append([u_num, l_num])
    
    return stem_list, buldge_list

def parse_buldge_duplex(stem_list, buldge_list, b_i):
    # b_i is the size of buldge
    parsed_buldge_list = [sum(x) if ((x[0] == 0 and x[1] <= b_i) or (x[1] == 0 and x[0] <= b_i)) else 1000 for x in buldge_list]
    longest_stem = stem_merger(stem_list, parsed_buldge_list)
    return longest_stem

def parse_symmetric_internal_loop_duplex(stem_list, buldge_list, si_i):
    # si_i is the size of loop
    parsed_si_list = [sum(x) if ((x[0] == x[1]) and (1 <= x[0] <= si_i) and (1 <= x[1] <= si_i)) else 1000 for x in buldge_list]
    longest_stem = stem_merger(stem_list, parsed_si_list)
    return longest_stem

def parse_Asymmetric_internal_loop_duplex(stem_list, buldge_list, ai_i):
    # ai_i is the size of loop
    parsed_ai_list = [sum(x) if ((x[0] != x[1]) and (1 <= x[0] <= ai_i) and (1 <= x[1] <= ai_i)) else 1000 for x in buldge_list]
    longest_stem = stem_merger(stem_list, parsed_ai_list)
    return longest_stem

def parse_rnahybrid_line(a,len_seq):
    upper_buldge = a.split(":")[7]
    stem = a.split(":")[8]
    lower_buldge = a.split(":")[10]
    folding_e = a.split(":")[4]
    
    stem_list, buldge_list = parse_stem_buldge_list(upper_buldge, stem, lower_buldge)
    longest_stem = max(stem_list)
    paired_density = float(sum(stem_list)) / float(len_seq)
    
    # The RNA structure type can be classified as buldge (b) symmetric internal loop (si) and asymmetirc internal loop (ai)
    # See The double-stranded-RNA-binding motif: interference and much more (Nature Reviews Molecular Cell Biology 5, 1013-1023 (December 2004)) for detail
    # Complete duplexes: longest_stem
    # Buldge: bi_longest_duplex 
    # Symmetric internal loop : si_longest_duplex where i = j & i >= 1 & j >= 1
    # Asymmetric internal loop: ai_longest_duplex where i > j & i >= 1 & j >= 1
    
    b1_longest_duplex = parse_buldge_duplex(stem_list, buldge_list, b_i = 1)
    b2_longest_duplex = parse_buldge_duplex(stem_list, buldge_list, b_i = 2)
    b3_longest_duplex = parse_buldge_duplex(stem_list, buldge_list, b_i = 3)
    b4_longest_duplex = parse_buldge_duplex(stem_list, buldge_list, b_i = 4)
    b5_longest_duplex = parse_buldge_duplex(stem_list, buldge_list, b_i = 5)
    b6_longest_duplex = parse_buldge_duplex(stem_list, buldge_list, b_i = 6)
    
    si1_longest_duplex = parse_symmetric_internal_loop_duplex(stem_list, buldge_list, si_i = 1)
    si2_longest_duplex = parse_symmetric_internal_loop_duplex(stem_list, buldge_list, si_i = 2)
    si3_longest_duplex = parse_symmetric_internal_loop_duplex(stem_list, buldge_list, si_i = 3)
    si4_longest_duplex = parse_symmetric_internal_loop_duplex(stem_list, buldge_list, si_i = 4)
    si5_longest_duplex = parse_symmetric_internal_loop_duplex(stem_list, buldge_list, si_i = 5)
    si6_longest_duplex = parse_symmetric_internal_loop_duplex(stem_list, buldge_list, si_i = 6)
    
    ai1_longest_duplex = parse_Asymmetric_internal_loop_duplex(stem_list, buldge_list, ai_i = 1)
    ai2_longest_duplex = parse_Asymmetric_internal_loop_duplex(stem_list, buldge_list, ai_i = 2)
    ai3_longest_duplex = parse_Asymmetric_internal_loop_duplex(stem_list, buldge_list, ai_i = 3)
    ai4_longest_duplex = parse_Asymmetric_internal_loop_duplex(stem_list, buldge_list, ai_i = 4)
    ai5_longest_duplex = parse_Asymmetric_internal_loop_duplex(stem_list, buldge_list, ai_i = 5)
    ai6_longest_duplex = parse_Asymmetric_internal_loop_duplex(stem_list, buldge_list, ai_i = 6)
    
    return [folding_e, longest_stem, b1_longest_duplex, b2_longest_duplex, b3_longest_duplex, b4_longest_duplex, b5_longest_duplex, b6_longest_duplex, si1_longest_duplex, si2_longest_duplex, si3_longest_duplex, si4_longest_duplex, si5_longest_duplex, si6_longest_duplex, ai1_longest_duplex, ai2_longest_duplex, ai3_longest_duplex, ai4_longest_duplex, ai5_longest_duplex, ai6_longest_duplex, paired_density]

def run_rnahybrid_analysis(inputfile, outputdir, hiCLIPdir, trim_l = 0, prefix = ""):
    temp_out = outputdir + "/temp" + prefix
    temp_file = temp_out + "/temp.txt"
    if not os.path.exists(temp_out):
        os.makedirs(temp_out)
    
    outputfile = outputdir + "/" + prefix + "max_length_stem" + str(trim_l)  + ".txt"
    
    Crich_fafile = outputdir + "/" + prefix + "C_rich" + str(trim_l)  + ".fa"
    Grich_fafile = outputdir + "/" + prefix + "G_rich" + str(trim_l)  + ".fa"
    
    inFile = open(inputfile,'r')
    outFile = open(outputfile, "w")
    outFile.write("rname\tfolding_e\tlongest_stem\tb1_longest_duplex\tb2_longest_duplex\tb3_longest_duplex\tb4_longest_duplex\tb5_longest_duplex\tb6_longest_duplex\tsi1_longest_duplex\tsi2_longest_duplex\tsi3_longest_duplex\tsi4_longest_duplex\tsi5_longest_duplex\tsi6_longest_duplex\tai1_longest_duplex\tai2_longest_duplex\tai3_longest_duplex\tai4_longest_duplex\tai5_longest_duplex\tai6_longest_duplex\tpaired_density\tpyrimidine_balance\tseq_left\tseq_right\n")
    
    
    CoutFile = open(Crich_fafile, "w")
    GoutFile = open(Grich_fafile, "w")
    
    lineCount = 0
    
    for line in inFile:
        lineCount = lineCount + 1
        
        if lineCount % 2 == 1:
            rname = line.strip()
            rname = rname.replace(">", "")
        else:
            seq = line.strip()
            seq_list = seq.split("&")
            
            if trim_l != 0:
                seq_list[0] = seq_list[0][trim_l:(-trim_l)]
                seq_list[1] = seq_list[1][trim_l:(-trim_l)]
            else:
                1 + 1
            
            Cseq, Gseq = C_enrichment(seq_list[0], seq_list[1], trim_l)
            if len(Cseq) == len(Gseq):
                CoutFile.write(">" + rname + "\n" + Cseq + "\n")
                GoutFile.write(">" + rname + "\n" + Gseq + "\n")
            else:
                None
            
            # pyrimidine_balance = pyrimidine_enrichment(seq_list[0], seq_list[1])
            pyrimidine_balance = 0
            len_seq = len(seq_list[0])
            
            RNAhybrid_cmd = "/lmb/home/yoichiro/RNAprogam/RNAhybrid/bin/RNAhybrid -n 100 -s 3utr_human -c " + seq_list[0] + " " + seq_list[1] + " > " + temp_file
            runUNIXCommands(RNAhybrid_cmd)
            
            with open(temp_file, 'r') as content_file:
                content = content_file.read()
                content = content.replace("\n", "")
                print content
                rnahybrid_out = parse_rnahybrid_line(content, len_seq)
                rnahybrid_out = map(str, rnahybrid_out)
                t_rnahybrid_out = "\t".join(rnahybrid_out)
                outFile.write(rname + "\t" + t_rnahybrid_out + "\t" + str(pyrimidine_balance) + "\t" + seq_list[0] + "\t" + seq_list[1] + "\n")
                
    print "The number of total reads processed..."
    print str(lineCount/2)
    
    inFile.close()
    outFile.close()
    CoutFile.close()
    GoutFile.close()

if __name__ == "__main__":
    from optparse import OptionParser

    parser = OptionParser()

    parser.add_option("-i", "--input", dest="inputDir",
                  help="Input directory which contains fastq files, the fastq files shoule be anarchived and the extension should be .fq", metavar="INPUT")
    parser.add_option("-o", "--output",
                  dest="outputDirectry",
                  help="Directry for output, have to spply")
    parser.add_option("-d", "--hiCLIP_project_dir", dest="hiCLIPDir",
                  help="Path to hiCLIP_project")
    parser.add_option("-p", "--prefix", dest="Prefix",
                  help="Prefix for output file")

    (options, args) = parser.parse_args()

    inputfile = options.inputDir
    outputdir = options.outputDirectry
    hiCLIPdir = options.hiCLIPDir
    prefix = options.Prefix
    
    run_rnahybrid_analysis(inputfile, outputdir, hiCLIPdir, trim_l = 0, prefix = prefix)
    run_rnahybrid_analysis(inputfile, outputdir, hiCLIPdir, trim_l = 5, prefix = prefix)
    run_rnahybrid_analysis(inputfile, outputdir, hiCLIPdir, trim_l = 10, prefix = prefix)
    run_rnahybrid_analysis(inputfile, outputdir, hiCLIPdir, trim_l = 15, prefix = prefix)


