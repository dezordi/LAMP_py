#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
Build in Python 3.6.5+

Dependencies:
argparse 1.1
Bio (biopython) 1.76
BLAST 2.6.0+

External tools:
GLAPD - Genome based LAMP primer designer.

Author:
Filipe Dezordi
zimmer.filipe@gmail.com
https://github.com/dezordi
"""

import argparse, os, re, subprocess, shlex
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast.Applications import NcbimakeblastdbCommandline

parser = argparse.ArgumentParser(description = 'This script receives a list of genomes and for each genome a set of unique k-length sequences are identified. In a final step the GLAPD tool generates n. primer sets for LAMP.')

parser.add_argument("-in", "--input", help="List with genome file names",  required=True)
parser.add_argument("-ln", "--length", help="length of queries",type=int, default=300)
parser.add_argument("-st", "--step", help="length of step for k-mer creation",type=int, default=150)
parser.add_argument("-p", "--threads", help="Threads for blastn analysis",type=int, default=1)
parser.add_argument("-gl","--glapd",help="Path of GLAPD tool",default = False)
parser.add_argument("-num","--number",help="the expected output number of LAMP primer sets", default = 10)
parser.add_argument("-rm","--removetmp",help="Remove temporary files generated through iterions?", default = True)

args = parser.parse_args()
read_file = args.input
query_length = args.length
threads = args.threads
step_length = args.step
glapd_path = args.glapd
num_lamp = args.number
temp_remove = args.removetmp

def getKmers(seq_name,sequence, size, window):
    """
    This function creates a set of kmers from an original sequence, with
    one ID for each kmer.

    Keyword arguments:
    seq_name - fasta sequence name
    sequence - nucleotide sequence
    size - kmer specific size
    window - length of each step of kmer creation
    """

    return [[seq_name+"_"+str(x),sequence[x:x+size].upper()] for x in range(0, len(sequence) - size, window)]

def makeblastdb(data,type):
    """
    This function creates blast databases.

    Keyword arguments:
    data - fasta file
    type - 'nucl' or 'prot'
    """

    clinedb = NcbimakeblastdbCommandline(dbtype=type, input_file = data)
    stdout, stderr = clinedb()
    return(print(f'Makeblastdb for {data}: DONE'))

def runblastn(query_file,database):
    """
    This function runs blastn.

    Keyword arguments:
    query_file - fasta file used as query
    database - pre-formated blast database
    """
    cline = NcbiblastnCommandline(query = query_file, db = database, out = query_file+'.blastn',
    outfmt = 6, evalue = 0.001, max_target_seqs = 1, num_threads = threads, word_size = 11)
    stdout, stderr = cline()
    return(print(f'BLASTn with {query_file} against {database}: DONE'))


if __name__ == '__main__':
    """
    Main Routine
    This block of code is executed, whenever the script
    is started from the command line.
    """

    with open(read_file,'r') as input_list: # read genome list
        lines = input_list.readlines() # read each line as one line
        databases = [] # list for databases
        print('='*120)
        print('Running Makeblastdb')
        print('='*120)
        for var_file in lines: # read each line as variable 'var_file'
            var_file = var_file.rstrip('\n') # remove the '\n' from file name
            remove_N = "sed -i -e '/^>/! s/[Nn]//g;/^$/d' "+ var_file # command-line for first step of GLAPD
            remove_N = shlex.split(remove_N) # Split the string using shell-like syntax.
            cmd_remove_N = subprocess.Popen(remove_N) # Running the firs step of GLAPD
            cmd_remove_N.wait()
            with open(var_file, "r", newline='') as input_handle, open(var_file+".fmt", "w") as output_handle: # read genome file, and creates an output
                sequences = SeqIO.parse(input_handle, "fasta") # creates a variable sequences with fasta sequences from input file
                count = SeqIO.write(sequences, output_handle, "fasta-2line") # creates a 2line formated fasta file for the kmer creation step
                makeblastdb(var_file,'nucl') # makes blast databases with input genome files
                databases.append(var_file) # append the names of databases in databses list


    with open(read_file,'r') as input_list: # read genome list
        queries = [] # list for queries
        lines = input_list.readlines() # read each line as line
        print('='*120)
        print('Creating k-mer files')
        print('='*120)
        for var_file in lines:  # read each line as variable 'var_file'
            var_file = var_file.rstrip('\n')+'.fmt' # variable for 2line formated fasta in place of original fasta
            with open(var_file.replace('.fmt', '')+".kmers","w") as output: #creates an output for kmer sequences
                seqs = [] # creates a list for sequences
                head = 0 # count for selection of header lines
                seq = 1 # count for selection of nucleotide sequence lines
                with open(var_file, 'r') as input_stream: #read 2line formated fasta
                    lines=input_stream.readlines() # read each line as line
                    for i in range(len(lines)): # creates a loop for read all lines
                        while seq and head < len(lines): # quick fix for read all headers and sequences
                            a = re.sub(r'\s.*','',lines[head]) # remove all characters of headers after contig code
                            X = a.rstrip('\n') # remove '\n' from headers
                            Y = getKmers(X,lines[seq], query_length, step_length) #get kmers
                            seq += 2 # count for next header
                            head += 2 # count for next sequence
                            seqs.append(Y) # put header and sequence in list of lists'seqs'
                for i in seqs: # for each list in list of lists 'seqs'
                    for x in i: # for each element in list 'i'
                        output.write(x[0]+'\n') # write header
                        output.write(x[1]+'\n') # write sequence
            print(f'K-mer file for {var_file}: DONE')
            queries.append(var_file.replace('.fmt','')+'.kmers') # put kmer file name in list 'queries'
            if temp_remove == True:
                os.remove(var_file) # remove 2line formated file

    uniq_list = []
    print('='*120)
    print('Running iterative BLASTs')
    print('='*120)
    y = '' #variable created to structure the names of temp files
    for i in queries: # for each kmer file
        if y == '': # in a first momment, Y receives the same name of I
            y = i
        var_range = 1 # used for iteration name
        for x in databases: # for each database
            if (y.find(x) != -1): # avoid running a blastn with a genoma and it's kmer file
                continue
            else:
                if os.path.getsize(i) ==0: # this line avoid blasts with empty query files.
                    continue
                else:
                    runblastn(i,x) # running blastn with i as query (kmer file) and x as database (genome file)
                    with open(y,'r') as input_list, open(y+'.blastn','r') as blast_out: #creates a read file for blast output
                        lines = input_list.readlines() # read lines of kmer file
                        all_blast_lines = blast_out.read() # read lines of blast output
                        no_match = [] # creates a list of no-matched sequences
                        for line in lines: # read each line of query file
                            if ">" in line: # for each header in query file
                                a = line.replace('>','').rstrip('\n') # remove '>' and '\n' for comparrisons against blast output file
                                if (all_blast_lines.find(a)) != -1: # if header of query file exist in blast output, nothing occours
                                    pass
                                else: # else, append header name in 'no_match' list
                                    no_match.append(a)
                        if temp_remove == True:
                            os.remove(i+'.blastn')
                        input_list.close() # close query file
                        blast_out.close() # close blast output
            if var_range == 1: # if is this the first iteration of blast, maintain, just add the number 1 on file name
                with open(i) as seqs, open(y+str(var_range), 'w') as result: # open query file and an output file named 'result'
                    record_dict = SeqIO.to_dict(SeqIO.parse(seqs, 'fasta')) # creates a dictionary with sequences of query file
                    result_records = [record_dict[id_] for id_ in no_match] # search for no_match sequences in query file
                    SeqIO.write(result_records, result, "fasta") # write an output file with queries that don't present match in previous analysis
                    if temp_remove == True:
                        os.remove(i) # remove query file
                    y = i # consider as query for next blastn analysis the file with sequences that don't present match
                    y = y+str(var_range)
                    i = y
                    var_range += 1
            elif var_range <= 10: # for the others iterations, the previos number of iteration is removed, then the number of another iteration is added
                with open(i) as seqs, open(y[:-1]+str(var_range), 'w') as result: # open query file and an output file named 'result'
                    record_dict = SeqIO.to_dict(SeqIO.parse(seqs, 'fasta')) # creates a dictionary with sequences of query file
                    result_records = [record_dict[id_] for id_ in no_match] # search for no_match sequences in query file
                    SeqIO.write(result_records, result, "fasta") # write an output file with queries that don't present match in previous analysis
                    if temp_remove == True:
                        os.remove(i) # remove query file
                    y = y[:-1]
                    y = y+str(var_range)
                    i = y
                    var_range += 1
            elif 10 < var_range <= 100:
                with open(i) as seqs, open(y[:-2]+str(var_range), 'w') as result: # open query file and an output file named 'result'
                    record_dict = SeqIO.to_dict(SeqIO.parse(seqs, 'fasta')) # creates a dictionary with sequences of query file
                    result_records = [record_dict[id_] for id_ in no_match] # search for no_match sequences in query file
                    SeqIO.write(result_records, result, "fasta") # write an output file with queries that don't present match in previous analysis
                    if temp_remove == True:
                        os.remove(i) # remove query file
                    y = y[:-2]
                    y = y+str(var_range)
                    i = y
                    var_range += 1
            elif 100 < var_range <= 1000:
                with open(i) as seqs, open(y[:-3]+str(var_range), 'w') as result: # open query file and an output file named 'result'
                    record_dict = SeqIO.to_dict(SeqIO.parse(seqs, 'fasta')) # creates a dictionary with sequences of query file
                    result_records = [record_dict[id_] for id_ in no_match] # search for no_match sequences in query file
                    SeqIO.write(result_records, result, "fasta") # write an output file with queries that don't present match in previous analysis
                    if temp_remove == True:
                        os.remove(i) # remove query file
                    y = y[:-3]
                    y = y+str(var_range)
                    i = y
                    var_range += 1
            elif 1000 < var_range <= 10000:
                with open(i) as seqs, open(y[:-4]+str(var_range), 'w') as result: # open query file and an output file named 'result'
                    record_dict = SeqIO.to_dict(SeqIO.parse(seqs, 'fasta')) # creates a dictionary with sequences of query file
                    result_records = [record_dict[id_] for id_ in no_match] # search for no_match sequences in query file
                    SeqIO.write(result_records, result, "fasta") # write an output file with queries that don't present match in previous analysis
                    if temp_remove == True:
                        os.remove(i) # remove query file
                    y = y[:-4]
                    y = y+str(var_range)
                    i = y
                    var_range += 1
            elif 10000 < var_range <= 100000:
                with open(i) as seqs, open(y[:-5]+str(var_range), 'w') as result: # open query file and an output file named 'result'
                    record_dict = SeqIO.to_dict(SeqIO.parse(seqs, 'fasta')) # creates a dictionary with sequences of query file
                    result_records = [record_dict[id_] for id_ in no_match] # search for no_match sequences in query file
                    SeqIO.write(result_records, result, "fasta") # write an output file with queries that don't present match in previous analysis
                    if temp_remove == True:
                        os.remove(i) # remove query file
                    y = y[:-5]
                    y = y+str(var_range)
                    i = y
                    var_range += 1
        b = re.sub(r'.kmers.*','.uniq',i) # write output with uniq sequences
        os.rename(i,b)
        if os.path.getsize(b) != 0: # this line avoid to make glapd analysis with empty files...
            uniq_list.append(b)
        y = ''
    print('='*120)
    print('Finished.\nResults with unique regions for each genome are present in *.uniq files')
    print('='*120)
    if glapd_path != False:
        print('='*120)
        print('Running GLAPD')
        print('='*120)
        for x in uniq_list: # for each file of unique kmers
            print(f'Running GLAPD for {x} file')
            unique_file_path = os.path.abspath(x) # atribute the path of uniq file to a variable
            unique_file_dir = unique_file_path.replace(x,'')
            os.chdir(glapd_path) # change directory to GLAPD path - informed by the user
            GLAPD_SINGLE = './Single -in '+unique_file_path+' -out '+x+'_glapd_single' # command-line for first step of GLAPD
            GLAPD_SINGLE = shlex.split(GLAPD_SINGLE) # Split the string using shell-like syntax.
            GLAPD_LAMP = './LAMP -in '+x+'_glapd_single'+' -ref '+unique_file_path+' -out '+unique_file_path+'_lamp_primersq.txt'+' -num '+str(num_lamp)  # command-line for second step of GLAPD
            GLAPD_LAMP = shlex.split(GLAPD_LAMP) # Split the string using shell-like syntax.
            c = subprocess.Popen(GLAPD_SINGLE) # Running the firs step of GLAPD
            c.wait() # wait until the end of first step of GLAPD
            d = subprocess.Popen(GLAPD_LAMP) # Running the second step of GLAPD
            d.wait() # wait until the end of second step of GLAPD
            print(f'GLADP for {x} file: DONE')
            os.chdir(unique_file_dir)
        print('='*120)
        print('Finished.\nResults with lamp primers each genome are present in *_lamp_primers.txt files')
        print('='*120)
