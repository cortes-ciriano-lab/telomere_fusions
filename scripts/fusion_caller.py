import fileinput
from difflib import SequenceMatcher
import regex
import os.path
from os import path
import sys
import re 
import argparse


def similar(a, b):
    return SequenceMatcher(None, a, b).ratio()


def initialize_parser():
    parser = argparse.ArgumentParser(
        description=(
            'This tool serves to detect telomere fusions in DNA sequencing data sets'
            'Type .py --help for '
            'further instructions.'
        )
    )
    parser.add_argument(
        '--mode',
        help='Action to perform',
        required=True,
        choices = ["callfusions", "extractmates", "summarise"]
    )
    parser.add_argument(
        '--readIDs',
        help='File containing the IDs for the reads with telomere fusions identified after running using mode "callfusions"',
        required=False,
        default=None
    )
    parser.add_argument(
        '--fusionsfile',
        help='File containing the candidate fusions identified when running the mode callfusions',
        required=False,
        default=None
    )
    parser.add_argument(
        '--matesfile',
        help='File generated when running the mode extractmates, containing the mates for the reads containing the candidate fusions identified when running the mode callfusions',
        required=False,
        default=None
    )
    parser.add_argument(
        '--alignmentinfo',
        help='File generated when running coverage_info.sh and containing read alignemnt information',
        required=False,
        default=" "
    )
    parser.add_argument(
        '--outprefix',
        help='Prefix for output files',
        required=False
    )
    parser.add_argument('reads', metavar='FILE', nargs='*', help='files to read, if empty, stdin is used')
    return parser



'''
Function to extract fusions from stdin
'''
def extract_fusions(args):
    read_names = []
    tot=0
    forwardtwice="TTAGGGTTAGGG"
    reversetwice="CCCTAACCCTAA"
    fusions_out = open("{}_fusions".format(args.outprefix),'w+')
    readIDs_out = open("{}_readIDs".format(args.outprefix),'w+')
    for line in fileinput.input(args.reads):
        line = line.split("\t")
        # matching forward
        matches_forward = regex.findall("("+forwardtwice+"){s<=2}", line[9])
        matches_reverse = regex.findall("("+reversetwice+"){s<=2}", line[9])
        if len(matches_forward)>0 and len(matches_reverse) >0:
            read_names.append(line[0])
            fusions_out.write("{}\t{}\t{}\n".format("\t".join(line).strip("\n"),len(matches_forward), len(matches_reverse)))
            readIDs_out.write(line[0])
            readIDs_out.write("\n")
        tot+=1
    fusions_out.write("{}".format(tot))
    fusions_out.close()
    readIDs_out.close()  


'''
Function to extract read mates
'''
def extract_mates(args):
    keep=[]
    if args.readIDs is not None:
        readIDs_file = args.readIDs
    else:
        readIDs_file = "{}_readIDs".format(args.outprefix)
    mates_out = open("{}_mates".format(args.outprefix),'w+')
    readIDs = []
    
    with open(readIDs_file,'r') as file:
        for line in file: 
            line = line.strip() #split("\t") #or some other preprocessing
            readIDs.append(line) #storing everything in memory!
    
    for line in sys.stdin:
        line = line.split("\t")
        if line[0] in readIDs:
            keep.append("\t".join(line))
    
    for i,x in enumerate(keep):
        mates_out.write(x)
    mates_out.close()


'''
Function to match read pairs
'''
def summarise(fusions_file, mates_file, args):
    read_names=[]
    output=[]
    forward="TTAGGG" 
    reverse="CCCTAA"
    forward_twice="TTAGGGTTAGGG"
    reverse_twice="CCCTAACCCTAA"
    forward_one_repeat ="TTAGGG"
    reverse_one_repeat ="CCCTAA"
    
    reads_visited=[]

    summary_out = open("{}_fusions_summary".format(args.outprefix),'w+')
    
    with open(fusions_file) as line_in:
        for line in line_in:
            line = line.split("\t")
            line = line[0:12]
    
            if len(line)>1:
                matches_forward0 = regex.findall("("+forward+")", line[9])
                matches_reverse0 = regex.findall("("+reverse+")", line[9])
                matches_forward_twice0 = regex.findall("("+forward_twice+")", line[9])
                matches_reverse_twice0 = regex.findall("("+reverse_twice+")", line[9])
                matches_forward1 = regex.findall("("+forward+"){s<=1}", line[9])
                matches_reverse1 = regex.findall("("+reverse+"){s<=1}", line[9])
                matches_forward2 = regex.findall("("+forward+"){s<=2}", line[9])
                matches_reverse2 = regex.findall("("+reverse+"){s<=2}", line[9])
               
                # multiple matches
                if (len(matches_forward_twice0) > 0 and len(matches_reverse_twice0) > 0): # or (len(matches_forward0) > 0 and len(matches_reverse1) > 0) or (len(matches_forward1) > 0 and len(matches_reverse0) > 0):
                    if (len(matches_forward0) > 1 and len(matches_reverse0) > 1):
                        match_span_forward = re.search(forward_twice,line[9]).span() 
                        match_span_reverse = re.search(reverse_twice,line[9]).span()
    
                    if type(match_span_forward) is tuple:
                        match_span_forward=[match_span_forward]
                    if type(match_span_reverse) is tuple:
                        match_span_reverse=[match_span_reverse]
                    mate=[""]*15 
                    if path.exists(mates_file):
                        with open(mates_file) as file_in:
                            for linenow in file_in:
                                linenow = linenow.split("\t")   
                                #if line[0] == linenow[0] and line[1] != linenow[2]:
                                if line[0] == linenow[0] and line[9] != linenow[9] and line[7] == linenow[3]:
                                    mate=linenow
                        read_names.append(line[0])
                        if match_span_forward[-1][0] < match_span_reverse[-1][0]:
                            middle_sequence = line[9][ match_span_forward[-1][1] : match_span_reverse[0][0] ]
                            middle_sequence = re.sub(forward_one_repeat, "_forward_", middle_sequence)
                            middle_sequence = re.sub(reverse_one_repeat, "_reverse_", middle_sequence)
                            # left
                            left = line[9][ 0 : match_span_forward[0][0] ]
                            # right
                            right = line[9][ match_span_reverse[-1][1] : len(line[9]) ]
    
                            matches_forward = regex.findall(forward, mate[9])
                            matches_reverse = regex.findall(reverse, mate[9])
                            match_telo="No\t0\t0"
                            if len(matches_forward)>1 or  len(matches_reverse) >1: 
                                match_telo="Yes\t{}\t{}".format(len(matches_forward), len(matches_reverse) )
    
                            if line[9][0:30] not in reads_visited:
                                output.append([line, len(matches_forward0), len(matches_reverse0), len(matches_forward1), len(matches_reverse1), len(matches_forward2), len(matches_reverse2), "end-to-head",left, middle_sequence, right, mate[1], mate[2],mate[3],mate[4],mate[5], mate[9], match_telo])
                                reads_visited.append(line[9][0:30])
                        else:
                            middle_sequence = line[9][ match_span_reverse[-1][1] : match_span_forward[0][0] ]
                            middle_sequence = re.sub(forward_one_repeat, "_forward_", middle_sequence)
                            middle_sequence = re.sub(reverse_one_repeat, "_reverse_", middle_sequence)
                            # left
                            left = line[9][ 0 : match_span_reverse[0][0] ]
                            # right
                            right = line[9][ match_span_forward[-1][1] : len(line[9]) ]
                            
                            matches_forward = regex.findall(forward, mate[9])
                            matches_reverse = regex.findall(reverse, mate[9])
                            match_telo="No\t0\t0"
                            if len(matches_forward)>0 or  len(matches_reverse) >0: 
                                match_telo="Yes\t{}\t{}".format(len(matches_forward), len(matches_reverse) )
                            if line[9][0:30] not in reads_visited:
                                output.append([line, len(matches_forward0), len(matches_reverse0), len(matches_forward1), len(matches_reverse1), len(matches_forward2), len(matches_reverse2), "outward",left, middle_sequence,right, mate[1], mate[2],mate[3],mate[4],mate[5], mate[9], match_telo])
                                reads_visited.append(line[9][0:30])
            else:
                global count
                count=line[0].strip("\n")
    # load read alignemnt information
    with open(args.alignmentinfo,'r') as file:
        for line in file:
            line=line.strip("\n")
            line=line.split(",")
            Path_file=line[0]
            File=line[1]
            Total_reads=line[2]
            Supplementary_reads=line[3]
            Duplicate_reads=line[4]
            Paired_reads=line[5]

    for i,x in enumerate(output):
        summary_out.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format("\t".join(x[0][:-2]).strip("\n"),x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8],x[9],x[10],x[11],x[12],x[13],x[14],x[15],x[16],x[17],count,args.outprefix,Path_file,File,Total_reads,Supplementary_reads,Duplicate_reads,Paired_reads))
    summary_out.close()
    
    
#---------------------------------------------------------------------------------------------------------------
'''
Run..
'''
def main():
    parser = initialize_parser()
    args = parser.parse_args()
    #print(args)

    if args.mode == "callfusions":
        extract_fusions(args)

    if args.mode == "extractmates":
        extract_mates(args)

    if args.mode == "summarise":
        if args.fusionsfile is None:
            fusions_file = "{}_fusions".format(args.outprefix)
        else:
            fusions_file = args.fusionsfile
        if args.matesfile is None:
            mates_file = "{}_mates".format(args.outprefix)
        else:
            mates_file = args.matesfile

        if not  path.exists(args.alignmentinfo):
            raise  FileNotFoundError("File with coverage information does not exist")
        if not path.exists(mates_file):
            raise ValueError("File with read mates does not exist")
        if not path.exists(fusions_file):
            raise ValueError("File with telomere fusions does not exist")

        summarise(fusions_file,mates_file,args) 

if __name__ == '__main__':
    main()


