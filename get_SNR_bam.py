#!/usr/bin/env python


# Import libraries
import os
import pysam
import argparse
from random import random
import errno
import time



############################################################
#                                                          #
#                         Functions                        #
#                                                          #
############################################################

# check if folder exists
def make_sure_path_exists(path):
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise


# CIGAR parse function
def readLength_CIGAR(cigar):
    """
    Reads the CIGAR string to take it into consideration when defining
    the read length
    """

    digitHold=""
    cigarChars=""
    numList=[]
    for char_c in cigar:
        if char_c.isdigit():
            digitHold+=str(char_c)
        else:
            numList.append(int(digitHold))
            cigarChars+=char_c
            digitHold=""

    result=-1 #subtract one because of first position
    for char in range(0, len(cigarChars)):
        if not cigarChars[char] in ["I","D"]:
            result+=numList[char]

    return result


# Correct strand info on extra fields
def setStrand(line):
    """
    Takes a split line, searches for the strand field and changes it
    """

    for i in range(0, len(line)):
        if line[i].startswith("XS:A:"):
            if line[i][-1]=="+":
                line[i]="XS:A:-"
            else:
                line[i]="XS:A:+"

    return line


# Remove temporary files
def removeTemp(bamid):
    filesInDir=os.listdir(os.getcwd())
    for f in filesInDir:
        if bamid in f:
            os.remove(f)




############################################################
#                                                          #
#                           MAIN                           #
#                                                          #
############################################################

# Input arguments
parser = argparse.ArgumentParser(description="""Returns the 
    coordinates for the last base of mNET-seq reads. 
    It considers the library to be secondstranded, 
    with the last base to be incorporated 
    in read 2 and the strand information in read 1. 
    Result will be filename_sorted.bam file and 
    respective bam index (bai).""")


parser.add_argument('-f', '--filepath', metavar='filepath', nargs='+',
                   help='One or more file paths, separated by spaces.')
parser.add_argument('-s', '--filename', metavar='filename', nargs='+',
                   help="""New file prefixes (no extension), separated by 
                   spaces. In the same number as the file paths.""")
parser.add_argument('-d','--outdir', dest='outDir', nargs='?', 
    default="./", help="Output directory. Defaults to './'.")
args = parser.parse_args()

files = args.filepath
filenames = args.filename
outdir=args.outDir

if len(files)!=len(filenames):
    raise NameError('Number of files and new file names disagree.')

make_sure_path_exists(outdir)


# Running for each bam file supplied
for bam, name in zip(files, filenames):
    start_time = time.time()
    bamid=str(int(random()*10000)) # id for the temp files

    #bam to sam
    infileHeader = pysam.AlignmentFile(bam, mode='rb').header
    infile = pysam.AlignmentFile(bam, mode='rb', header=infileHeader)
    outfile = pysam.AlignmentFile(name+"_temp"+bamid+".sam", "wh", 
        header=infileHeader, template=infile)
    for s in infile:
        outfile.write(s)
    infile.close()
    outfile.close()

    #get read 2 5' with read 1 directionality
    result=open(name+"_snr_temp"+bamid+".sam",'w')
    flags=["147", "163"]

    # Get coordinates
    with open(name+"_temp"+bamid+".sam") as infile:
        for line in infile:
            if not line.startswith("@"):
                line=line.strip("\n").split("\t")

                if line[1] in flags and line[6]=="=" and \
                    (not "I" in line[5] and not "D" in line[5]):

                    line[8]=str(int(line[8])/abs(int(line[8])))
                    #field 8 can be negative

                    if line[1]==flags[0]:
                        #second in pair (147), reverse
                        line[1]="99"
                        line=setStrand(line)
                        add_to_start=readLength_CIGAR(line[5])
                        line[3]=str(int(line[3])+add_to_start)
                        line[9]=line[9][-1]
                        line[10]=line[10][-1]
                        line[5]="1M"

                    elif line[1]==flags[1]:
                        #second in pair (163), forward
                        line[1]="83"
                        line=setStrand(line)
                        line[9]=line[9][0]
                        line[10]=line[10][0]
                        line[5]="1M"

                    else:
                        print line

                    result.write("\t".join(line)+"\n")
                    
            else: #write header
                result.write(line)

    result.close()
    

    
    #new sam to bam, then sort and index
    infile = pysam.AlignmentFile(name+"_snr_temp"+bamid+".sam", "r")
    outfile = pysam.AlignmentFile(name+"_snr_temp"+bamid+".bam", 
        "wb", template=infile)
    for s in infile:
        outfile.write(s)
    infile.close()
    outfile.close()

    sorted_file=pysam.sort(name+"_snr_temp"+bamid+".bam", 
        outdir+name+"_sorted")
    sorted_file=pysam.index(outdir+name+"_sorted.bam")
    
    removeTemp(bamid)

    print name+" done!"
    print("--- %s seconds ---" % round(time.time() - start_time, 2))

print "Finished!"
