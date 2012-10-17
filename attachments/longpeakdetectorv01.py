"""
PROGRAM: longpeakdetectorv01.py
GOAL: detects long peaks based on absolute length of peak
INVOCATION: python longpeakdetectorv01.py infile outfile size
infile and outfile should both have full path (if necessary)
infile is gff format
outfile is  gff format (.gff)
size is absolute integer value
MARK BIEDA

STATUS 
TESTING: NONE

EXAMPLE
***note that python must be invoked from commandline***
This will not directly find python
SO:
python longpeakdetectorv01.py lin9peaks.bed.gff line9peaks.bed.gff.longpeaks.gff 10000

September 16, 2012


"""

import sys

f = open(sys.argv[1], "r")
outf = open(sys.argv[2],"w")
minsize=int(sys.argv[3])
outlines=[]
for line in f:
        pieces=line.split("\t")
	totsize=1+ long(pieces[4])-long(pieces[3])
	if totsize>=minsize:
                outf.write(line)

f.close()
outf.close()
