
"""
PROGRAM: macsbed2gffv05.py
GOAL: will convert macs bed format to gff format. Note macs format is like UCSC format, only fairly short.
INVOCATION: python macsbed2gffv05.py infile outfile
infile and outfile should both have full path (if necessary)
infile is MACS program output bed format (.bed)
outfile is GFF format (.gff)
MARK BIEDA

STATUS COMMAND LINE IMPLEMENTED
TESTING: minor testing

EXAMPLE
***note that python must be invoked from commandline***
This will not directly find python
SO:
C:\Python25\python.exe macsbed2gffv05.py testfile.bed gff_file.gff

September 16, 2012
September 14, 2012
August 29, 2012
AUGUST 13, 2008

"""

#import os
import sys


def macsbed2gff(bline):
	#this is for macs format which is ucsc format
	bline2=bline.rstrip()
	bp=bline2.split("\t")
	
	descrip=bp[3] #add other stuff from bed onto this
	bp[1]=str(long(bp[1]) + 1) #bed first coord is zero-based from UCSC - see UCSC def. so must adjust by one
	gffline="\t".join([bp[0],"MACS","macsbedline2gff",bp[1],bp[2],bp[4],"+","0",descrip])
	gffline += '\n'
	return(gffline)


#mydir=r'C:\Users\markb\Desktop\CURRENT_PROJECTS\Kepler2\testdata'
#fulldir=mydir + os.sep



#f = open(fulldir+"broad2bed.bed", "r")
#outf = open(fulldir+"bed2gff.gff","w")

f = open(sys.argv[1], "r")
outf = open(sys.argv[2],"w")
for line in f:
        if (line[0:5] != "track"):
                gffline=macsbed2gff(line)
                outf.write(gffline)

f.close()
outf.close()

"""

#for getting from command line
#probably must be near top of script

import sys

myargs=sys.argv

infilenm=myargs[1]
outfilenm=myargs[2]

"""

"""

#this is approx code for handling a file...
f = open("myFile.bed", "r")
    for line in f:
		gffline=macsbedline2gffline(line)
		print gffline
		
#maybe this but not sure of syntax
#need python black book!
f = open("myFile.bed", "r")
outf = open("outfile.gff","w")
    for line in f:
		gffline=macsbedline2gffline(line)
		outf.write(gffline)


"""

