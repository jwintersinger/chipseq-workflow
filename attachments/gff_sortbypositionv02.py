
"""
PROGRAM: gff_sortbypositionv02.py
#Mark Bieda May 31, 2008 modified for better comments June 2009
#this simple program sorts a gff file and outputs a few lines
#a sort of test
don't know current status
INPUT:
see f line for file location

OUTPUT:
just a few lines to terminal; could be changed to file easily

TESTING:
none

"""

def makesortkey(x):
    pieces=x.split()
    chrom=pieces[0]
    startp_s=str(pieces[3])
    startp_sform='0'*(9-len(startp_s)) + startp_s
    sortkey=chrom+'_'+startp_sform
    return (sortkey)
    

# 1. open and read file
f=open(r'C:\Documents and Settings\user-clk\Desktop\CURRENT_SCIENCE\Dong Stuff\File1.gff.zfsort', 'r')
b=f.readlines()


# 2. remove comments
z=[]
for x in b:
    if x[0]<>'#':
        z.append(x)

del b
b=z

# 3. sort (note that b will be the ultimate sorted list)
b.sort(key=makesortkey)

# 4. print stuff - this is the semi-lame output
print
print 'second'
print b[0:5]



#OLDER STUFF

"""
justvals=[]
for x in b:
    #y=x.rstrip("\n")
    val=y.split()[3]
    val=float(val)
    justvals.append(val)
    total=total+val
    
#note h on the next line - change to b if necessary!!
print total/len(h)
print "minimum is ",
print min(justvals)
print "max is ",
print max(justvals)
#now sort the list to do the 95th percentile thing
justvals.sort()
njustvals=len(justvals)

#testing
p40=int(0.40*njustvals)-1
val40=justvals[p40]
print "40th percentile is ", val40
#commented out temp
"""

"""
p95=int(0.95*njustvals)-1
val95=justvals[p95]
print "95th percentile is ", val95
"""

