#!/usr/bin/python

import sys

def fileNoExtension(s) :
	i = s.rfind('.')
	if i == -1 :
		return s
	return s[:i]

def filePath(fname, path=True) :
	i = fname.rfind('\\')
	if i == -1 : i = fname.rfind('/')
	if path == True :
		if i == -1 : return ''
		else : return fname[:i+1]
	else :
		if i == -1 : return fname
		else : return fname[i+1:]

def generateFileName(f1) :
	i = f1.rfind('.')
	if i == -1 :
		i = len(f1)
	return f1[:i] + '_SORTED' + f1[i:]

def gffMakeSortKey_value(x) :
	# for gff files, splits at '\t'
	pieces=x.split('\t')

	"""
	chrm=pieces[0]
	startp_s=str(pieces[3])
	startp_sform='0'*(9-len(startp_s)) + startp_s
	sortkey=chrm+'_'+startp_sform
	"""
	sortkey=float(pieces[5]) #make sure is float value
	return (sortkey)

def usage() :
	py = filePath(sys.argv[0],path=False)
	print '\nUsage:'
	print py, '<filename> [<new file>]'
	print '<filename> is the name of the file to be sorted'
	print '<new file> (optional) is the name of the sorted file being created\n'

if len(sys.argv) < 2 :
	usage()
	sys.exit()

fname = sys.argv[1]
if len(sys.argv) > 2 :
	newFile = sys.argv[2]
else :
	newFile = generateFileName(fname)

f = open(fname,'r')
lines = f.readlines()
f.close()

# separate comments
data = []
comments = []
[data.append(x) if x[0] != '#' else comments.append(x) for x in lines]
del lines

# sort gff
data.sort(key=gffMakeSortKey_value)


f = open(newFile,'w')
f.writelines(comments + data)
f.close()

print 'File created: ', newFile
