
"""
PROGRAM: convertwinpathtoR.py
GOAL: take a windows path on the command line, convert to an R-acceptable path
To do this, the program simply subsitutes / for \
it prints the result

"""

import sys

slash="\\"
original=sys.argv[1] #just takes string on command line

newone=original.replace(slash,"/")
newone2=newone.rstrip()
sys.stdout.write(newone2)


