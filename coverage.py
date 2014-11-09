#!/usr/bin/python
###########################################################################
# filename      : coverage.py
# date written  : 01/08/2013
# written by    : THChew (teonghan@gmail.com)
# description   : check fasta and discard any with value more than
#                 specified
# last update   : 01/08/2013    - version 0.1
#                 30/10/2013    - version 0.2 coverage + length filtering
###########################################################################
# split length function
# original author: Abu Zahed Jony
def split_by_length(s,block_size):
    w=[]
    n=len(s)
    for i in range(0,n,block_size):
        w.append(s[i:i+block_size])
    return w

import sys
from Bio import SeqIO

if len(sys.argv)!=int(6):
    print ""
    print "Invalid argument"
    print "Usage: python coverage.py fasta_file lower_value upper_value lower_length upper_length"
    print ""
    print "lower_value  : lowest coverage value to keep (exclusive), -1 if no lower_value"
    print "upper_value  : highest coverage value to keep (exclusive), -1 if no upper_value"
    print "lower_length : lowest length value to keep (exclusive), -1 if no lower_length"
    print "upper_length : highest length to keep (exclusive), -1 if no upper_length"
    print ""
    print "Example      : python extract.py input 11 13 40000 -1"
    print ""
    sys.exit()

else:

    fasta=str(sys.argv[1])
    
    try:
        lower_value=float(sys.argv[2])
    except ValueError:
        print 'lower_value is not a number'
        print ""
        sys.exit()

    try:
        upper_value=float(sys.argv[3])
    except ValueError:
        print 'upper_value is not a number'
        print ""
        sys.exit()

    if lower_value == -1 and upper_value==-1:
        print 'warning : coverage value will not be filtered (lower_value = upper_value = -1)'
    
    elif lower_value!=-1 and upper_value!=-1 and upper_value < lower_value:
        print 'value error : lower_value is larger than upper_value'
        print ""
        sys.exit()
        
    elif lower_value!=-1 and upper_value!=-1 and upper_value == lower_value:
        print 'value error : lower_value is equal to upper_value'
        print ""
        sys.exit()

    try:
        lower_length=float(sys.argv[4])
    except ValueError:
        print 'lower_length is not a number'
        print ""
        sys.exit()

    try:
        upper_length=float(sys.argv[5])
    except ValueError:
        print 'upper_length is not a number'
        print ""
        sys.exit()

    if lower_length == -1 and upper_length==-1:
        print 'warning : length value will not be filtered (lower_length = upper_length = -1)'

    elif lower_length!=-1 and upper_length!=-1 and upper_length < lower_length:
        print 'value error : lower_length is larger than upper_length'
        print ""
        sys.exit()
        
    elif lower_length!=-1 and upper_length!=-1 and upper_length == lower_length:
        print 'value error : lower_length is equal to upper_length'
        print ""
        sys.exit()
        
    out=fasta+'.out'

    #print lower_value
    #print upper_value
    #print lower_length
    #print upper_length
    
    tmp=[]
    
    for seq_record in SeqIO.parse(fasta, "fasta"):
        name=seq_record.name
        START=name.index('cov_')
        END=name[START+len('cov_'):].index('_')+START+len('cov_')
        VALUE=float(name[START:END].split('_')[1])

        START=name.index('length_')
        END=name[START+len('length_'):].index('_')+START+len('length_')
        LENGTH=float(name[START:END].split('_')[1])

        name='>'+name
        sequence=seq_record.seq.tostring()
        sequence=split_by_length(sequence,80)

        tmp.append([name,sequence,VALUE,LENGTH])

    for i in range(0,len(tmp)):
        VALUE=tmp[i][2]
        
        if lower_value==-1 and upper_value==-1:
            pass
                
        elif lower_value==-1 and upper_value!=-1:
            if not VALUE < upper_value:
                tmp[i]=['']

        elif lower_value!=-1 and upper_value==-1:
            if not VALUE > lower_value:
                tmp[i]=['']
                
        elif lower_value!=-1 and upper_value!=-1 and upper_value > lower_value:
            if not lower_value < VALUE < upper_value:
                tmp[i]=['']

    for i in range(0,len(tmp)):
        if tmp[i]!=['']:
            VALUE=tmp[i][3]
        
            if lower_length==-1 and upper_length==-1:
                pass
                    
            elif lower_length==-1 and upper_length!=-1:
                if not VALUE < upper_length:
                    tmp[i]=['']

            elif lower_length!=-1 and upper_length==-1:
                if not VALUE > lower_length:
                    tmp[i]=['']
                    
            elif lower_length!=-1 and upper_length!=-1 and upper_length > lower_length:
                if not lower_length < VALUE < upper_length:
                    tmp[i]=['']

f=open(out,"w")
for i in tmp:
    if i!=['']:
        f.write(i[0]+'\n')
        f.write('\n'.join(i[1])+'\n')
f.close()
