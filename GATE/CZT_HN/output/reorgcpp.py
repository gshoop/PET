#!/usr/bin/python  

import sys

#fname = sys.argv[1]
#with open(fname,'r') as r:
#    for kk,line in enumerate(r):
#        l = line.split(' ')[1:]
#        print line
#        if kk == 1:
#            break
#with open('datatest.dat','w') as w:
#    w.write(' '+'123\n')
#    w.write(' '+'456')
fin = sys.argv[1]
fout = fin.replace('.lst','.txt')
tn=0
sn=0
with open(fin,'r') as r, open(fout,'w') as w:
    for kk,line in enumerate(r):
#        print line
#        if kk == 10:
#            break
        if kk==0:
            continue

        l = line.split()
        tn+=1
        x1=float(l[1])
        y1=float(l[2])
        z1=float(l[3])
        x2=float(l[12])
        y2=float(l[13])
        z2=float(l[14])
#        aa=(x2-x1)**2+(y2-y1)**2
#        bb=2.0*x1*(x2-x1)+2.0*y1*(y2-y1)
#        cc=x1**2+y1**2-10000.
#        delta=bb**2-4.*aa*cc
#        if delta<0:
#            continue
#        elif delta==0:
#            ts=-bb/(2.*aa)
#            z=z1+ts*(z2-z1)
#            if z>108. or z < -108.:
#                continue
#        else:
#            ts1=(-bb+delta**0.5)/(2.*aa)
#            ts2=(-bb-delta**0.5)/(2.*aa)
#            zz1=z1+ts1*(z2-z1)
#            zz2=z1+ts2*(z2-z1)
#            if (zz1>108. or zz1<-108.) and (zz2>108. or zz2<-108.):
#                continue

        sn+=1
#        w.write(l[1]+' '+l[2]+' '+l[3]+' '+l[12]+' '+l[13]+' '+l[14]+'\n')
        w.write(str("%.2f"%x1)+' '+str("%.2f"%y1)+' '+str("%.2f"%z1)+' '+str("%.2f"%x2)+' '+str("%.2f"%y2)+' '+str("%.2f"%z2)+'\n')
#print 'Total number of LOR is {0}. Number of LOR that is within FOV is {1}.'.format(tn,sn)
