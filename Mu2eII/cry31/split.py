#!/usr/bin/python
#------------------------------------------------------------------------------
# call: 
# python split.py /pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cosm0s00b0.s4_helix_filter/outstage/35162167/00
#
# prints a list, splits input by ~2000 events per job
#------------------------------------------------------------------------------

import sys, glob

#------------------------------------------------------------------------------
def split(fn):
    nmax = 2000;

    lines = open(fn).readlines()

    nev = len(lines);

    # print(" nev = ",nev)

    ev = {}
    for l in lines:
        sr = int(l.split()[8])
        if ( sr in ev.keys()):
            ev[sr] += 1;
        else:
            ev[sr] = 1
    
    #    ev.sort()

    # print(ev)
    # assume sorted

    n0   = -1;
    i    = 0;
    sr1  = -1

    #    print('keys:', ev.keys())

    # print('keys sorted:', sorted(ev.iterkeys(),key=int))

    xx = []

    for sr in ev.keys():
        if (n0 < 0):
            sr0 = sr
            n0  = 0

        if (n0 > nmax):
            # print ('i,sr0,n0 = ',i,sr0,n0)
            xx.append((sr0,n0))
            sr0 = sr
            n0  = 0
            i   = i+1;

        n0 = n0 + ev[sr];

    if (n0 > 0):
        # print ('i,sr1,n0 = ', i,sr0,n0);
        xx.append((sr0,n0))
            
    # print('xx : ',xx)

    print('("%s",  %5i, '%(fn, nev), end='')

    for i in range(0,len(xx)):
        print("%5i, %5i"%(int(xx[i][0]),int(xx[i][1])),end='');
        if (i < len(xx)-1):
            print(', ',end='')

    print('),')

#------------------------------------------------------------------------------
# main program, just make a GridSubmit instance and call its methods
#------------------------------------------------------------------------------
if (__name__ == '__main__'):
    
    # dir = '/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cosm0s00b0.s4_helix_filter/outstage/35162167/00'

    dir = sys.argv[1]
    
    sd_list = glob.glob(dir+'/*')
    sd_list.sort()

    # print(sd_list)
    # fn  = dir+'/00000/'+'dig.mu2e.cosm0s41b0.Mu2eII.001002_00000000.art.event_list'
    for d in sd_list:
        fn = glob.glob(d+'/*.event_list')[0]
        # print(fn)
        split(fn)

    sys.exit(0);
