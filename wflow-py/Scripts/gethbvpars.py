# -*- coding: utf-8 -*-
"""

gethbvpars:
    gets the HBV catchment parameters from a hbv model. It assules the default
    parameters are stored in the root (the basin) and that tuned parameters
    are stored in each catchment.
    
    syntax:
        gethbvpars -p pathtobasin -o outputfilename 

"""

import os.path
import getopt
import sys


sep = ","
csvfile = "test.csv"


def usage(*args):
    sys.stdout = sys.stderr
    for msg in args:
        print(msg)
    print(__doc__)
    sys.exit(0)


def readpar(fname, skip):
    a = {}
    f = open(fname, "rb")
    if skip:
        x = f.readline()
    x = f.readlines()
    f.close()

    for l in x:
        ll = [c for c in l if c not in "'"].split()
        if len(ll) > 0:
            a[ll[0]] = ll[1]

    return a


def readbas(fname):
    a = []
    f = open(fname, "rb")
    x = f.readline()
    x = f.readlines()
    f.close()

    for l in x:
        ll = [c for c in l if c not in "'\\"].split()
        if len(ll) > 0:
            if ll[0] == "basindir":
                a.append(ll[1])

    return a


basin = ""
catch = {}

try:
    opts, args = getopt.getopt(sys.argv[1:], "o:p:h")
except getopt.error as msg:
    usage(msg)

for o, a in opts:
    if o == "-p":
        basin = a
    if o == "-o":
        csvfile = a
    if o == "-h":
        usage()


# read basin structure and order
basstruc = readbas(basin + "/basin.par")
# read default parameters
baspar = readpar(basin + "/rmod.par", 0)

for ddri in basstruc:
    pfile = basin + "/" + ddri + "/bmod.par"
    if os.path.exists(pfile):
        xx = readpar(pfile, 1)
        catch[os.path.basename(ddri)] = xx


f = open(csvfile, "w")
i = 0
print("Id,Name", end=" ", file=f)
for ppar in baspar:
    print(sep + ppar, end=" ", file=f)
print("", file=f)


# for c in catch:
for ii in range(0, len(basstruc) - 1):
    i = i + 1
    c = basstruc[ii]
    print(str(i) + sep + c, end=" ", file=f)
    for ppar in baspar:
        if ppar in catch[c]:
            print(sep + catch[c][ppar], end=" ", file=f)
        else:
            print(sep + baspar[ppar], end=" ", file=f)
    print("", file=f)


f.close()
