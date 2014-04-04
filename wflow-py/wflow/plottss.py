#!/usr/bin/python

#
# Wflow is Free software, see below:
# 
# Copyright (c) J. Schellekens 2005-2013
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

Open an interactive plot of one or more tss files
syntax:

    plottss.py [-T title][-L][-S][-C colums] file1 file2 filen .....
    
    -L if specified omit the legend
    -S make a subplot for each file
    -C specify colums to plot (e.g. -C "1,2,3,6" plots column 1,2,3 and 6)

"""

try:
    from wflow.pcrut import *
except:
    from pcrut import *
import getopt
from pylab import *



def usage(*args):
    sys.stdout = sys.stderr
    for msg in args: print msg
    print __doc__
    sys.exit(0)



def main(argv=None):
    
    
    makelegend = True
    subplots = False
    cols =":"
    
    if argv is None:
        argv = sys.argv[1:]
        if len(argv) == 0:
            usage()
            return     
    
    plottitle = ""
    opts, args = getopt.getopt(argv, 'T:LSC:')
    
    for o, a in opts:
        if o == '-T': plottitle = a  
        if o == '-L': makelegend = False
        if o == '-S': subplots = True
        if o == '-C': cols = a


    nrplotfiles = len(args)
    
    if subplots:
        # Two subplots, the axes array is 1-d
        f, axarr = plt.subplots(nrplotfiles, sharex=True)
        plotnr = 0
        for plotfile in args:    
            dat, x = readtss(plotfile)
            axarr[plotnr].plot(eval("dat[" + cols +"]"),label=plotfile)

            if makelegend:
                axarr[plotnr].legend()
            plotnr = plotnr + 1
    else:
        f, axarr = plt.subplots()
        for plotfile in args:    
            dat, x = readtss(plotfile)
            axarr.plot(eval("dat[" + cols +"]"),label=plotfile)
        if makelegend:
            legend()
        
    title(plottitle)
    xlabel("Time")
    ylabel("Value")
    show()
    
if __name__ == "__main__":
    main()