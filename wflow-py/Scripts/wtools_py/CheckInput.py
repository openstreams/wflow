# -*- coding: utf-8 -*-
"""
Created on Tue Mar 17 19:08:52 2015

@author: tollena
"""

import wflowtools_lib as wt
import getopt
import os
import shutil
import sys

def Usage():
    print('')             
    print('Usage: CheckInput -h [-r river] ')
    print '-r   rivers to be checked for network (optional) (ESRI Shapefile)'
    print('')

def main():
    argv = sys.argv
    #argv = ['x','-r', 'river_selec.shp']
    
    
    rivershp= None
    
    resultdir = 'check\\'
    if os.path.isdir(resultdir):
        shutil.rmtree(resultdir)
    os.makedirs(resultdir)
    
    try:
        opts, args = getopt.getopt(argv[1:], 'hr:')
    except getopt.error:
        print 'fout'
        Usage()
        sys.exit(1)
    
    for o, a in opts:
        if o == '-h':
         Usage()
         sys.exit()        
        elif o == '-r': rivershp = a
        
    if rivershp == None:
        print 'input file or extent need to be specified'
        Usage()
        sys.exit(1)
    
    shapes = wt.Reach2Nodes(rivershp,4326,0.0001,resultdir)
if __name__ == "__main__":
    main()
