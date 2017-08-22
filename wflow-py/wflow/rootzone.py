#-Function to calculate rootzone runoff
def RootRunoff(pcr, rainfrac, rootwater, rootsat):
    rootrunoff = pcr.ifthenelse(rainfrac > 0, pcr.max(rootwater - rootsat, 0), 0)
    return rootrunoff

#-Function to calculate rootzone drainage
def RootDrainage(pcr, rootwater, rootdrain, rootfield, rootsat, drainvel, rootTT):
    rootexcess = pcr.max(rootwater - rootfield, 0)
    rootexcessfrac = rootexcess / (rootsat - rootfield)
    rootlat = rootexcessfrac * drainvel
    rootdrainage = pcr.max(pcr.min(rootwater, rootlat * (1-pcr.exp(-1/rootTT)) + rootdrain * pcr.exp(-1/rootTT)), 0)
    return rootdrainage

#-Function to calculate rootzone percolation
def RootPercolation(pcr, rootwater, subwater, rootfield, rootTT, subsat):
    rootexcess = pcr.max(rootwater - rootfield, 0)
    rootperc = rootexcess * (1 - pcr.exp(-1 / rootTT))
    rootperc = pcr.ifthenelse(subwater >= subsat, 0, pcr.min(subsat - subwater, rootperc))
    rootperc = pcr.max(pcr.min(rootperc, rootwater), 0)
    return rootperc
