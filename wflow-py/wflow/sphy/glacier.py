print('glacier module imported')

# -Function to calculate melt from clean ice or debris covered glaciers
def GlacCDMelt(pcr, temp, ddf, glacfrac):
    glacdmelt = pcr.max(0, temp) * ddf * glacfrac
    return glacdmelt


# -Total glacier melt
def GMelt(glaccimelt, glacdcmelt):
    glacmelt = glaccimelt + glacdcmelt
    return glacmelt


# -Function to calculate runoff from glaciers
def GlacR(glacf, gmelt, glacfrac):
    glacr = glacf * gmelt * glacfrac
    return glacr


# -Function to calculate glacier percolation to groundwater
def GPerc(glacf, gmelt, glacfrac):
    gperc = (1 - glacf) * gmelt * glacfrac
    return gperc
