print 'routing module imported'

def ROUT(pcr, q, oldq, flowdir, kx):
    rr = (q * 0.001 * pcr.cellarea()) / (24*3600)
    ra = pcr.accuflux(flowdir, rr)
    ra = (1 - kx) * ra + kx * oldq
    return ra
