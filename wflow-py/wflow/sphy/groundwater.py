print "groundwater module imported"

# -Function to calculate groundwater recharge
def GroundWaterRecharge(pcr, deltagw, gwrecharge, subperc, glacperc):
    gwseep = (1 - pcr.exp(-1 / deltagw)) * (subperc + glacperc)
    gwrecharge = (pcr.exp(-1 / deltagw) * gwrecharge) + gwseep
    return gwrecharge


# -Function to calculate baseflow
def BaseFlow(pcr, gw, baser, gwrecharge, basethresh, alphagw):
    baser = pcr.ifthenelse(
        gw <= basethresh,
        0,
        (baser * pcr.exp(-alphagw) + gwrecharge * (1 - pcr.exp(-alphagw))),
    )
    return baser


# -Function to calculate the groundwater height, taken from the bottom of the gw layer (zero reference)
def HLevel(pcr, Hgw, alphagw, gwrecharge, yield_gw):
    Hgw = (Hgw * pcr.exp(-alphagw)) + (
        (gwrecharge * (1 - pcr.exp(-alphagw))) / (800 * yield_gw * alphagw)
    )
    return Hgw
