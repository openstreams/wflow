# -Function to calculate the potential evapotranspiration
def ETpot(etr, kc):
    etpot = etr * kc
    return etpot


# -Function to calculate the actual evapotranspiration
def ETact(pcr, etpot, rootwater, rootsat, etreddry, rainfrac):
    etredwet = pcr.ifthenelse(rootwater >= rootsat, pcr.scalar(0), 1)
    etact = pcr.ifthenelse(
        rainfrac > 0, pcr.min(etpot * etreddry * etredwet, rootwater), 0
    )
    return etact
