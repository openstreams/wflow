print("snow module imported")

# -Function to calculate the potential snow melt
def PotSnowMelt(pcr, temp, ddfs):
    melt = pcr.max(0, temp) * ddfs
    return melt


# -Function to calculate the actual snow melt
def ActSnowMelt(pcr, snowstore, potmelt):
    melt = pcr.min(snowstore, potmelt)
    return melt


# -Function that updates the snow storage
def SnowStoreUpdate(pcr, snowstore, snow, actmelt, temp, snowwatstore):
    snowstore = (
        snowstore
        + snow
        - actmelt
        + pcr.ifthenelse(temp < 0, pcr.scalar(snowwatstore), 0)
    )
    return snowstore


# -Function that determines the maximum amount of water that can be stored in the snowpack
def MaxSnowWatStorage(snowsc, snowstore):
    maxsnowwatstore = snowsc * snowstore
    return maxsnowwatstore


# -Function to calculate the actual snow water storage
def SnowWatStorage(pcr, temp, maxsnowwatstore, snowwatstore, actmelt, rain):
    snowwatstore = pcr.ifthenelse(
        temp < 0, 0, pcr.min(maxsnowwatstore, snowwatstore + actmelt + rain)
    )
    return snowwatstore


# -Function to calculate the total snow storage (snowstore + snowwatstore)
def TotSnowStorage(snowstore, snowwatstore, snowfrac, rainfrac):
    totalsnowstore = (snowstore + snowwatstore) * (snowfrac + rainfrac)
    return totalsnowstore


# -Function to calculate runoff from snow
def SnowR(pcr, snowwatstore, maxsnowwatstore, actmelt, rain, oldsnowwatstore, snowfrac):
    snowr = pcr.ifthenelse(
        snowwatstore == maxsnowwatstore,
        (((actmelt + rain) - (snowwatstore - oldsnowwatstore)) * snowfrac),
        0,
    )
    return snowr
