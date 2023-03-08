# -Function that returns crop factor (Kc) and maximum storage (Smax)
def Veg_function(
    pcr, ndvi, fpar_max, fpar_min, lai_max, ndvi_min, ndvi_max, kc_min, kc_max
):
    SR = (1 + ndvi) / (1 - ndvi)
    SR_max = (1 + ndvi_max) / (1 - ndvi_max)
    SR_min = (1 + ndvi_min) / (1 - ndvi_min)
    FPAR = pcr.min(
        (((SR - SR_min) * (fpar_max - fpar_min)) / (SR_max - SR_min)) + 0.001, 0.95
    )
    LAI = lai_max * pcr.log10(1 - FPAR) / pcr.log10(1 - fpar_max)
    Smax = 0.935 + 0.498 * LAI - 0.00575 * (LAI**2)
    Kc = kc_min + (kc_max - kc_min) * pcr.max(
        pcr.min((ndvi - ndvi_min) / (ndvi_max - ndvi_min), 1), 0
    )
    return Kc, Smax


# -Function that returns the interception, precipitation throughfall, and remaining storage
def Inter_function(pcr, S, Smax, Etr):
    PreT = pcr.max(0, S - Smax)
    S = S - PreT
    Int = pcr.min(1.5 * Etr, S)
    S = S - Int
    return Int, PreT, S
