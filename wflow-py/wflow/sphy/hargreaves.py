#-Extraterrestrial radiation    
def extrarad(self, pcr):
    DayNo = self.wf_supplyJulianDOY() # timecalc.julian(self)[0]
    LatRad = self.Lat * (self.pi / 180)
    dr = 1 + 0.033 * pcr.cos((2 * self.pi * DayNo) /  365)
    delta = 0.409 * pcr.sin(((2 * self.pi * DayNo) / 365) - 1.39)
    omegas = pcr.acos(-1 * pcr.tan(LatRad) * pcr.tan(delta))
    Ra = ((24 * 60) / self.pi) * self.Gsc * dr * (pcr.scalar(omegas) * pcr.sin(LatRad) * pcr.sin(delta) +\
        pcr.cos(LatRad) * pcr.cos(delta) * pcr.sin(omegas))
    return Ra    

#-Modified Hargreaves for calculation of ETref
def Hargreaves(pcr, ra, temp, tempmax, tempmin):
    ETref = pcr.max(0.0023 * 0.408 * ra * (temp + 17.8) * (pcr.max(tempmax - tempmin, 0))**0.5, 0)
    return ETref