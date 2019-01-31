# -Advanced reservoir
def QAdv(self, pcr):
    DayNo = self.timecalc.julian(self)[0]
    # -determine if it is flood or dry season
    S1 = pcr.ifthenelse(
        self.ResFlStart < self.ResFlEnd,
        pcr.ifthenelse(
            DayNo >= self.ResFlStart,
            pcr.ifthenelse(DayNo <= self.ResFlEnd, pcr.boolean(1), pcr.boolean(0)),
            pcr.boolean(0),
        ),
        pcr.ifthenelse(
            DayNo >= self.ResFlEnd,
            pcr.ifthenelse(DayNo >= self.ResFlStart, pcr.boolean(1), pcr.boolean(0)),
            pcr.ifthenelse(
                DayNo <= self.ResFlEnd,
                pcr.ifthenelse(
                    DayNo <= self.ResFlStart, pcr.boolean(1), pcr.boolean(0)
                ),
                pcr.boolean(0),
            ),
        ),
    )

    S_avail = pcr.max(self.StorRES - self.ResPVOL, 0)
    Q = pcr.max(
        pcr.ifthenelse(
            S1,
            self.ResMaxFl * S_avail / (self.ResEVOL - self.ResPVOL),
            self.ResDemFl * S_avail / (self.ResEVOL - self.ResPVOL),
        ),
        self.StorRES - self.ResEVOL,
    )
    return Q


# -Simple reservoir
def QSimple(self, pcr):
    Q = pcr.max(
        self.ResKr * self.StorRES * (self.StorRES / self.ResSmax) ** 1.5,
        self.StorRES - self.ResSmax,
    )
    return Q


# -Calculates reservoir outflow and the fraction to release, depending on the type of reservoir (simple or advanced)
def QRes(self, pcr):
    if self.ResSimple and self.ResAdvanced:
        Qout = pcr.ifthenelse(
            self.ResFunc == 1,
            QSimple(self, pcr),
            pcr.ifthenelse(self.ResFunc == 2, QAdv(self, pcr), 0),
        )
    elif self.ResSimple:
        Qout = pcr.ifthenelse(self.ResFunc == 1, QSimple(self, pcr), 0)
    else:
        Qout = pcr.ifthenelse(self.ResFunc == 2, QAdv(self, pcr), 0)

    return Qout
