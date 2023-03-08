# -Function that updates the lake storage and lake level given a measured lake level. If no lake
# level is measured, then the actual storage is not updated with a measured level. The function
# returns the updated storage and lake level
def UpdateLakeHStore(self, pcr, pcrm):
    # -buffer actual storage
    OldStorage = self.StorRES
    # -Check if measured lake levels area available
    try:
        LakeLevel = pcr.readmap(pcrm.generateNameT(self.LLevel, self.counter))
        Level = True
    except:
        Level = False
    if Level:
        # -update lake storage according to measured water level
        self.StorRES = pcr.ifthenelse(
            self.UpdateLakeLevel,
            pcr.ifthenelse(
                pcr.defined(LakeLevel),
                pcr.ifthenelse(
                    self.LakeSH_Func == 1,
                    self.LakeSH_exp_a * pcr.exp(self.LakeSH_exp_b * LakeLevel),
                    pcr.ifthenelse(
                        self.LakeSH_Func == 2,
                        self.LakeSH_pol_a1 * LakeLevel + self.LakeSH_pol_b,
                        pcr.ifthenelse(
                            self.LakeSH_Func == 3,
                            (self.LakeSH_pol_a2 * LakeLevel**2)
                            + self.LakeSH_pol_a1 * LakeLevel
                            + self.LakeSH_pol_b,
                            (self.LakeSH_pol_a3 * LakeLevel**3)
                            + (self.LakeSH_pol_a2 * LakeLevel**2)
                            + (self.LakeSH_pol_a1 * LakeLevel + self.LakeSH_pol_b),
                        ),
                    ),
                ),
                self.StorRES,
            ),
            self.StorRES,
        )
        # prevent storage becoming negative for whatever reason
        self.StorRES = pcr.max(self.StorRES, 0)
        # -Update the lake level based on the storage for lakes where no levels are measured
        LakeLevel = pcr.ifthenelse(
            self.UpdateLakeLevel,
            pcr.ifthenelse(
                pcr.defined(LakeLevel),
                LakeLevel,
                pcr.ifthenelse(
                    self.LakeHS_Func == 1,
                    self.LakeHS_exp_a * pcr.exp(self.LakeHS_exp_b * self.StorRES),
                    pcr.ifthenelse(
                        self.LakeHS_Func == 2,
                        self.LakeHS_pol_a1 * self.StorRES + self.LakeHS_pol_b,
                        pcr.ifthenelse(
                            self.LakeHS_Func == 3,
                            (self.LakeHS_pol_a2 * self.StorRES**2)
                            + self.LakeHS_pol_a1 * self.StorRES
                            + self.LakeHS_pol_b,
                            (self.LakeHS_pol_a3 * self.StorRES**3)
                            + (self.LakeHS_pol_a2 * self.StorRES**2)
                            + self.LakeHS_pol_a1 * self.StorRES
                            + self.LakeHS_pol_b,
                        ),
                    ),
                ),
            ),
            pcr.ifthenelse(
                self.LakeHS_Func == 1,
                self.LakeHS_exp_a * pcr.exp(self.LakeHS_exp_b * self.StorRES),
                pcr.ifthenelse(
                    self.LakeHS_Func == 2,
                    self.LakeHS_pol_a1 * self.StorRES + self.LakeHS_pol_b,
                    pcr.ifthenelse(
                        self.LakeHS_Func == 3,
                        (self.LakeHS_pol_a2 * self.StorRES**2)
                        + self.LakeHS_pol_a1 * self.StorRES
                        + self.LakeHS_pol_b,
                        (self.LakeHS_pol_a3 * self.StorRES**3)
                        + (self.LakeHS_pol_a2 * self.StorRES**2)
                        + self.LakeHS_pol_a1 * self.StorRES
                        + self.LakeHS_pol_b,
                    ),
                ),
            ),
        )

    else:
        # if no lake level map is available, then calculate the h based on storages
        LakeLevel = pcr.ifthenelse(
            self.LakeHS_Func == 1,
            self.LakeHS_exp_a * pcr.exp(self.LakeHS_exp_b * self.StorRES),
            pcr.ifthenelse(
                self.LakeHS_Func == 2,
                self.LakeHS_pol_a1 * self.StorRES + self.LakeHS_pol_b,
                pcr.ifthenelse(
                    self.LakeHS_Func == 3,
                    (self.LakeHS_pol_a2 * self.StorRES**2)
                    + self.LakeHS_pol_a1 * self.StorRES
                    + self.LakeHS_pol_b,
                    (self.LakeHS_pol_a3 * self.StorRES**3)
                    + (self.LakeHS_pol_a2 * self.StorRES**2)
                    + self.LakeHS_pol_a1 * self.StorRES
                    + self.LakeHS_pol_b,
                ),
            ),
        )
    self.StorRES = pcr.ifthenelse(self.LakeID != 0, self.StorRES, OldStorage)
    return LakeLevel, self.StorRES


# -function that calculates the fraction of lake storage that is available for routing, and the lake outflow
def QLake(self, pcr, LakeLevel):
    Qout = pcr.ifthenelse(
        self.LakeQH_Func == 1,
        self.LakeQH_exp_a * pcr.exp(self.LakeQH_exp_b * LakeLevel),
        pcr.ifthenelse(
            self.LakeQH_Func == 2,
            self.LakeQH_pol_a1 * LakeLevel + self.LakeQH_pol_b,
            pcr.ifthenelse(
                self.LakeQH_Func == 3,
                (self.LakeQH_pol_a2 * LakeLevel**2)
                + self.LakeQH_pol_a1 * LakeLevel
                + self.LakeQH_pol_b,
                (self.LakeQH_pol_a3 * LakeLevel**3)
                + (self.LakeQH_pol_a2 * LakeLevel**2)
                + self.LakeQH_pol_a1 * LakeLevel
                + self.LakeQH_pol_b,
            ),
        ),
    )
    Qout = pcr.max(0, Qout)
    Qout = Qout * 3600 * 24  # -convert to m3/d
    Qout = pcr.cover(Qout, 0)  # -for non-lake cells, Qout is zero
    return Qout
