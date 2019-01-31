# -Function to rout the specific runoff
def ROUT(self, pcr, rvolume, qold, qout, sres):
    # Calculate the discharge Q (m3/d)
    Q = pcr.accufractionflux(self.FlowDir, rvolume, self.QFRAC)
    # Re-calculate Q, based on qold en kx, and assign Qout for cells being lake/reservoir
    Q = pcr.ifthenelse(
        self.QFRAC == 0, qout, (1 - self.kx) * Q + (qold * 24 * 3600) * self.kx
    )
    # Only calculate inflow for lake/reservoir cells
    Qin = pcr.ifthenelse(self.QFRAC == 0, pcr.upstream(self.FlowDir, Q), 0)
    sres = sres - qout + Qin
    Q = Q / (24 * 3600)  # -only convert Q to m3/s
    return sres, Q, Qin
