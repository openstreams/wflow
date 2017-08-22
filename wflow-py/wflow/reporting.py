#-Function to report the output
def REPM(self, pcr, tot, var, fname, outops, TSS=False, MAP=False):
    if outops == 'Day':
        if TSS:
            TSS.sample(var)
        if MAP:
            self.report(var, self.outpath + fname)
        tot = 0
    elif outops == 'Month':
        dim = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
        if self.calendar.isleap(self.curdate.year):
            dim[1] = 29
        else:
            dim[1] = 28
        tot = tot + var
        if self.curdate.day == dim[self.curdate.month-1]:
            if TSS:
                TSS.sample(tot)
            if MAP:
                self.report(tot, self.outpath + fname + 'M')
            tot = 0
    elif outops == 'Year':
        if self.calendar.isleap(self.curdate.year):
            ydays = 366
        else:
            ydays = 365
        tot = tot + var
        if self.timecalc.julian(self)[0] == ydays:
            if TSS:
                TSS.sample(tot)
            if MAP:
                self.report(tot, self.outpath + fname + 'Y')
            tot = 0
    else:
        tot = tot + var
        if self.curdate == self.enddate:
            pcr.report(tot, self.outpath + fname + '.map')
            tot = 0
    return tot
 
#-Function to initialise the reporting
def reporting(self, pcr, tot, var):
    for outops in ['Day','Month','Year','Final']:
        try:
            TSS = eval('self.' + tot + '_' + outops + 'TS')
            try:
                MAP = eval('self.' + tot + '_' + outops + '_map')
                setattr(self, tot + '_'+outops, REPM(self, pcr, eval('self.'+tot+'_'+outops), var, eval('self.'+tot+'_fname'), outops, TSS, MAP))
            except:
                setattr(self, tot + '_'+outops, REPM(self, pcr, eval('self.'+tot+'_'+outops), var, eval('self.'+tot+'_fname'), outops, TSS))
        except:
            try:
                MAP = eval('self.' + tot + '_' + outops + '_map')
                setattr(self, tot + '_'+outops, REPM(self, pcr, eval('self.'+tot+'_'+outops), var, eval('self.'+tot+'_fname'), outops, False, MAP))
            except:
                pass    
    
