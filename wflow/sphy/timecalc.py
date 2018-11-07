# -Function to return the julian day of the year
def julian(self):
    y = self.curdate.year
    start = self.datetime.datetime(y, 1, 1).toordinal()
    current = self.curdate.toordinal()
    day = current - start + 1
    return day, 1


# -Function to calculate the number of timesteps for the model run
def timesteps(self):
    nrTimeSteps = (self.enddate - self.startdate).days + 1
    print(
        "Running SPHY for "
        + str(self.startdate.day)
        + "-"
        + str(self.startdate.month)
        + "-"
        + str(self.startdate.year)
        + " through "
        + str(self.enddate.day)
        + "-"
        + str(self.enddate.month)
        + "-"
        + str(self.enddate.year)
    )
    print("with " + str(nrTimeSteps) + " daily timesteps")
    return nrTimeSteps
