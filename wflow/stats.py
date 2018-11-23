#!/usr/local/bin/python
#
# Created on July 10, 2000
#  by Keith Cherkauer
#
# This python script computes several standard statistics on arrays
# of values
#
# Functions include:
#  get_mean
#  get_median
#  get_var
#  get_stdev
#  get_skew
#  get_sum
#  get_min
#  get_max
#  get_count_over_threshold
#  get_quantile
#  get_running_average
#  get_running_slope
#  get_bias
#  get_root_mean_square
#  get_mean_absolute_error
#  get_max_absolute_error
#  get_nash_sutcliffe
#  get_peak_diff
#  get_number_of_sign_changes
#  get_peak_threshold_diff
#  get_covariance
#  get_correlation
#  get_cross_correlation
#  filter_threshold
#  get_days
#  get_last_day
#  get_first_day
#  get_box_plot_parameters
#


from math import fabs
from math import sqrt

from numpy import isnan

NoDataVal = -999
SmallValue = 1.0e-10


def get_mean(values, N="", NoData=NoDataVal, Skip=""):
    """This function computes the mean or average of an array of values
       after filtering out the NoData values.  It returns both the mean
       value and the number of valid data points used.  The mean value
       is set to the NoData value if there are no valid data points.
       If Skip is set, then values equal to skip are included in the count
       of active data, but not included in the calculation of the mean.
       An example of when this would be used is for computing average
       snow cover, where 0 indicates a valid measurement but does not
       contribute to a meaningful measurement of snow cover."""
    if not N:
        N = len(values)
    mean = 0
    Nact = 0
    Nskip = 0
    for i in range(N):
        if values[i] != NoData and not np.isnan(values[i]):
            if Skip and values[i] == Skip:
                Nskip = Nskip + 1
            else:
                mean = mean + values[i]
            Nact = Nact + 1
    if Nact - Nskip > 0:
        mean = mean / (Nact - Nskip)
    else:
        mean = NoData
    return (mean, Nact)


def get_median(values, N="", NoData=NoDataVal):
    """This function computes the median of an array of values
       after filtering out the NoData values.  It returns both the median
       value and the number of valid data points used.  The median value
       is set to the NoData value if there are no valid data points."""
    if not N:
        N = len(values)
    new_value = []
    Nact = 0
    for i in range(N):
        if values[i] != NoData and not np.isnan(values[i]):
            new_value = new_value + [values[i]]
            Nact = Nact + 1
    if Nact > 0:
        new_value.sort()
        if Nact % 2 == 0:
            median = (new_value[int(Nact / 2)] + new_value[int(Nact / 2)]) / 2.0
        else:
            median = new_value[int(Nact / 2)]
    else:
        median = NoData
    return (median, Nact)


def get_var(values, N="", mean="", NoData=NoDataVal):
    """This function computes the variance of an array of values after
       filtering out the NoData values.  The mean value of the array
       must be provided to the routine.  It returns both the variance
       value and the number of valid data points used.  The variance
       is set to the NoData value if there are no valid data points."""
    if not N:
        N = len(values)
    if not mean:
        mean = get_mean(values, N, NoData)[0]
    var = 0
    Nact = 0
    for i in range(N):
        if values[i] != NoData and not np.isnan(values[i]):
            var = var + (values[i] - mean) * (values[i] - mean)
            Nact = Nact + 1
    if Nact > 1:
        var = var / (Nact - 1)
    else:
        var = NoData
    return (var, Nact)


def get_stdev(values, N="", mean="", NoData=NoDataVal):
    """This function computes the standard deviation of an array of
       values after filtering out the NoData values.  The mean of the
       array must be provided to the routine.  It returns both
       the standard deviation value and the number of valid data
       points used.  The standard deviation is set to the NoData value
       if there are no valid data points."""
    if not N:
        N = len(values)
    if not mean:
        mean = get_mean(values, N=N, NoData=NoData)[0]
    stdev = 0
    Nact = 0
    for i in range(N):
        if values[i] != NoData and not np.isnan(values[i]):
            stdev = stdev + (values[i] - mean) * (values[i] - mean)
            Nact = Nact + 1
    if Nact > 1:
        stdev = stdev / (Nact - 1)
        stdev = math.sqrt(stdev)
    else:
        stdev = NoData
    return (stdev, Nact)


def get_skew(values, N="", mean="", stdev="", NoData=NoDataVal):
    """This function computes the skewness of an array of values after
       filtering out the NoData values.  The mean and standard deviation
       of the array must be provided to the routine.  It returns both
       the skewness value and the number of valid data points used.  The
       skewness is set to the NoData value if there are no valid data
       points."""
    if not N:
        N = len(values)
    if not mean:
        mean = get_mean(values, N, NoData)[0]
    if not stdev:
        stdev = get_stdev(values, N, mean, NoData)[0]
    skew = 0
    Nact = 0
    for i in range(N):
        if values[i] != NoData and not np.isnan(values[i]):
            skew = skew + (values[i] - mean) ** 3
            Nact = Nact + 1
    if (stdev ** 3 * (Nact - 1) * (Nact - 2)) != 0:
        skew = (skew * Nact) / (stdev ** 3 * (Nact - 1) * (Nact - 2))
    else:
        skew = NoData
    return (skew, Nact)


def get_sum(values, N="", NoData=NoDataVal):
    """This function computes the sum of an array of values after
       filtering out the NoData values.  It returns both the sum value
       and the number of valid data points used.  The sum is set to
       the NoData value if there are no valid data points."""
    if not N:
        N = len(values)
    sum = 0
    Nact = 0
    for i in range(N):
        if values[i] != NoData and not np.isnan(values[i]):
            sum = sum + values[i]
            Nact = Nact + 1
    if Nact == 0:
        sum = NoData
    return (sum, Nact)


def get_min(values, N="", NoData=NoDataVal):
    """This function finds the minimum value of an array after
       filtering out the NoData values.  It returns both the
       minimum value and the number of valid data points used.
       The minimum is set to the NoData value if there are no
       valid data points."""
    if not N:
        N = len(values)
    pos = 0
    while pos < N and values[pos] == NoData:
        pos = pos + 1
    if pos < N:
        Nact = 1
        min = values[pos]
        minpos = pos
        for i in range(pos, N):
            if values[i] != NoData and not np.isnan(values[i]):
                if values[i] < min:
                    min = values[i]
                    minpos = i
                Nact = Nact + 1
        if Nact == 0:
            min = NoData
    else:
        min = NoData
        minpos = NoData
        Nact = 0
    return (min, Nact, minpos)


def get_max(values, N="", NoData=NoDataVal):
    """This function finds the maximum value of an array after
       filtering out the NoData values.  It returns both the
       maximum value and the number of valid data points used.
       The maximum is set to the NoData value if there are no
       valid data points."""
    if not N:
        N = len(values)
    pos = 0
    while pos < N and values[pos] == NoData:
        pos = pos + 1
    if pos < N:
        max = values[pos]
        maxpos = 0
        Nact = 0
        for i in range(pos, N):
            if values[i] != NoData and not np.isnan(values[i]):
                if values[i] > max:
                    max = values[i]
                    maxpos = i
                Nact = Nact + 1
        if Nact == 0:
            max = NoData
    else:
        max = NoData
        maxpos = NoData
        Nact = 0
    return (max, Nact, maxpos)


def get_count_over_threshold(values, threshold, N="", NoData=NoDataVal):
    """This function determines the number of values that are equal to
    or exceed the given threshold.  Values equal to NoData are not
    included in the count and the number of valid values is returned
    along with the over threshold count."""
    if not N:
        N = len(values)
    count = 0
    Nact = 0
    for i in range(N):
        if values[i] != NoData and not np.isnan(values[i]):
            if values[i] >= threshold:
                count = count + 1
            Nact = Nact + 1
    if Nact == 0:
        count = NoData
    return (count, Nact)


def get_quantile(values, quantile, N="", NoData=NoDataVal):
    """This function selects the numeric value representing the
    requested quantile (0-1) from the original data set.  Data values
    are sorted low to high then the Weibull plotting function
    is used to determine the quantile value.  Values equal to
    NoData are not included in the process and the number of valid
    values is returned along with the quantile value."""
    if not N:
        N = len(values)
    tmpvalues = values
    tmpvalues.sort()
    while len(tmpvalues) > 0 and tmpvalues[0] == NoData:
        del tmpvalues[0]
    Nact = len(values)
    if Nact == 0:
        return (NoData, 0)
    else:
        quantile = int(quantile * (Nact + 1) + 0.5)
        return (tmpvalues[quantile], Nact)


def get_running_average(values, Navg, N="", NoData=NoDataVal):
    """This function computes a running average of Navg items
    and reports the results as an array of length N.  The returned
    array is padded at the start and end with NoData so that the
    average values are centered."""
    if not N:
        N = len(values)
    Nsets = N - Navg + 1
    avgvalues = [NoData] * (N)
    Nact = 0
    for i in range(Nsets):
        avgvalues[i + Navg / 2] = get_mean(values[i : (i + Navg)])[0]
        if avgvalues[i + Navg / 2] != NoData:
            Nact = Nact + 1
    if Nact == 0:
        avgvalues = NoData
    return (avgvalues, Nact)


def get_running_slope(values, Nslope=2, N="", NoData=NoDataVal):
    """This function computes running slopes between values at 0 and Nslope
    and reports the results as an array of length N.  The returned
    array is padded at the end with NoData so that the
    average values are centered."""
    if not N:
        N = len(values)
    slopes = [NoData] * (N)
    Nact = 0
    Nsets = N - Nslope
    for i in range(Nsets):
        if values[i] == NoData or values[i + Nslope] == NoData:
            slopes[i + Nslope - 1] = NoData
        else:
            slopes[i + Nslope - 1] = (values[i + Nslope] - values[i]) / float(Nslope)
        if slopes[i + Nslope - 1] != NoData:
            Nact = Nact + 1
    if Nact == 0:
        slopes = NoData
    return (slopes, Nact)


def get_bias(Avalues, Bvalues, N="", NoData=NoDataVal):
    """This function computes the bias between two arrays of length N.
    If using with streamflow, than Avalues should be the observed data
    and Bvalues are the simulated values.  The bias and the number of
    comparisons between actual data values (no NoData values) are
    returned."""
    if not N:
        N = len(Avalues)
    bias = 0
    Nact = 0
    for i in range(N):
        if (Avalues[i] != NoData and not np.isnan(Avalues[i])) and (
            Bvalues[i] != NoData and not np.isnan(Avalues[i])
        ):
            bias = bias + (Avalues[i] - Bvalues[i])
            Nact = Nact + 1
    if Nact == 0:
        bias = NoData
    else:
        bias = bias / Nact
    return (bias, Nact)


def get_root_mean_square(Avalues, Bvalues, N="", NoData=NoDataVal):
    """This function computes the root mean square error between two
    arrays of length N.  If using with streamflow, than Avalues should
    be the observed data and Bvalues are the simulated values.  The
    root mean squared error and the number of comparisons between actual
    data values (no NoData values) are returned."""
    if not N:
        N = len(Avalues)
    rmse = 0
    Nact = 0
    for i in range(N):
        if (Avalues[i] != NoData and not np.isnan(Avalues[i])) and (
            Bvalues[i] != NoData and not np.isnan(Avalues[i])
        ):
            rmse = rmse + (Avalues[i] - Bvalues[i]) * (Avalues[i] - Bvalues[i])
            Nact = Nact + 1
    if Nact == 0:
        rmse = NoData
    else:
        rmse = math.sqrt(rmse / Nact)
    return (rmse, Nact)


def get_mean_absolute_error(Avalues, Bvalues, N="", NoData=NoDataVal):
    """This function computes the mean absolute error between two
    arrays of length N.  If using with streamflow, than Avalues should
    be the observed data and Bvalues are the simulated values.  The
    mean absolute error and the number of comparisons between actual
    data values (no NoData values) are returned."""
    if not N:
        N = len(Avalues)
    abserr = 0
    Nact = 0
    for i in range(N):
        if (Avalues[i] != NoData and not np.isnan(Avalues[i])) and (
            Bvalues[i] != NoData and not np.isnan(Avalues[i])
        ):
            abserr = abserr + math.fabs(Avalues[i] - Bvalues[i])
            Nact = Nact + 1
    if Nact == 0:
        abserr = NoData
    else:
        abserr = math.sqrt(abserr / Nact)
    return (abserr, Nact)


def get_max_absolute_error(Avalues, Bvalues, N="", NoData=NoDataVal):
    """This function computes the maximum absolute error between two
    arrays of length N.  If using with streamflow, than Avalues should
    be the observed data and Bvalues are the simulated values.  The
    maximum absolute error and the number of comparisons between actual
    data values (no NoData values) are returned."""
    if not N:
        N = len(Avalues)
    absmax = []
    Nact = 0
    for i in range(N):
        if (Avalues[i] != NoData and not np.isnan(Avalues[i])) and (
            Bvalues[i] != NoData and not np.isnan(Avalues[i])
        ):
            absmax = absmax + [fabs(Avalues[i] - Bvalues[i])]
            Nact = Nact + 1
    if Nact == 0:
        absmax = NoData
    else:
        absmax = get_max(absmax)[0]
    return (absmax, Nact)


def get_nash_sutcliffe(Avalues, Bvalues, N="", NoData=NoDataVal):
    """This function computes the nash-sutcliffe R^2 between two
    arrays of length N.  If using with streamflow, than Avalues should
    be the observed data and Bvalues are the simulated values.  The
    nash-sutcliffe R^2 and the number of comparisons between actual
    data values (no NoData values) are returned."""
    if not N:
        N = len(Avalues)
    num = 0
    denom = 0
    Nact = 0
    mean_val = get_mean(Avalues, NoData=NoData)[0]
    for i in range(N):
        if (Avalues[i] != NoData and not np.isnan(Avalues[i])) and (
            Bvalues[i] != NoData and not np.isnan(Avalues[i])
        ):
            num = num + (Avalues[i] - Bvalues[i]) * (Avalues[i] - Bvalues[i])
            denom = denom + (Avalues[i] - mean_val) * (Avalues[i] - mean_val)
            Nact = Nact + 1
    if Nact == 0 or denom == 0:
        NS = NoData
    else:
        NS = 1 - (num / Nact) / (denom / Nact)
    return (NS, Nact)


def get_peak_diff(Avalues, Bvalues, N="", NoData=NoDataVal):
    """This function computes the peak differences between two
    arrays of length N.  If using with streamflow, than Avalues should
    be the observed data and Bvalues are the simulated values.  The
    peak differences and the number of comparisons between actual
    data values (no NoData values) are returned."""
    max_A, Aact, Rec = get_max(Avalues, NoData=NoData)
    max_B, Bact, Rec = get_max(Bvalues, NoData=NoData)
    if max_A == NoData or max_B == NoData:
        Nact = 0
        return (NoData, Nact)
    else:
        Nact = (Aact + Bact) / 2
        return (math.fabs(max_A - max_B), Nact)


def get_number_of_sign_changes(Avalues, Bvalues, N="", NoData=NoDataVal):
    """This function computes the number of sign changes between two
    arrays of length N.  If using with streamflow, than Avalues should
    be the observed data and Bvalues are the simulated values.  The
    number of sign changes (times -1) and the number of comparisons
    between actual data values (no NoData values) are returned."""
    if not N:
        N = len(Avalues)
    if Avalues[0] != Bvalues[0]:
        sign = (Avalues[0] - Bvalues[0]) / math.fabs(Avalues[0] - Bvalues[0])
    else:
        sign = 1
    NSC = 0
    Nact = 0
    for i in range(N):
        if (Avalues[i] != NoData and not np.isnan(Avalues[i])) and (
            Bvalues[i] != NoData and not np.isnan(Avalues[i])
        ):
            if Avalues[i] != Bvalues[i]:
                curr_sign = (Avalues[i] - Bvalues[i]) / math.fabs(Avalues[i] - Bvalues[i])
            else:
                curr_sign = sign
            if curr_sign != sign:
                NSC = NSC + 1
                sign = curr_sign
            Nact = Nact + 1
    if Nact == 0:
        NSC = NoData
    return (-NSC, Nact)


def get_peak_threshold_diff(Avalues, Bvalues, threshold, N="", NoData=NoDataVal):
    """This functions computes the difference in the number of peaks
    over threshold between two arrays of N values.  If using with
    streamflow, than Avalues should be the observed data and Bvalues
    are the simulated values.  The absolute difference in peaks over
    threshold and the number of comparisons between actual
    data values (no NoData values) are returned."""
    if not N:
        N = len(Avalues)
    Apeaks, Aact = filter_threshold(Avalues, threshold, NoData=NoData)
    Bpeaks, Bact = filter_threshold(Bvalues, threshold, NoData=NoData)
    if Aact == 0 or Bact == 0:
        return (NoData, 0)
    else:
        return (math.fabs(Apeaks - Bpeaks), (Aact + Bact) / 2)


def get_covariance(Avalues, Bvalues, N="", NoData=NoDataVal):
    """This functions computes the covariance between two arrays of
    N values.  The covariance between the two series and the number
    of comparisons between actual data values (no NoData values) are
    returned."""
    if not N:
        N = len(Avalues)
    Asum = 0
    Bsum = 0
    ABsum = 0
    Nact = 0
    for idx in range(N):
        if Avalues[idx] != NoData and Bvalues[idx] != NoData:
            Asum = Asum + Avalues[idx]
            Bsum = Bsum + Bvalues[idx]
            ABsum = ABsum + Avalues[idx] * Bvalues[idx]
            Nact = Nact + 1
    if Nact < 2:
        return (NoData, Nact)
    COV = (ABsum - (Asum * Bsum) / Nact) / (Nact - 1)
    return (COV, Nact)


def get_correlation(
    Avalues, Bvalues, Amean="", Bmean="", Astdev="", Bstdev="", N="", NoData=NoDataVal
):
    """This functions computes the correlation coefficient between two
    arrays of N values.  The correlation coeffificent between the two
    series and the number of comparisons between actual data values
    (no NoData values) are returned."""
    if not N:
        N = len(Avalues)
    Aact = Bact = N
    COV, Nact = get_covariance(Avalues, Bvalues, N, NoData)
    if not Amean:
        Amean = get_mean(Avalues, N, NoData)[0]
    if not Astdev:
        Astdev, Aact = get_stdev(Avalues, N, Amean, NoData)
    if not Bmean:
        Bmean = get_mean(Bvalues, N, NoData)[0]
    if not Bstdev:
        Bstdev, Bact = get_stdev(Bvalues, N, Bmean, NoData)
    if Aact == 0 or Bact == 0:
        return (NoData, Nact)
    if Astdev < SmallValue or Bstdev < SmallValue:
        return (NoData, Nact)
    return (COV / Astdev / Bstdev, Nact)


def get_cross_correlation(
    Avalues,
    Bvalues,
    lag,
    Amean="",
    Bmean="",
    Astdev="",
    Bstdev="",
    N="",
    NoData=NoDataVal,
):
    """This functions computes the cross-correlation between two data
    series using the given lagto shift values in Bvalues forward (negative)
    and backwards (positive) versus in Avalues.  The cross-correlation
    coefficient and the number of comparisons between actual
    data values (no NoData values) are returned."""
    if not N:
        N = len(Avalues)
    if abs(lag) >= N:
        return (NoData, 0)
    if lag > 0:
        Astart = lag
        Aend = len(Avalues)
        Bstart = 0
        Bend = len(Bvalues) - lag
    elif lag < 0:
        Astart = 0
        Aend = len(Avalues) + lag
        Bstart = -lag
        Bend = len(Bvalues)
    else:
        Astart = 0
        Aend = len(Avalues)
        Bstart = 0
        Bend = len(Bvalues)
    N = len(Avalues[Astart:Aend])
    if not Amean:
        Amean, Aact = get_mean(Avalues[Astart:Aend], N, NoData)
    if not Astdev:
        Astdev, Aact = get_stdev(Avalues[Astart:Aend], N, Amean, NoData)
    if not Bmean:
        Bmean, Bact = get_mean(Bvalues[Bstart:Bend], N, NoData)
    if not Bstdev:
        Bstdev, Bact = get_stdev(Bvalues[Bstart:Bend], N, Bmean, NoData)
    if Aact == 0 or Bact == 0:
        return (NoData, 0)
    COR, Nact = get_correlation(
        Avalues[Astart:Aend],
        Bvalues[Bstart:Bend],
        Amean,
        Bmean,
        Astdev,
        Bstdev,
        N,
        NoData,
    )

    return (COR, Nact)


def filter_threshold(values, threshold, FILTER="ABOVE", N="", NoData=NoDataVal):
    """This function counts the number of peaks above a threshold in an
    array of length N.  The number of peaks and the number of actual
    data values (no NoData values) are returned."""
    if not N:
        N = len(values)
    Npeaks = 0
    Nact = 0
    for i in range(N):
        if values[i] != NoData and not np.isnan(values[i]):
            if FILTER == "ABOVE":
                if values[i] >= threshold:
                    Npeaks = Npeaks + 1
            else:
                if values[i] <= threshold:
                    Npeaks = Npeaks + 1
            Nact = Nact + 1

    return (Npeaks, Nact)


def get_days(values, N="", NoData=NoDataVal, Thres=1.0):
    """This function computes the number of time steps an array of
       values is above a threshold (i.e. the number of snow covered days."""
    if not N:
        N = len(values)
    days = 0
    Nact = 0
    for i in range(N):
        if values[i] >= Thres and values[i] != NoData:
            days = days + 1
            Nact = Nact + 1
    if Nact == 0:
        days = NoData
    return (days, Nact)


def get_last_day(values, N="", NoData=NoDataVal, Thres=1.0):
    """This function determines the index of the last time the
       array of values exceeds the threshold (i.e. last day of snow)."""
    if not N:
        N = len(values)
    last_day = 0
    Nact = 0
    for i in range(N - 1, -1, -1):
        if values[i] != NoData and values[i] > Thres:
            last_day = i
            break
    return last_day


def get_first_day(values, N="", NoData=NoDataVal, Thres=1.0):
    """This function determines the index of the first time the
       array of values exceeds the threshold (i.e. last day of snow)."""
    if not N:
        N = len(values)
    first_day = 0
    Nact = 0
    for i in range(N):
        if values[i] != NoData and values[i] > Thres:
            first_day = i
            break
    return first_day


def get_box_plot_parameters(values, N="", NoData=NoDataVal):
    """This function estimates the values required to create a box
    and whiskers plot.  Values equal to NoData are not included in
    the process and the number of valid values is returned along
    with the quantile value.  Modified November 8, 2007 to output
    values in the order required by GMT (median, min, 25%, 75%, max)."""
    tmpvalues = values
    # sort data and strip no data values
    tmpvalues.sort()
    try:
        while tmpvalues[0] == NoData and len(tmpvalues) > 0:
            del tmpvalues[0]
    except IndexError as errstr:
        return (
            NoData,
            NoData,
            NoData,
            NoData,
            NoData,
            NoData,
            NoData,
            [NoData],
            [NoData],
        )

    Nact = len(values)
    # get median
    SubSetVals = [0] * 2
    if Nact > 0:
        if Nact % 2 == 0:
            median = (tmpvalues[int(Nact / 2)] + tmpvalues[int(Nact / 2)]) / 2.0
            SubSetVals[0] = tmpvalues[: int(Nact / 2) + 1]
            SubSetVals[1] = tmpvalues[int(Nact / 2) + 1 :]
        else:
            median = tmpvalues[int(Nact / 2)]
            SubSetVals[0] = tmpvalues[: int(Nact / 2)]
            SubSetVals[1] = tmpvalues[int(Nact / 2) + 1 :]
    else:
        return (
            NoData,
            NoData,
            NoData,
            NoData,
            NoData,
            NoData,
            NoData,
            [NoData],
            [NoData],
        )
    N = Nact

    # get mean
    mean = get_mean(tmpvalues, NoData=NoData)[0]

    # get quartiles
    quartiles = [0] * 2
    for subset in range(2):
        Nact = len(SubSetVals[subset])
        if Nact > 0:
            if Nact % 2 == 0:
                quartiles[subset] = (
                    SubSetVals[subset][int(Nact / 2)]
                    + SubSetVals[subset][int(Nact / 2)]
                ) / 2.0
            else:
                quartiles[subset] = SubSetVals[subset][int(Nact / 2)]
        else:
            return (NoData, NoData, NoData, NoData, NoData, NoData, [NoData], [NoData])

    # compute interquartile range
    IQR = quartiles[1] - quartiles[0]

    # find outliers
    ExtremeOutliers = []
    MildOutliers = []
    for idx in range(Nact - 1, -1, -1):
        if (
            tmpvalues[idx] < quartiles[0] - 3.0 * IQR
            or tmpvalues[idx] > quartiles[1] + 3.0 * IQR
        ):
            ExtremeOutliers = ExtremeOutliers + [tmpvalues[idx]]
            del tmpvalues[idx]
        elif (
            tmpvalues[idx] < quartiles[0] - 1.5 * IQR
            or tmpvalues[idx] > quartiles[1] + 1.5 * IQR
        ):
            MildOutliers = MildOutliers + [tmpvalues[idx]]
            del tmpvalues[idx]

    # find minimum and maximum
    MinVal = tmpvalues[0]
    MaxVal = tmpvalues[-1]

    return (
        median,
        MinVal,
        quartiles[0],
        quartiles[1],
        MaxVal,
        mean,
        N,
        MildOutliers,
        ExtremeOutliers,
    )
