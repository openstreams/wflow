"""

bmi2runner - runs multiple (2) linked bmi models

Usage:

        bmi2runner -c configfile

"""

import wflow.wflow_bmi_combined as wfbmi
import getopt
import sys
import logging
import wflow.wflow_bmi_combined as bmi
import wflow.pcrut as pcrut
import wflow.wflow_adapt as wfa



def usage(*args):
    sys.stdout = sys.stderr
    """Way"""
    for msg in args: print msg
    print __doc__
    sys.exit(0)


def main(argv=None):
    """
    Perform command line execution of the model.
    """

    configfile = "bmirunner.ini"
    LogFileName = "bmirunner.log"
    loglevel = logging.DEBUG

    if argv is None:
        argv = sys.argv[1:]
        if len(argv) == 0:
            usage()
            return
    ########################################################################
    ## Process command-line options                                        #
    ########################################################################
    try:
        opts, args = getopt.getopt(argv, 'c:')
    except getopt.error, msg:
        usage(msg)

    for o, a in opts:
        if o == '-c': configfile = a
        if o == '-l': exec "loglevel = logging." + a


    logje = pcrut.setlogger('bmi2runner','bmi2runner',thelevel=loglevel)
    bmiobj = bmi.wflowbmi_csdms()
    bmiobj.initialize(configfile)
    start = bmiobj.get_start_time()
    end = bmiobj.get_end_time()
    ts = bmiobj.get_time_step()
    curtime = bmiobj.get_current_time()

    while curtime < end:
        bmiobj.update_until(curtime + ts)
        curtime = bmiobj.get_current_time()

    bmiobj.finalize()



if __name__ == "__main__":
    main()




