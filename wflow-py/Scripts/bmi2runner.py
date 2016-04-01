"""

bmi2runner - runs multiple linked bmi models

Usage:

        bmi2runner -c configfile -l loglevel

        -l: loglevel (most be one of DEBUG, WARNING, ERROR)

Example ini file:

::

    [models]
    wflow_sbm=wflow_sbm/wflow_sbm_comb.ini
    wflow_routing=wflow_routing/wflow_routing_comb.ini

    [exchanges]
    # From_model/var -> To_model/var
    wflow_sbm@InwaterMM=wflow_routing@IW

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


    if argv is None:
        argv = sys.argv[1:]
        if len(argv) == 0:
            usage()
            return
    ########################################################################
    ## Process command-line options                                        #
    ########################################################################
    try:
        opts, args = getopt.getopt(argv, 'c:l:')
    except getopt.error, msg:
        usage(msg)

    loglevel=logging.WARN
    for o, a in opts:
        if o == '-c': configfile = a
        if o == '-l': exec "loglevel = logging." + a


    combilogger = pcrut.setlogger('bmi2runner.log','bmi2runner_logging',thelevel=loglevel)
    # Construct object and initilize the models
    combilogger.info('Starting combined bmi object')
    bmiobj = bmi.wflowbmi_csdms()
    bmiobj.initialize_config(configfile,loglevel=loglevel)
    bmiobj.initialize_model()
    start = bmiobj.get_start_time()
    end = bmiobj.get_end_time()
    bmiobj.set_start_time(start)
    bmiobj.set_end_time(end)
    # Get time for the loop

    ts = bmiobj.get_time_step()
    curtime = bmiobj.get_current_time()
    # Loop over the time duration

    while curtime < end:
        combilogger.info("time is: " + str(curtime))
        bmiobj.update_until(curtime + ts)
        curtime = bmiobj.get_current_time()

    bmiobj.finalize()
    combilogger.info('Finishing run')



if __name__ == "__main__":
    main()




