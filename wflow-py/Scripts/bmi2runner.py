"""

bmi2runner - runs multiple linked bmi models

Usage:

        bmi2runner -c configfile

Example ini file:

::

    [models]
    wflow_sbm=wflow_sbm@wflow_sbm_comb.ini
    wflow_routing=wflow_routing@wflow_routing_comb.ini

    [exchanges]
    # From_model/var -> To_model/var
    wflow_sbm@InwaterMM=wflow_routing@IW

"""

import wflow.wflow_bmi_combined as wfbmi
import getopt
import sys
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
        opts, args = getopt.getopt(argv, 'c:')
    except getopt.error, msg:
        usage(msg)

    for o, a in opts:
        if o == '-c': configfile = a


    # Construct object and initilize the models
    bmiobj = bmi.wflowbmi_csdms()
    bmiobj.initialize_config(configfile)
    bmiobj.initialize_model()
    # Get time for the loop
    start = bmiobj.get_start_time()
    end = bmiobj.get_end_time()
    ts = bmiobj.get_time_step()
    curtime = bmiobj.get_current_time()
    # Loop over the time duration
    while curtime < end:
        bmiobj.update_until(curtime + ts)
        curtime = bmiobj.get_current_time()

    bmiobj.finalize()



if __name__ == "__main__":
    main()




