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
import wflow.pcrut as pcrut
import sys
import getopt
import logging


def usage(*args):
    sys.stdout = sys.stderr
    """Way"""
    for msg in args:
        print(msg)
    print(__doc__)
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
        opts, args = getopt.getopt(argv, "c:l:",['version'])
    except getopt.error as msg:
        usage(msg)

    loglevel = logging.WARN
    for o, a in opts:
        if o == "-c":
            configfile = a
        if o == "-l":
            exec("loglevel = logging." + a)
        if o == "--version":
            import wflow
            print("wflow version: ", wflow.__version__)
            sys.exit(0)

    combilogger = pcrut.setlogger(
        "bmi2runner.log", "bmi2runner_logging", thelevel=loglevel
    )

    # Construct object and initilize the models
    combilogger.info("Starting combined bmi object")
    bmiobj = wfbmi.wflowbmi_csdms()

    bmiobj.initialize_config(configfile, loglevel=loglevel)
    bmiobj.initialize_model()

    # Get and set start and end times
    start = bmiobj.get_start_time()
    end = bmiobj.get_end_time()
    bmiobj.set_start_time(start)
    bmiobj.set_end_time(end)

    # Update models (if necessary) to start time
    bmiobj.update_to_start_time(start)

    # Number of steps to run models
    ts = bmiobj.get_time_step()
    steps = int((end - start) / ts + 1)

    cts = bmiobj.currenttimestep
    # Loop over the time duration
    while cts < steps:
        combilogger.info("time is: " + str(bmiobj.get_current_time()))
        bmiobj.update()
        cts = bmiobj.currenttimestep

    bmiobj.finalize()
    combilogger.info("Finishing run")


if __name__ == "__main__":
    main()
