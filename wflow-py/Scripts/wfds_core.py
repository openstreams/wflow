"""
Usage:
  call_functionality <engine> <configfile> [-o <outputvar>...] [-g <globalvar>...] [--interval <interval>] [--disable-logger] [--pause] [--mpi <method>] [--track <server>] [--port <port>] [--bmi-class]
  call_functionality -h | --help

Positional arguments:
  engine model  engine
  configfile    model configuration file

Optional arguments:
  -h, --help               show this help message and exit
  --interval <interval>    publish results every <interval> timesteps
  -o <outputvar>           output variables, will be broadcasted each <interval> timestep
  -g <globalvar>           global variables, will be send if requested
  --disable-logger         do not inject logger into the BMI library
  --pause                  start in paused mode, send update messages to progress
  --mpi <method>           communicate with mpi nodes using one of the methods: root (communicate with rank 0), all (one socket per rank)
  --port <port>            "random" or integer base port, port is computed as req/rep = port + rank*3 + 0, push/pull = port + rank*3 + 1, pub/sub = port + rank*3 + 2 [default: random]
  --track <tracker>        server to subscribe to for tracking
  --bmi-class              when used - the engine is assumed to be the full name of a Python class that implements bmi [default: bmi.wrapper.BMIWrapper]

"""
import bmi
import mmi.runner


import sys
import docopt

# from mmi.runner import runner


def call_functionality(arguments):
    print "hello"
    print arguments
    mmi.runner.runner(arguments)


arguments = docopt.docopt(__doc__)
call_functionality(arguments)
# call_functionality(sys.argv[1:])
