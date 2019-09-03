"""
Converts tss files to pi_xml

usage:

    tss2xml.py -I tssfile [-s timestep in secs][-X xmlfilename][-S startdate]


-s timestep of tss file in seconds
-S startdatetime ('%Y-%m-%d %H:%M:%S')
-I input tssfile
-X xmlfile (defaylt tssfile + .xml)
"""

import wflow.wflow_adapt as wf
import wflow.pcrut as pcrut
import getopt
import sys


def usage(*args):
    sys.stdout = sys.stderr
    for msg in args:
        print(msg)
    print(__doc__)
    sys.exit(0)


def main(argv=None):
    """
    Perform command line execution of the script.
    """
    tssfile = "input.tss"
    xmlfile = tssfile + ".xml"
    timestepsecs = 86400
    parameter = tssfile.split(".")[0]
    startdatestr = "1970-01-01 00:00:00"
    startdate = wf.datetime.strptime(startdatestr, "%Y-%m-%d %H:%M:%S")

    if argv is None:
        argv = sys.argv[1:]
        if len(argv) == 0:
            usage()
            return

    try:
        opts, args = getopt.getopt(argv, "X:I:S",['version'])
    except getopt.error as msg:
        pcrut.usage(msg)

    for o, a in opts:
        if o == "-X":
            xmlfile = a
        if o == "-I":
            tssfile = a
            xmlfile = tssfile + ".xml"
        if o == "-s":
            timestepsecs = a
        if o == "-s":
            timestepsecs = a
        if o == "--version":
            import wflow
            print("wflow version: ", wflow.__version__)
            sys.exit(0)

    wf.tss_topixml(tssfile, xmlfile, "wflow", parameter, startdate, timestepsecs)


if __name__ == "__main__":
    main()
