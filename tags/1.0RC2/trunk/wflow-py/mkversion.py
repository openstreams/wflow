import os

thisversion = "$Id: mkversion.py 550 2012-11-27 19:17:46Z schelle $"
vers=os.popen('svnversion -n').read()
manualversion = "1.0RC2"
a = open("_version.py","w")

a.write("VERSION=\"" + vers + " " + thisversion + "\"\n")
a.write("MVERSION=\"" + manualversion +"\"\n")
