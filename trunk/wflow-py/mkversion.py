import os

thisversion = "$Id: mkversion.py 550 2012-11-27 19:17:46Z schelle $"
vers=os.popen('svnversion -n').read()

a = open("0version.py","w")

a.write("VERSION=\"" + vers + " " + thisversion + "\"")
