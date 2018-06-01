from WflowDeltashell.Shortcuts import *
from WflowDeltashell.plotcsv import *


def OpenDoc(url):
    import Libraries.StandardFunctions as sf
    from DelftTools.Utils import Url

    murl = Url(url, url)
    sf.OpenView(murl)


def notimplemented():
	print("Not implemented yet...")

name = "Web Documentation"
tabName = "Wflow-Tools"
groupName = "Internet"

CreateShortcutButton(
    name,
    groupName,
    tabName,
    lambda: OpenDoc("http://wflow.readthedocs.io/en/latest/"),
    None,
)

name = "Github"
tabName = "Wflow-Tools"
groupName = "Internet"

CreateShortcutButton(
    name,
    groupName,
    tabName,
    lambda: OpenDoc("http://github.com/openstreams/wflow"),
    None,
)

name = "Plotcsv"
tabName = "Wflow-Tools"
groupName = "Plots"

CreateShortcutButton(name, groupName, tabName, lambda: plotit(getcsvname()), None)

name = "Netcdf Input"
tabName = "Wflow-Tools"
groupName = "Conversion"

CreateShortcutButton(name, groupName, tabName, lambda: notimplemented(), None)

# RemoveShortcut(name,groupName,tabName)
# RemoveShortcutsTab(tabName)
