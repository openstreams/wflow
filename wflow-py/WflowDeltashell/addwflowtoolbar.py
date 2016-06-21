from Scripts.Shortcuts import *

def OpenDoc(url):
	import webbrowser
	webbrowser.open(url)


name = "Web Documentation"
tabName = "Wflow-Tools"
groupName = "Wflow"

CreateShortcutButton(name,groupName,tabName, lambda: OpenDoc("http://wflow.readthedocs.io/en/latest/"), None)

name = "Github"
tabName = "Wflow-Tools"
groupName = "Wflow"

CreateShortcutButton(name,groupName,tabName, lambda: OpenDoc("http://github.com/openstreams/wflow"), None)

#RemoveShortcut(name,groupName,tabName)
#RemoveShortcutsTab(tabName)