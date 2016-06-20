from Scripts.UI_Examples.Shortcuts import *

def OpenDoc():
	import webbrowser
	webbrowser.open("http://wflow.readthedocs.io/en/latest/")


name = "Wflow-Tools"
tabName = "Wflow-Tools"
groupName = "Wflow"

CreateShortcutButton(name,groupName,tabName, lambda: OpenDoc(), None)

RemoveShortcut(name,groupName,tabName)
RemoveShortcutsTab(tabName)