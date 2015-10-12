import os
import sys

def set_path():
	print os.environ['WS_HOME']
	global WS_LOC_SHELL 
	WS_LOC_SHEL= os.environ['WS_HOME']+'/DISTIQUE/src/shell'
	global WS_LOC_FASTME
	WS_LOC_FASTME = os.environ['WS_HOME']+'fastme-2.1.4/src'
	return
