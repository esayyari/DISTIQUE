#!/bin/bash
if [ -z $WS_HOME ]; then
  echo you need to have setup.sh under $WS_HOME/global/src/shell/setup.sh from main repository
  echo you need to set $WS_HOME or else nothing else works.
  echo set $WS_HOME to the parent directory where the 'global' direcotry resides
  exit 1
fi
source $WS_HOME/global/src/shell/setup.sh
export WS_LOC_HOME=$WS_HOME/DISTIQ
export WS_LOC_LIB=$WS_LOC_HOME/lib
export WS_LOC_BIN=$WS_LOC_HOME/bin
export WS_LOC_SH=$WS_LOC_HOME/src/shell
export WS_LOC_PUTIL=$WS_LOC_HOME/src/utils
export WS_LOC_TS=$WS_LOC_HOME/tests
export PATH=$PATH:$WS_LOC_SH:$WS_LOC_BIN:$WS_LOC_TS:$WS_LOC_PUTIL

