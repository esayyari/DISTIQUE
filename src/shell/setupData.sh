#!/bin/bash

DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

source $DIR/setup.sh

export WS_LOC_DT=$WS_HOME/data
export PATH=$PATH:$WS_LOC_DT
