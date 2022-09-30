#!/bin/sh

export LD_LIBRARY_PATH=./lib::CfgManager/lib/:$LD_LIBRARY_PATH
export DYLD_LIBRARY_PATH=./lib:CfgManager/lib/:$DYLD_LIBRARY_PATH
export ROOT_INCLUDE_PATH=./interface:CfgManager/interface/:$ROOT_INCLUDE_PATH
