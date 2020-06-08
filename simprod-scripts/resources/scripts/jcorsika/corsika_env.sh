#!/bin/env bash
MY_PATH="`dirname "${BASH_SOURCE[0]}"`"
dir="`( cd \"$MY_PATH\" && pwd )`"
echo export FLUPRO=${dir}/fluka
echo export LD_LIBRARY_PATH=${dir}/lib:$LD_LIBRARY_PATH

