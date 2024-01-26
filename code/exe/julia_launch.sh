#!/bin/bash

if [ $# -lt 2 ]
then
    echo "Usage: ./julia_launch [experiment_folder] [namelist] [julia_opts]"
    echo " - experiment_folder: folder name in exe/ to start the experiment"
    echo " - namelist: path to the namelist YAML file"
    echo " - julia_opts: additional julia options  (OPTIONAL)"
    exit 1
fi

namelist="$PWD/$2"

cd "$(dirname $0)/$1"
echo "$namelist"

# JULIA options

#number of thread julia may use (may be set to auto)
threads=16

julia_options="-t$threads ${@:3}"


julia $julia_options -- src/main.jl "$namelist" "./PSOM_LIB.so"
