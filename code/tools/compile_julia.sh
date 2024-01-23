#!/bin/bash

cd "$(dirname $0)/../"


if [ "$1" = "" ] || [ "$1" = "model" ] || [ "$1" = "expe_template" ]
then
    enddir="nh"
else
    enddir="nh_$1"
fi

echo "$enddir"

mkdir "exe/$enddir" 2> "/dev/null"
mkdir "exe/$enddir/inc" 2> "/dev/null"
mkdir "exe/$enddir/src" 2> "/dev/null"

cp -fp model/inc/*.jl "exe/$enddir/inc/"
cp -fp model/src/*.jl "exe/$enddir/src/"


if [ "$1" != "" ]
then
    cp -fp $1/inc/*.jl "exe/$enddir/inc/" 2> /dev/null
    cp -fp $1/src/*.jl "exe/$enddir/src/" 2> /dev/null
fi