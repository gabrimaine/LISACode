#!/bin/bash

################################################################
##### Script for setting a working directory with LISACode #####
################################################################


if [[ "$1" = "--help" ]]
	then DispHelp='1'
else if [[ "$1" = "-h" ]]
		then DispHelp='1'
		else DispHelp='0'
	fi
fi
if [[ "$DispHelp" = "1" ]]
then
	echo "Script for setting a working directory with LISACode :"
	echo "------------------------------------------------------"
	echo "   - argument 1 : Path of working directory [required]"
	echo "   - argument 2 : Some id number. example : 4 for Cfg4 "
	exit
fi

if test $# -lt 1
then
	echo "Error : 1 argument (path of working directory is required) required"
	exit
fi


## Copy the executable
mkdir $1/bin$2
cp Main/Exe/LISACode2 $1/bin$2/.
cp Main/Exe/LC2Sensitivty $1/bin$2/.
cp Main/Exe/LC2PSD $1/bin$2/.
cp Main/Exe/LC2MeanFiles $1/bin$2/.
cp Main/Exe/LC2SNR $1/bin$2/.
cp Main/Exe/LC2FIM $1/bin$2/.
cp Main/Exe/LC2SNRFIMMulti $1/bin$2/.
cp Main/Exe/LC2MCMC $1/bin$2/.
cp Noise/Test/LC2TestNoise $1/bin$2/.
cp GW/Test/LC2GWbary $1/bin$2/.
cp Orbits/Test/LC2Orbits $1/bin$2/.
cp ArmResponse/Test/LC2GWArmResp $1/bin$2/.
chmod u+x Main/Exe/*.py
cp Main/Exe/makeTDI-lisacode2.py $1/bin$2/.
#for f in $(find $1/bin$2/.); do sed "s#WORKDIR/Cfg#$1/Cfg$2#g" $f > LCtmp ; mv -f LCtmp $f ; sed "s#WORKDIR/bin#$1/bin$2#g" $f > LCtmp ; mv -f LCtmp $f ; done
for f in $1/bin$2/* ; do sed "s#WORKDIR/Cfg#$1/Cfg$2#g" $f > LCtmp ; mv -f LCtmp $f ; sed "s#WORKDIR/bin#$1/bin$2#g" $f > LCtmp ; mv -f LCtmp $f ; done
chmod u+x $1/bin$2/*

## Copy configuration files
echo "Copy ConfigFiles directory to $1/Cfg$2 ..."
cp -rf ConfigFiles $1/Cfg$2

## Remove all .svn directory (until level 5)
echo "Remove all .svn directory (until level 5) ..."
rm -rf $1/Cfg$2/.svn
rm -rf $1/Cfg$2/*/.svn
rm -rf $1/Cfg$2/*/*/.svn
rm -rf $1/Cfg$2/*/*/*/.svn
rm -rf $1/Cfg$2/*/*/*/*/.svn
rm -rf $1/Cfg$2/*/*/*/*/*/.svn


## Change WORKDIR by effective working directory
echo "Change WORKDIR by effective working directory $2 in all files ..."
#for f in $(find $1/Cfg$2/.); do sed "s#WORKDIR/Cfg#$1/Cfg$2#g" $f > LCtmp ; mv -f LCtmp $f ; done
for f in $1/Cfg$2/* ; do sed "s#WORKDIR/Cfg#$1/Cfg$2#g" $f > LCtmp ; mv -f LCtmp $f ; done
for f in $1/Cfg$2/*/* ; do sed "s#WORKDIR/Cfg#$1/Cfg$2#g" $f > LCtmp ; mv -f LCtmp $f ; done
for f in $1/Cfg$2/*/*/* ; do sed "s#WORKDIR/Cfg#$1/Cfg$2#g" $f > LCtmp ; mv -f LCtmp $f ; done
for f in $1/Cfg$2/*/*/*/* ; do sed "s#WORKDIR/Cfg#$1/Cfg$2#g" $f > LCtmp ; mv -f LCtmp $f ; done
rm $1/Cfg$2/LCtmp
rm $1/Cfg$2/*/LCtmp
rm $1/Cfg$2/*/*/LCtmp
rm $1/Cfg$2/*/*/*/LCtmp