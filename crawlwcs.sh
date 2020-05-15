#!/bin/bash


base=/archive/engineering
cameras="ak01"
sites="lsc"
dates="2019????"

inputselection="*-g00.fits.fz"

NCPU=30

# see https://stackoverflow.com/questions/51256738/multiple-instances-of-python-running-simultaneously-limited-to-35
export OPENBLAS_NUM_THREADS=1

for site in $sites; do
 for camera in $cameras; do

  sitecameras=`find ${base}/${site}  -maxdepth 1 -type d -wholename "*/$camera"`

  for sitecamera in $sitecameras; do

   directories=`find "${sitecamera}" -maxdepth 1 -type d  -wholename "*/${dates}" `

   for day in $directories; do



     searchpath=${day}/raw/${inputselection}
     searchpath=`ls $searchpath | shuf -n 400 |  xargs  echo`
     sem -j $NCPU python SIPWCSOptimizer.py --minmatched 40 --refcat2 /AstroCatalogs/Atlas-refcat2/refcat2.db   --inputfiles $searchpath

   done

  done

 done
done

sem --wait

