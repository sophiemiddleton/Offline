#!/bin/bash

# default geometry definition: Mu2eG4/geom/geom_2019_PhaseI_hayman_v2.txt
#  previous one              : JobConfig/common/geom_baseline.txt

# oldgeom="Mu2eG4/geom/geom_2019_PhaseI_hayman_v2.txt"
oldgeom="Mu2eG4/geom/geom_common_su2020.txt"
newgeom="JobConfig/common/geom_baseline.txt"

for f in `find JobConfig -name \*.fcl`; do
    x=`grep $oldgeom $f`
    if [ ".$x" != "." ] ; then
	echo --- $f
	f1=$f.tmp
	cat $f | sed "s#$oldgeom#$newgeom#" >| $f1
	mv $f1 $f
    fi
done
