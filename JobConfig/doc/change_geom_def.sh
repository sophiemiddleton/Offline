#!/bin/bash

# default geometry definition: Mu2eG4/geom/geom_2019_PhaseI_hayman_v2.txt
#  previous one              : JobConfig/common/geom_baseline.txt

oldgeom="JobConfig/common/geom_baseline.txt"
newgeom="Mu2eG4/geom/geom_2019_PhaseI_hayman_v2.txt"

for f in `find JobConfig -name \*.fcl`; do
    x=`grep JobConfig/common/geom_baseline.txt $f`
    if [ ".$x" != "." ] ; then
	echo --- $f
	f1=$f.tmp
	cat $f | sed "s#$oldgeom#$newgeom#" >| $f1
	mv $f1 $f
    fi
done
