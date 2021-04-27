#!/usr/bin/bash

grid_id=$1

idir=/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/$grid_id

declare -a files=(
"/pnfs/mu2e/tape/phy-sim/dig/mu2e/cosm0s41b0/Mu2eII/art/69/fb/dig.mu2e.cosm0s41b0.Mu2eII.001002_00000000.art" 
"/pnfs/mu2e/tape/phy-sim/dig/mu2e/cosm0s41b0/Mu2eII/art/25/2a/dig.mu2e.cosm0s41b0.Mu2eII.001002_00000502.art" 
"/pnfs/mu2e/tape/phy-sim/dig/mu2e/cosm0s41b0/Mu2eII/art/1d/b5/dig.mu2e.cosm0s41b0.Mu2eII.001002_00001003.art" 
"/pnfs/mu2e/tape/phy-sim/dig/mu2e/cosm0s41b0/Mu2eII/art/7e/e1/dig.mu2e.cosm0s41b0.Mu2eII.001002_00001503.art" 
"/pnfs/mu2e/tape/phy-sim/dig/mu2e/cosm0s41b0/Mu2eII/art/98/12/dig.mu2e.cosm0s41b0.Mu2eII.001002_00002006.art" 
"/pnfs/mu2e/tape/phy-sim/dig/mu2e/cosm0s41b0/Mu2eII/art/48/79/dig.mu2e.cosm0s41b0.Mu2eII.001002_00002507.art" 
"/pnfs/mu2e/tape/phy-sim/dig/mu2e/cosm0s41b0/Mu2eII/art/c6/4a/dig.mu2e.cosm0s41b0.Mu2eII.001002_00003007.art" 
"/pnfs/mu2e/tape/phy-sim/dig/mu2e/cosm0s41b0/Mu2eII/art/a9/bf/dig.mu2e.cosm0s41b0.Mu2eII.001002_00003507.art" 
"/pnfs/mu2e/tape/phy-sim/dig/mu2e/cosm0s41b0/Mu2eII/art/07/5f/dig.mu2e.cosm0s41b0.Mu2eII.001002_00004007.art" 
"/pnfs/mu2e/tape/phy-sim/dig/mu2e/cosm0s41b0/Mu2eII/art/3a/f3/dig.mu2e.cosm0s41b0.Mu2eII.001002_00004508.art" 
"/pnfs/mu2e/tape/phy-sim/dig/mu2e/cosm0s41b0/Mu2eII/art/58/f9/dig.mu2e.cosm0s41b0.Mu2eII.001002_00005009.art" 
"/pnfs/mu2e/tape/phy-sim/dig/mu2e/cosm0s41b0/Mu2eII/art/7e/f5/dig.mu2e.cosm0s41b0.Mu2eII.001002_00005510.art" 
"/pnfs/mu2e/tape/phy-sim/dig/mu2e/cosm0s41b0/Mu2eII/art/b5/72/dig.mu2e.cosm0s41b0.Mu2eII.001002_00006011.art" 
"/pnfs/mu2e/tape/phy-sim/dig/mu2e/cosm0s41b0/Mu2eII/art/a9/ec/dig.mu2e.cosm0s41b0.Mu2eII.001002_00006511.art" 
"/pnfs/mu2e/tape/phy-sim/dig/mu2e/cosm0s41b0/Mu2eII/art/9a/11/dig.mu2e.cosm0s41b0.Mu2eII.001002_00007011.art" 
"/pnfs/mu2e/tape/phy-sim/dig/mu2e/cosm0s41b0/Mu2eII/art/93/89/dig.mu2e.cosm0s41b0.Mu2eII.001002_00007512.art" 
"/pnfs/mu2e/tape/phy-sim/dig/mu2e/cosm0s41b0/Mu2eII/art/eb/e6/dig.mu2e.cosm0s41b0.Mu2eII.001002_00008013.art" 
"/pnfs/mu2e/tape/phy-sim/dig/mu2e/cosm0s41b0/Mu2eII/art/bf/d2/dig.mu2e.cosm0s41b0.Mu2eII.001002_00008515.art" 
"/pnfs/mu2e/tape/phy-sim/dig/mu2e/cosm0s41b0/Mu2eII/art/73/19/dig.mu2e.cosm0s41b0.Mu2eII.001002_00009016.art" 
"/pnfs/mu2e/tape/phy-sim/dig/mu2e/cosm0s41b0/Mu2eII/art/cc/97/dig.mu2e.cosm0s41b0.Mu2eII.001002_00009518.art" 
"/pnfs/mu2e/tape/phy-sim/dig/mu2e/cosm0s41b0/Mu2eII/art/92/ff/dig.mu2e.cosm0s41b0.Mu2eII.001002_00010018.art" 
"/pnfs/mu2e/tape/phy-sim/dig/mu2e/cosm0s41b0/Mu2eII/art/fd/8d/dig.mu2e.cosm0s41b0.Mu2eII.001002_00010518.art" 
"/pnfs/mu2e/tape/phy-sim/dig/mu2e/cosm0s41b0/Mu2eII/art/41/6f/dig.mu2e.cosm0s41b0.Mu2eII.001002_00011019.art" 
"/pnfs/mu2e/tape/phy-sim/dig/mu2e/cosm0s41b0/Mu2eII/art/72/50/dig.mu2e.cosm0s41b0.Mu2eII.001002_00011519.art" 
"/pnfs/mu2e/tape/phy-sim/dig/mu2e/cosm0s41b0/Mu2eII/art/98/19/dig.mu2e.cosm0s41b0.Mu2eII.001002_00012019.art" 
"/pnfs/mu2e/tape/phy-sim/dig/mu2e/cosm0s41b0/Mu2eII/art/93/f4/dig.mu2e.cosm0s41b0.Mu2eII.001002_00012519.art" 
"/pnfs/mu2e/tape/phy-sim/dig/mu2e/cosm0s41b0/Mu2eII/art/73/11/dig.mu2e.cosm0s41b0.Mu2eII.001002_00013019.art" 
"/pnfs/mu2e/tape/phy-sim/dig/mu2e/cosm0s41b0/Mu2eII/art/0a/cf/dig.mu2e.cosm0s41b0.Mu2eII.001002_00013519.art" 
"/pnfs/mu2e/tape/phy-sim/dig/mu2e/cosm0s41b0/Mu2eII/art/cf/60/dig.mu2e.cosm0s41b0.Mu2eII.001002_00014019.art" 
"/pnfs/mu2e/tape/phy-sim/dig/mu2e/cosm0s41b0/Mu2eII/art/12/db/dig.mu2e.cosm0s41b0.Mu2eII.001002_00014519.art" 
"/pnfs/mu2e/tape/phy-sim/dig/mu2e/cosm0s41b0/Mu2eII/art/4c/6e/dig.mu2e.cosm0s41b0.Mu2eII.001002_00015019.art" 
"/pnfs/mu2e/tape/phy-sim/dig/mu2e/cosm0s41b0/Mu2eII/art/c7/7e/dig.mu2e.cosm0s41b0.Mu2eII.001002_00015519.art" 
"/pnfs/mu2e/tape/phy-sim/dig/mu2e/cosm0s41b0/Mu2eII/art/33/75/dig.mu2e.cosm0s41b0.Mu2eII.001002_00016020.art" 
"/pnfs/mu2e/tape/phy-sim/dig/mu2e/cosm0s41b0/Mu2eII/art/69/cd/dig.mu2e.cosm0s41b0.Mu2eII.001002_00016520.art" 
"/pnfs/mu2e/tape/phy-sim/dig/mu2e/cosm0s41b0/Mu2eII/art/fd/d4/dig.mu2e.cosm0s41b0.Mu2eII.001002_00017021.art" 
"/pnfs/mu2e/tape/phy-sim/dig/mu2e/cosm0s41b0/Mu2eII/art/56/bd/dig.mu2e.cosm0s41b0.Mu2eII.001002_00017521.art" 
"/pnfs/mu2e/tape/phy-sim/dig/mu2e/cosm0s41b0/Mu2eII/art/ef/e2/dig.mu2e.cosm0s41b0.Mu2eII.001002_00018021.art" 
"/pnfs/mu2e/tape/phy-sim/dig/mu2e/cosm0s41b0/Mu2eII/art/9c/9c/dig.mu2e.cosm0s41b0.Mu2eII.001002_00018523.art" 
"/pnfs/mu2e/tape/phy-sim/dig/mu2e/cosm0s41b0/Mu2eII/art/c4/a9/dig.mu2e.cosm0s41b0.Mu2eII.001002_00019026.art" 
"/pnfs/mu2e/tape/phy-sim/dig/mu2e/cosm0s41b0/Mu2eII/art/6e/b7/dig.mu2e.cosm0s41b0.Mu2eII.001002_00019528.art" 
"/pnfs/mu2e/tape/phy-sim/dig/mu2e/cosm0s41b0/Mu2eII/art/d3/78/dig.mu2e.cosm0s41b0.Mu2eII.001002_00020029.art" 
"/pnfs/mu2e/tape/phy-sim/dig/mu2e/cosm0s41b0/Mu2eII/art/0e/b5/dig.mu2e.cosm0s41b0.Mu2eII.001002_00020530.art" 
"/pnfs/mu2e/tape/phy-sim/dig/mu2e/cosm0s41b0/Mu2eII/art/5e/d7/dig.mu2e.cosm0s41b0.Mu2eII.001002_00021030.art" 
"/pnfs/mu2e/tape/phy-sim/dig/mu2e/cosm0s41b0/Mu2eII/art/f6/77/dig.mu2e.cosm0s41b0.Mu2eII.001002_00021530.art" 
"/pnfs/mu2e/tape/phy-sim/dig/mu2e/cosm0s41b0/Mu2eII/art/c4/61/dig.mu2e.cosm0s41b0.Mu2eII.001002_00022031.art" 
"/pnfs/mu2e/tape/phy-sim/dig/mu2e/cosm0s41b0/Mu2eII/art/c6/39/dig.mu2e.cosm0s41b0.Mu2eII.001002_00022532.art" 
"/pnfs/mu2e/tape/phy-sim/dig/mu2e/cosm0s41b0/Mu2eII/art/50/d2/dig.mu2e.cosm0s41b0.Mu2eII.001002_00023033.art" 
"/pnfs/mu2e/tape/phy-sim/dig/mu2e/cosm0s41b0/Mu2eII/art/a0/40/dig.mu2e.cosm0s41b0.Mu2eII.001002_00023533.art" 
"/pnfs/mu2e/tape/phy-sim/dig/mu2e/cosm0s41b0/Mu2eII/art/6d/97/dig.mu2e.cosm0s41b0.Mu2eII.001002_00024034.art" 
"/pnfs/mu2e/tape/phy-sim/dig/mu2e/cosm0s41b0/Mu2eII/art/c6/7e/dig.mu2e.cosm0s41b0.Mu2eII.001002_00024536.art" 
"/pnfs/mu2e/tape/phy-sim/dig/mu2e/cosm0s41b0/Mu2eII/art/6d/d1/dig.mu2e.cosm0s41b0.Mu2eII.001002_00025039.art" 
"/pnfs/mu2e/tape/phy-sim/dig/mu2e/cosm0s41b0/Mu2eII/art/7e/27/dig.mu2e.cosm0s41b0.Mu2eII.001002_00025540.art" 
"/pnfs/mu2e/tape/phy-sim/dig/mu2e/cosm0s41b0/Mu2eII/art/af/4f/dig.mu2e.cosm0s41b0.Mu2eII.001002_00026041.art" 
"/pnfs/mu2e/tape/phy-sim/dig/mu2e/cosm0s41b0/Mu2eII/art/b0/47/dig.mu2e.cosm0s41b0.Mu2eII.001002_00026543.art" 
"/pnfs/mu2e/tape/phy-sim/dig/mu2e/cosm0s41b0/Mu2eII/art/74/b9/dig.mu2e.cosm0s41b0.Mu2eII.001002_00027043.art" 
"/pnfs/mu2e/tape/phy-sim/dig/mu2e/cosm0s41b0/Mu2eII/art/ba/f3/dig.mu2e.cosm0s41b0.Mu2eII.001002_00027543.art" 
"/pnfs/mu2e/tape/phy-sim/dig/mu2e/cosm0s41b0/Mu2eII/art/43/28/dig.mu2e.cosm0s41b0.Mu2eII.001002_00028045.art" 
"/pnfs/mu2e/tape/phy-sim/dig/mu2e/cosm0s41b0/Mu2eII/art/72/39/dig.mu2e.cosm0s41b0.Mu2eII.001002_00028545.art" 
"/pnfs/mu2e/tape/phy-sim/dig/mu2e/cosm0s41b0/Mu2eII/art/9d/b1/dig.mu2e.cosm0s41b0.Mu2eII.001002_00029048.art" 
"/pnfs/mu2e/tape/phy-sim/dig/mu2e/cosm0s41b0/Mu2eII/art/1e/4e/dig.mu2e.cosm0s41b0.Mu2eII.001002_00029552.art" 
"/pnfs/mu2e/tape/phy-sim/dig/mu2e/cosm0s41b0/Mu2eII/art/79/ca/dig.mu2e.cosm0s41b0.Mu2eII.001002_00030052.art" 
"/pnfs/mu2e/tape/phy-sim/dig/mu2e/cosm0s41b0/Mu2eII/art/05/a8/dig.mu2e.cosm0s41b0.Mu2eII.001002_00030552.art" 
"/pnfs/mu2e/tape/phy-sim/dig/mu2e/cosm0s41b0/Mu2eII/art/ab/65/dig.mu2e.cosm0s41b0.Mu2eII.001002_00031052.art" 
"/pnfs/mu2e/tape/phy-sim/dig/mu2e/cosm0s41b0/Mu2eII/art/e0/72/dig.mu2e.cosm0s41b0.Mu2eII.001002_00031552.art" 
"/pnfs/mu2e/tape/phy-sim/dig/mu2e/cosm0s41b0/Mu2eII/art/b0/72/dig.mu2e.cosm0s41b0.Mu2eII.001002_00032052.art" 
"/pnfs/mu2e/tape/phy-sim/dig/mu2e/cosm0s41b0/Mu2eII/art/01/b3/dig.mu2e.cosm0s41b0.Mu2eII.001002_00032552.art" 
"/pnfs/mu2e/tape/phy-sim/dig/mu2e/cosm0s41b0/Mu2eII/art/36/40/dig.mu2e.cosm0s41b0.Mu2eII.001002_00033052.art" 
"/pnfs/mu2e/tape/phy-sim/dig/mu2e/cosm0s41b0/Mu2eII/art/cf/9f/dig.mu2e.cosm0s41b0.Mu2eII.001002_00033555.art" 
"/pnfs/mu2e/tape/phy-sim/dig/mu2e/cosm0s41b0/Mu2eII/art/b1/f0/dig.mu2e.cosm0s41b0.Mu2eII.001002_00034055.art" 
"/pnfs/mu2e/tape/phy-sim/dig/mu2e/cosm0s41b0/Mu2eII/art/b0/f8/dig.mu2e.cosm0s41b0.Mu2eII.001002_00034555.art" 
"/pnfs/mu2e/tape/phy-sim/dig/mu2e/cosm0s41b0/Mu2eII/art/b7/f6/dig.mu2e.cosm0s41b0.Mu2eII.001002_00035055.art" 
"/pnfs/mu2e/tape/phy-sim/dig/mu2e/cosm0s41b0/Mu2eII/art/54/f1/dig.mu2e.cosm0s41b0.Mu2eII.001002_00035555.art" 
"/pnfs/mu2e/tape/phy-sim/dig/mu2e/cosm0s41b0/Mu2eII/art/10/69/dig.mu2e.cosm0s41b0.Mu2eII.001002_00036055.art" 
"/pnfs/mu2e/tape/phy-sim/dig/mu2e/cosm0s41b0/Mu2eII/art/89/0d/dig.mu2e.cosm0s41b0.Mu2eII.001002_00036556.art" 
"/pnfs/mu2e/tape/phy-sim/dig/mu2e/cosm0s41b0/Mu2eII/art/80/0e/dig.mu2e.cosm0s41b0.Mu2eII.001002_00037057.art" 
"/pnfs/mu2e/tape/phy-sim/dig/mu2e/cosm0s41b0/Mu2eII/art/18/12/dig.mu2e.cosm0s41b0.Mu2eII.001002_00037560.art" 
"/pnfs/mu2e/tape/phy-sim/dig/mu2e/cosm0s41b0/Mu2eII/art/c9/22/dig.mu2e.cosm0s41b0.Mu2eII.001002_00038061.art" 
"/pnfs/mu2e/tape/phy-sim/dig/mu2e/cosm0s41b0/Mu2eII/art/7c/ac/dig.mu2e.cosm0s41b0.Mu2eII.001002_00038563.art" 
"/pnfs/mu2e/tape/phy-sim/dig/mu2e/cosm0s41b0/Mu2eII/art/b4/c5/dig.mu2e.cosm0s41b0.Mu2eII.001002_00039065.art" 
"/pnfs/mu2e/tape/phy-sim/dig/mu2e/cosm0s41b0/Mu2eII/art/9a/a0/dig.mu2e.cosm0s41b0.Mu2eII.001002_00039566.art" 
)

odir=/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cosm0s00b0.s4_helix_filter/outstage/35162167

nf=${#files[@]}
# for f in "${files[@]}" ; do
for i in `seq 0 $((nf-1))` ; do
    od=$odir/00/`printf "%05i" $i`

    f=${files[$i]}
    dn=`dirname  $f`
    bn=`basename $f`

    echo $od     $bn

    mu2e -c list.fcl -s $f  2>&1 | grep  "Begin processing the"  >| $od/$bn.event_list
done 
