#!/bin/bash
## This shell script generates and copies across all the FABM files needed to compile and run ROMS-FABM
## for A20.
## Example usage:
## ./make_fabm_roms_a20_v3_npzd_franks.sh fabm_roms_a20_v3_npzd_Franks.py fabm_npzd_Franks.yaml
##
## Phil Wallhead 29/04/2020

#./make_fabm_roms_a20_v3_npzd_Franks.sh $arg1 $arg2

python /cluster/projects/nn9630k/A20/fabm_python/$1 /cluster/projects/nn9630k/A20/fabm_python/$2
#cp -p /cluster/projects/nn9630k/A20/fabm_python/rfabm_* /cluster/projects/nn9630k/A20/Run2/Src_modify/.
cp -p /cluster/projects/nn9630k/A20/fabm_python/rfabm_mod.h /cluster/projects/nn9630k/A20/Run3/Src_modify/.
cp -p /cluster/projects/nn9630k/A20/fabm_python/rfabm_inp.h /cluster/projects/nn9630k/A20/Run3/Src_modify/.
cp -p /cluster/projects/nn9630k/A20/fabm_python/rfabm_var.h /cluster/projects/nn9630k/A20/Run3/Src_modify/.
cp -p /cluster/projects/nn9630k/A20/fabm_python/rfabm_a20_v3_npzd_franks.in /cluster/projects/nn9630k/A20/Run3/.
cp -p /cluster/projects/nn9630k/A20/fabm_python/varinfo_a20_v3_fabm_npzd_franks.dat /cluster/projects/nn9630k/A20/Run3/.
cp -p /cluster/projects/nn9630k/A20/fabm_python/$2 /cluster/projects/nn9630k/A20/Run3/fabm.yaml

