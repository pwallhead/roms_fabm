GENERAL OVERVIEW
================

The ROMS-FABM coupling has been developed to allow use of the FABM framework for biogeochemical modelling
within the ROMS regional ocean modelling system, see:
https://github.com/fabm-model/fabm/wiki
https://www.myroms.org/

NOTE: It has been developed and tested so far ONLY FOR MPI PARALLEL configurations of ROMS,
      and only for Hedstrom, Rutgers, and COAWST branches of the ROMS code.
NOTE: The Hedstrom ROMS code is no longer maintained, and therefore the FABM coupling code may cease
      to be maintained in the near future. Users are advised to switch to the Rutgers or COAWST code.

To set up ROMS with FABM coupling you need to:

1) Download and install the FABM with ROMS support
2) Download and install the FABM python front end
3) Adapt the FABM model specification file fabm.yaml for the chosen biological model
4) Adapt the python script fabm_roms_...py for the chosen model
5) Generate the ROMS biological model/input files using fabm_roms_...py and fabm.yaml
6) Compile ROMS with FABM coupling using the python-generated biological model files
7) Run ROMS using the python-generated biological input file

Steps 1) and 2) have extensive guidance online, see:
https://github.com/fabm-model/fabm/wiki

Steps 3) - 5) are best done in a dedicated folder fabm_python/
For testing, an example fabm.yaml file for the npzd_Franks model (FABM version) is provided in example_scripts/.
The user should then adapt the example python script fabm_roms_a20_v3_npzd_Franks.py to their specific
ROMS model configuration (the example file is adapted for our ROMS model "A20_v3").
Note that this python script encodes specific choices of boundary condition types
and output types (Hout, Aout, Sout, Dout) for the biological ROMS input file rfabm.in.
These choices can later be modified before submitting the ROMS job.
(See "Using pyfabm to generate ROMS files" below for more details.)
The user can compare a coupled run using the FABM (ROMS-FABM-npzdFranks) with a direct coupling
using the ROMS files for npzd_Franks. We have found very close agreement (discrepancies < 1e-8 mmolN/m3)
and a negligible increase in total wall clock time for the coupled simulation using FABM. 

Step 6) uses the folder Src_<version>_modify/ which contains the ROMS source files 
that have been modified to accommodate the FABM.
These are best incorporated at compile time via the build script (see example_scripts/build_bash_fram_A20.sh).
Otherwise, they could be used to modify files in the appropriate source folders.
If a different/incompatible version of ROMS is to be used:
1) The modifications in Src_<version>_modify files will need to be transferred to the new source files.
For guidance see Src_<version>_modify/notes/roms_<version>_fabm_notes.txt
2) The base files for python-generated files may need to be updated.
3) The python script fabm_roms_....py may need to be updated if line numbers have changed.
NOTE 1: The build script depends on FABM installation folders FABM_INCDIR and FABM_LIBDIR.
        These must be adapted to specify where the user has installed FABM. 
NOTE 2: Files makefile and Linux-<fortran compiler>.mk are platform/compiler dependent
        so they may need to be adapted to the user's specific platform/compiler.
See "Compiling ROMS + FABM-model using Src_<version>_modify" below for more details.

Step 7) will involve a job script. Here I recommend to include commands to copy all input files,
compilation files, module loading scripts, and compiled ROMS binary to the job/scratch directory,
so that the run can be easily and exactly reproduced if necessary. 



Using pyfabm to generate ROMS files
-----------------------------------

0. Install the FABM and FABM python front end
1. Copy your fabm_appname.yaml file (e.g. from ~/git/ersem-edge/testcases) to a python directory e.g. fabm_python/
2. Copy fabm_appname.yaml to runtime directory as fabm.yaml
(I recommend adding commented-out notes at the top of fabm.yaml to document provenance)
3. Generate (rfabm_mod.h, rfabm_var.h, rfabm_inp.h, varinfo_rfabm.dat, rfabm.in, fabm_appname.cdl) using:
python fabm_roms_appname.py fabm_appname.yaml
4. Copy (rfabm_mod.h, rfabm_var.h, rfabm_inp.h) to Src_modify in compile directory (or to ROMS/Nonlinear/Biology in source directory)
5. Copy (varinfo...dat, rfabm.in) to runtime directory

Example:
cp -p rfabm_* ~/Run/Src_modify/.
cp -p rfabm.in ~/Run/.
cp -p varinfo_rfabm.dat ~/Run/.

NOTE: Boundary condition options and source/sink switches in rfabm.in will need to be adapted to the specific ROMS application
      This can be done 1) by hand in the Runtime directory or 2) by adjusting fabm_roms_appname.py
      I use approach 2), and therefore have application-tailored fabm_roms_appname.py files, e.g: fabm_roms_a20_v3_npzd_Franks.py

NOTE: To save time you can do all the cp-ing by sourcing a single shell script, e.g. see: 
      example_scripts/make_fabm_roms_a20_v3_npzd_franks.sh
      (users will need to adapt this script to their folder structures)

TIP: If you get a Segmentation Fault, it could be because you have parameters in the fabm.yaml that
     are not used by the specified modules. A quicker test (than fabm_roms...py) is for example:
       python fabm_describe_model_pwa.py fabm.yaml > tmp.cdl
     If on the other hand you have parameters or couplings that are needed by the specified modules
     and not provided un the fabm.yaml, then you should get a more informative error message. 



Compiling ROMS + FABM-model using Src_<version>_modify
------------------------------------------------------

First, copy the files in Src_<version>_modify/ to a folder Src_modify/ in your compile directory.
Files in Src_modify/ will be used by the build script to substitute for unmodified files in the (pristine) ROMS source directory.
Python-generated source files for ROMS (e.g. rfabm_mod.h) and runtime input files for ROMS
are made using fabm_roms_appname.py (see above).
These are subsequently copied to Src_modify/ (source files) and <Run folder>/
The code is finally compiled by running the build script, e.g.:
./build_bash_fram_A20.sh -j 16

WARNING: Do not attempt to compile multiple versions at the same time!
         This may cause problems due to the exchange of files with the pristine ROMS source directory.

NOTE 1: The ROMS-FABM coupling code may require certain cpp options to be activated in the ROMS header file:
        (see e.g. a20_v3_fabm_npzd_franks.h):
BIOLOGY             /* MUST be defined for use of FABM */
DEBUGFABM           /* use for RFABM-specific debugging output */
DIAGNOSTICS         /* MUST be defined for outputting FABM diagnostics (specified in rfabm.in) */
DIAGNOSTICS_BIO     /* MUST be defined for outputting FABM diagnostics (specified in rfabm.in) */
MASKING             /* MUST be defined to provide input to fabm_set_mask */
ANA_SPFLUX          /* MUST be defined, or surface bgc fluxes provided in forcing files */
ANA_BPFLUX          /* MUST be defined, or bottom bgc fluxes provided in forcing files */
SHORTWAVE           /* MAY be required in order to provide light forcing for FABM model */
FABM_ADYTRACER      /* use to provide a light attenuation tracer via ROMS (input ADY_0, NOT YET TESTED) */
FABM_ASILT          /* use to provide a 3D forcing field for absorption due to silt (input Asilt, NOT YET TESTED) */
FABM_PCO2ATM        /* use to provide atmospheric pCO2 forcing via ROMS (input xCO2atm) */
FABM_N3ATMDEPO      /* use to provide atmospheric deposition flux of oxidized nitrogen via ROMS (input N3atmd) */
FABM_N4ATMDEPO      /* use to provide atmospheric deposition flux of reduced nitrogen via ROMS (input N4atmd) */
FABM_AICE           /* use to provide fractional ice area from ROMS internal ice model (NOT YET TESTED) */
FABM_TSS            /* use to provide Total Suspended Sediments concentration from input file(s) as a forcing (input tss) */
FABM_TSS_ONLINE     /* use to provide Total Suspended Sediments concentration calculated online by ROMS to FABM (NOT YET TESTED) */ 
FABM_NONNEG_S       /* use to cap salinity input to FABM at zero PSU (RECOMMENDED IF SALINITY NOT CAPPED WITHIN FABM MODEL) */
FABM_CHECK_STATE    /* use to cap bgc variable input to FABM (RECOMMENDED) */
FABM_INITIAL        /* use to set all bgc initial conditions to FABM default values (simple constants) */
FABM_INITIAL_SB     /* use to set initial conditions to FABM defaults only for surface/bottom attached variables */

In addition there is also code to allow mass inputs without fluid input, under cpp TS_ISOURCE.
This can be used to simulate inputs from e.g. fish farms, or waste water treatment plants (WWTPs).
TS_ISOURCE is independent of the FABM (thanks to John Wilkin at Rutgers for guidance) and has been successfully TESTED,
but is presently only coded for the Rutgers branch (ARANGO).

There is also some code under a cpp FABM_ISOURCES to provide bgc sources via the FABM, but this is also
NOT YET TESTED and may be deleted in future updates.

FABM_TSS and FABM_TSS_ONLINE options are presently only coded for Rutgers (ARANGO) and COAWST branches.

NOTE 2: If you have a VERY complex FABM model, you may need to increase maximum array size parameters in ROMS.
        We have had to increase maximum size parameters in:
          mod_coupler.F, mod_ncparam.F, lbc.F, tadv.F, base_rfabm_inp.h
        We have also had to increase some dimensions to permit reading of forcings split over MANY files:
          read_phypar.F/Cval (dimension 200->500) (Hedstrom) or value of inp_decode/nCval, 
          def_info/string (dimension 1024->8192 (Hedstrom) or 4096->8192)

NOTE 3: Modified code is provided for ccsm_flux.F to enable the CCSM_ICE_SHORTWAVE parameterization.
        However this code is still under development: it is NOT YET TESTED and NOT READY FOR USE.
        (CCSM_ICE_SHORTWAVE has been tested only with BULK_FLUXES cpp).

NOTE 4: If you are changing between ROMS versions, make sure that all input files are compatible with
        the new version.  For FABM-coupled configurations you should:
        i) Rerun the python script with the version switched 'roms_branch' modified.
           This should ensure that all input files created by the python script are compatible with the new version.
           If you are using an application-adapted base_varinfo.dat file, make sure that the base file used by
           the python script is compatible with your application and the new version.
           For example, cf. base_varinfo_a20_v3_hedstrom.dat and base_varinfo_a20_v3_arango.dat in example_scripts/.
           Note that you may need to adapt the line_split variables to ensure that the FABM entries
           are correctly inserted in the base_varinfo.dat.
        ii) Ensure that the header file is compatible with the new version (some differences in cpps between ROMS branches).
        iii) Modify your build script to use the new source code, and recompile ROMS.
        iv) Check that your ocean.in and stations.in files are compatible with the new version, before running.
            Note that there are major differences in the ocean.in files between Hedstrom and Arango branches;
            if you do not account for these you will likely see a segmentation fault early in the run.
            For example, cf. ocean_a20_v3_fabm_npzd_franks_hedstrom.in and ocean_a20_v3_fabm_npzd_franks_arango.in in example_scripts/.