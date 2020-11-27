#!/usr/bin/env python

#This script creates model-dependent .h files for A20 ROMS model with FABM coupling to npzd_Franks biology.
#
#Example usage: python fabm_roms_a20_v3_npzd_Franks.py fabm.yaml
#Output files: rfabm_mod.h, rfabm_inp.h, rfabm_var.h, varinfo_a20_v3.dat, 
#              rfabm_a20_npzd_Franks.in, fabm_npzd_Franks.cdl
#
#Note: If the ROMS version is updated, you may need to update the base files
#      and therefore also the line_split variables in this script (see below)
#
##Phil Wallhead 31/10/2020


import sys, numpy, os

if len(sys.argv)<2:
   print('YAML file name not specified, assuming fabm.yaml')
   yamlfile = 'fabm.yaml'
else:
   yamlfile = sys.argv[1]

try:
   import pyfabm
except ImportError:
   print('Unable to load pyfabm. Please build and install FABM with FABMHOST=python.')
   sys.exit(1)


## Basic simulation options (ROMS branch etc.)
roms_branch = 0        #0 for Hedstrom, 1 for Arango
default_Dout = 0       #1 to use default Dout, 0 for no biological diagnostics
                       #npzd_Franks currently has no diagnostics, must set default_Dout = 0
reduced_DiaBio = 1     #1 to reduce the size of DiaBio2d and DiaBio3d by subsetting the diagnostics at compilation
                       #Otherwise, the allocated array sizes will be set by the total number of possible diagnostics
                       #which may be a burden on memory for complex models with many diagnostics. 


# Create model object from YAML file.
model = pyfabm.Model(yamlfile)

outname = 'rfabm'
outfile1 = outname + '_mod.h'
outfile2 = outname + '_inp.h'
outfile3 = outname + '_var.h'
if roms_branch==0:
   outfile4 = 'varinfo_a20_v3_fabm_npzd_franks.dat'
if roms_branch==1:
   outfile4 = 'varinfo_a20_v3_fabm_npzd_franks.dat'
outfile5 = outname + '_a20_v3_npzd_franks.in'
outfile6 = yamlfile[:len(yamlfile)-5] + '.cdl'

nspace = 6
lmax = 30   #lmax is maximum allowed character length for variable names (<=40)
            #Any names longer than lmax will be replaced by "BS17", "SA34", "BA27", etc.

NBT = len(model.interior_state_variables)
NSAT = len(model.surface_state_variables)
NBAT = len(model.bottom_state_variables)


## APPLICATION SPECIFIC PARAMETERS FOR A20
#Define default boundary conditions for forced and unforced variables
#(application specific, can be modified at run-time)
TNU2str = '120.0d0'
TNU4str = '0.0d0'
AKT_BAKstr = '1.0d-6'
TNUDGstr = '30.0d0'

if roms_branch==1:
    Hadvectionstr = 'HS'
    Vadvectionstr = 'HS'

bc_forced = 'RadNud  Clo     RadNud  RadNud '
bc_unforced = 'Gra     Clo     Gra     Gra    '
isforced_bc = numpy.array([0]*NBT) #Test case, use gradient (unforced) BCs

default_LtracerSrc = 0 #1 to use default LtracerSrc, 0 for no point sources
LtracerSrc = 'T F F F' #Default LtracerSrc

default_LtracerSponge = 0 #1 to use default LtracerSponge, 0 for no sponge layers
LtracerSponge = str(NBT)+'*T' #Default LtracerSponge

#Here we define the "base" varinfo.dat file that is used for physics-only runs
if roms_branch==0:
   brystr = 'ocean_time'
   base_varinfo_name = 'base_varinfo_a20_v3.dat'
if roms_branch==1:
   brystr = 'bry_time'
   base_varinfo_name = 'base_varinfo_a20_v3.dat'

#Below are some example default diagnostics taken from an ERSEM model run.
#For the npzd_Franks model, there are currently no diagnostics, so these are all unused.

#Default oxygen diagnostics (application-dependent)
defDoutO2_2d = ['O2_fair']
defDoutO2_3d = []
#defDoutO2_3d = ['O2_o_sms','O2_o_w','O2_eO2mO2','O2_osat','O2_AOU',
#                'P1_O2o_sms','P2_O2o_sms','P3_O2o_sms','P4_O2o_sms','P5_O2o_sms','P6_O2o_sms',
#                'pel_nit_O2o_sms','Z4_O2o_sms','Z5_O2o_sms','Z6_O2o_sms']

#Default DIC diagnostics (application-dependent)
defDoutDIC_2d = ['O3_fair']
defDoutDIC_3d = []
#defDoutDIC_3d = ['O3_c_sms','O3_c_w','O3_pH','O3_Om_cal','O3_Om_arg',
#                 'P1_netPI','P1_fO3PIc','P1_fPIO3c','P2_netPI','P2_fO3PIc','P2_fPIO3c',
#                 'P3_netPI','P3_fO3PIc','P3_fPIO3c','P4_netPI','P4_fO3PIc','P4_fPIO3c',
#                 'P5_netPI','P5_fO3PIc','P5_fPIO3c','P6_netPI','P6_fO3PIc','P6_fPIO3c',
#                 'Z4_fZIO3c','Z5_fZIO3c','Z6_fZIO3c','B1_fB1O3c']

#Default TA diagnostics (application-dependent)
defDoutTA_2d = []
defDoutTA_3d = []
#defDoutTA_3d = ['O3_TA_sms','O3_TA_w']

#Default DOM diagnostics (application-dependent)
defDoutDOM_2d = []
defDoutDOM_3d = []
#defDoutDOM_3d = ['B1_fB1R1c','B1_fB1R2c','B1_fB1R3c','B1_fB1R1n','B1_fB1R1p',
#                 'B1_fNR1B1c','B1_fNR2B1c','B1_fNR3B1c','B1_fNR1B1n','B1_fNR1B1p',
#                 'B1_fCR1B1c','B1_fCR2B1c','B1_fCR3B1c','B1_fCR1B1n','B1_fCR1B1p',
#                 'P1_fPIR1c','P1_fPIR2c','P2_fPIR1c','P2_fPIR2c','P3_fPIR1c','P3_fPIR2c',
#                 'P4_fPIR1c','P4_fPIR2c','P5_fPIR1c','P5_fPIR2c','P6_fPIR1c','P6_fPIR2c',
#                 'Z4_fZIRDc','Z5_fZIRDc','Z6_fZIRDc']

#Default Aggregrated variable diagnostics (application-dependent)
defDoutAgg_2d = []
defDoutAgg_3d = ['light_Chl','light_POC','B1_DOC']

#Default Light-related diagnostics (application-dependent)
defDoutLight_2d = ['zenithAngle_zenithA','daylength_day_length']
defDoutLight_3d = ['light_PAR0','light_PAR01','light_PAR02','light_PAR03','light_PAR04',
                   'light_PAR05','light_PAR06','light_UV01','light_UV02']
# defDoutLight_3d = ['light_iopABSph1','light_iopABSph2','light_iopABSph3','light_iopABSph4',
#                    'light_iopABSph5','light_iopABSph6','light_iopABSUVph1','light_iopABSUVph2',
#                    'light_iopABS1','light_iopABS2','light_iopABS3','light_iopABS4',
#                    'light_iopABS5','light_iopABS6','light_iopABSUV1','light_iopABSUV2',
#                    'light_iopBBSp1','light_iopBBSp2','light_iopBBSp3','light_iopBBSp4',
#                    'light_iopBBSp5','light_iopBBSp6','light_iopBBSUVp1','light_iopBBSUVp2',
#                    'light_iopBBS1','light_iopBBS2','light_iopBBS3','light_iopBBS4',
#                    'light_iopBBS5','light_iopBBS6','light_iopBBSUV1','light_iopBBSUV2',
#                    'light_K01','light_K02','light_K03','light_K04',
#                    'light_K05','light_K06','light_K0UV1','light_K0UV2',
#                    'light_PAR0','light_PAR01','light_PAR02','light_PAR03','light_PAR04',
#                    'light_PAR05','light_PAR06','light_UV01','light_UV02',
#                    'light_fCR1O3c','light_fCR2O3c','light_fCR3O3c','light_fCR1NInphoto',
#                    'light_fCR1NR1c','light_fCR2NR1c','light_fCR3NR1c',
#                    'P1_Ae','P2_Ae','P3_Ae','P4_Ae','P5_Ae','P6_Ae']

#Default Other diagnostics (application-dependent)
defDoutOther_2d = []
defDoutOther_3d = []
# defDoutOther_3d = ['B1_c_sms','P1_c_sms','P1_c_w','P2_c_sms','P2_c_w','P3_c_sms','P3_c_w',
#                    'P4_c_sms','P4_c_w','P5_c_sms','P5_c_w','P6_c_sms','P6_c_w',
#                    'Z4_c_sms','Z5_c_sms','Z6_c_sms',
#                    'L2_c_sms','L2_c_w','L2_RainR','L2_O3c_sms',
#                    'R4_c_sms','R4_c_w','R6_c_sms','R6_c_w','R8_c_sms','R8_c_w']

if default_Dout==1:
   #Make diagnostic list variables, subsetting if required
   if reduced_DiaBio==1:
      DiaBio_2d_name = (defDoutO2_2d + defDoutDIC_2d + defDoutTA_2d + defDoutDOM_2d + 
                        defDoutAgg_2d + defDoutLight_2d + defDoutOther_2d)
      DiaBio_3d_name = (defDoutO2_3d + defDoutDIC_3d + defDoutTA_3d + defDoutDOM_3d + 
                        defDoutAgg_3d + defDoutLight_3d + defDoutOther_3d)

      NDbio2d = len(DiaBio_2d_name)
      NDbio3d = len(DiaBio_3d_name)
      idfabmD2 = [0] * NDbio2d
      idfabmD3 = [0] * NDbio3d
      DiaBio_2d = [0] * NDbio2d
      DiaBio_3d = [0] * NDbio3d
      for i,v in enumerate(DiaBio_2d_name):
         idfabmD2[i] = next((j for j,x in enumerate(model.horizontal_diagnostic_variables,1) if x.output_name==v), 0)
         if idfabmD2[i]>0:
            DiaBio_2d[i] = model.horizontal_diagnostic_variables[idfabmD2[i]-1]
         else:
            print('Variable name ' + v + ' not found in model.horizontal_diagnostic_variables.output_name')
            sys.exit(1)
      for i,v in enumerate(DiaBio_3d_name):
         idfabmD3[i] = next((j for j,x in enumerate(model.interior_diagnostic_variables,1) if x.output_name==v), 0)
         if idfabmD3[i]>0:
            DiaBio_3d[i] = model.interior_diagnostic_variables[idfabmD3[i]-1]
         else:
            print('Variable name ' + v + ' not found in model.interior_diagnostic_variables.output_name')
            sys.exit(1)
   else:
      NDbio2d = len(model.horizontal_diagnostic_variables)
      NDbio3d = len(model.interior_diagnostic_variables)
      idfabmD2 = range(1,NDbio2d)
      idfabmD3 = range(1,NDbio3d)
      DiaBio_2d = model.horizontal_diagnostic_variables
      DiaBio_3d = model.interior_diagnostic_variables
else:
   NDbio2d = 0
   NDbio3d = 0
   idfabmD2 = []
   idfabmD3 = []
   DiaBio_2d = []
   DiaBio_3d = []









## outname_mod.h
#  This module declares and allocates the biological tracer indices as internal ROMS variables.
with open(outfile1,'w') as f:

   f.write('\n! ROMS input file generated by fabm_roms_a20_v3_npzd_Franks.py\n')
   f.write('!\n! Phil Wallhead, NIVA (2015)\n!\n')
   f.write(' '*nspace + 'USE mod_param\n')
   f.write('!\n')
   f.write(' '*nspace + 'implicit none\n\n\n')


   #Here we declare the biological tracer variable indices
   #These are named according to i<short variable name from .yaml file>
   f.write('!\n!  Set (pelagic) biological tracer identification indices.\n!\n')
   f.write(' '*nspace + 'integer, allocatable :: idbio(:)  ! Biological tracers\n')
   for i,variable in enumerate(model.interior_state_variables,1):
      l1 = len(variable.output_name)
      if l1>lmax:
         vstr = 'BS' + str(i)
         l1 = len(vstr)
      else:
         vstr = variable.output_name

      f.write(' '*nspace + 'integer :: i' + vstr +
              ' '*(lmax-l1) + '! ' + variable.long_path + ' (' + variable.units + ')\n')

   if NSAT>0:
      f.write('!\n!  Set surface-attached state variable identification indices.\n!\n')
      for i,variable in enumerate(model.surface_state_variables,1):
         l1 = len(variable.output_name)
         if l1>lmax:
            vstr = 'SA' + str(i)
            l1 = len(vstr)
         else:
            vstr = variable.output_name

         f.write(' '*nspace + 'integer :: i' + vstr +
                 ' '*(lmax-l1) + '! ' + variable.long_path + ' (' + variable.units + ')\n')

   if NBAT>0:
      f.write('!\n!  Set benthic state variable identification indices.\n!\n')
      for i,variable in enumerate(model.bottom_state_variables,1):
         l1 = len(variable.output_name)
         if l1>lmax:
            vstr = 'BA' + str(i)
            l1 = len(vstr)
         else:
            vstr = variable.output_name

         f.write(' '*nspace + 'integer :: i' + variable.output_name +
                 ' '*(lmax-l1) + '! ' + variable.long_path + ' (' + variable.units + ')\n')

   f.write('\n#if defined DIAGNOSTICS && defined DIAGNOSTICS_BIO\n')
   f.write('!\n!  Biological 2D diagnostic variable IDs.\n!\n')
   f.write(' '*nspace + 'integer, allocatable :: iDbio2(:)\n\n')
   for i,variable in enumerate(DiaBio_2d,1):
      l1 = len(variable.output_name)
      if l1>lmax:
         vstr = 'HD' + str(i)
         l1 = len(vstr)
      else:
         vstr = variable.output_name

      f.write(' '*nspace + 'integer :: i' + vstr + ' = ' + str(i) + '\n')

   f.write('!\n!  Biological 3D diagnostic variable IDs.\n!\n')
   f.write(' '*nspace + 'integer, allocatable :: iDbio3(:)\n\n')
   for i,variable in enumerate(DiaBio_3d,1):
      l1 = len(variable.output_name)
      if l1>lmax:
         vstr = 'BD' + str(i)
         l1 = len(vstr)
      else:
         vstr = variable.output_name

      f.write(' '*nspace + 'integer :: i' + vstr + ' = ' + str(i) + '\n')
   f.write('#endif')

   f.write('!\n!  Biological 2D diagnostic variable indices within FABM.\n!\n')
   f.write(' '*nspace + 'integer, allocatable :: idfabmD2(:)\n\n')

   f.write('!\n!  Biological 3D diagnostic variable indices within FABM.\n!\n')
   f.write(' '*nspace + 'integer, allocatable :: idfabmD3(:)\n\n')

   f.write('\n' + ' '*nspace + 'integer, allocatable :: icheckmax(:)\n')
   f.write(' '*nspace + 'real(r8), allocatable :: dBdt1max(:)\n')
   f.write('\n\n' + ' '*nspace + 'CONTAINS\n')
   f.write(' '*nspace + 'SUBROUTINE initialize_biology\n')
   f.write('!\n!\n!  This routine sets several variables needed by the biology model.\n')
   f.write('!\n!  It allocates and assigns biological tracer indices.\n!\n!\n')
   f.write('!\n!  Local variable declarations\n!\n')
   f.write(' '*nspace + 'integer :: i, ic\n')
   f.write('!\n!  Determine number of biological tracers.\n!\n')
   f.write(' '*nspace + 'NBT=' + str(NBT) +
           '   !Number of pelagic biological tracers\n')
   f.write(' '*nspace + 'NSAT=' + str(NSAT) +
           '   !Number of surface-attached biological state variables (Ngrids vector assignment)\n')
   f.write(' '*nspace + 'NBAT=' + str(NBAT) +
           '   !Number of bottom-attached biological state variables (Ngrids vector assignment)\n')
   f.write('!Note: NBT, NSAT, and NBAT are all declared in mod_param\n')
   f.write('!      This is so that other subroutines to access these parameters via mod_param\n')
   f.write('!      However, we need to use (NSAT, NBAT) already below, hence these are initialized here\n\n')

   f.write('\n#if defined DIAGNOSTICS && defined DIAGNOSTICS_BIO\n')
   f.write('!\n!-----------------------------------------------------------------------\n')
   f.write('! Set sources and sinks biology diagnostic parameters.\n')
   f.write('!-----------------------------------------------------------------------\n')
   f.write('!\n! Set number of diagnostic terms.\n!\n')
   f.write(' '*nspace + 'NDbio2d=' + str(NDbio2d) + '\n')
   f.write(' '*nspace + 'NDbio3d=' + str(NDbio3d) + '\n')
   f.write('#endif\n\n')

   f.write('!\n!  Allocate various module variables.\n!\n')
   f.write(' '*nspace + 'IF (.not.allocated(icheckmax)) THEN\n')
   f.write(' '*nspace + '  allocate ( icheckmax(Ngrids) )\n')
   f.write(' '*nspace + 'END IF\n')
   f.write(' '*nspace + 'IF (.not.allocated(dBdt1max)) THEN\n')
   f.write(' '*nspace + '  allocate ( dBdt1max(Ngrids) )\n')
   f.write(' '*nspace + 'END IF\n')

   #Allocate pelagic biological tracer index vector
   f.write('!\n!  Allocate pelagic biological tracer index vector.\n!\n')
   f.write('!  Note: these are indices within the NetCDF index vector idTvar\n!\n')
   f.write(' '*nspace + 'IF (.not.allocated(idbio)) THEN\n' + ' '*nspace +
           '  allocate ( idbio(NBT) )\n' + ' '*nspace + 'END IF\n')

   #Allocate biological diagnostics NetCDF index vectors
   f.write('\n#if defined DIAGNOSTICS && defined DIAGNOSTICS_BIO\n')
   f.write('!\n!  Allocate biological diagnostics NetCDF index vectors iDbio2 and iDbio3.\n')
   f.write('!  Note: these are analagous to idTvar, not subindex vectors like idbio above;\n')
   f.write('!        they are assigned directly in rfabm_var.h, using the individual subindices set above.\n!\n')
   f.write(' '*nspace + 'IF (.not.allocated(iDbio2)) THEN\n' + ' '*nspace +
           '  allocate ( iDbio2(NDbio2d) )\n' + ' '*nspace + 'END IF\n\n')
   f.write(' '*nspace + 'IF (.not.allocated(iDbio3)) THEN\n' + ' '*nspace +
           '  allocate ( iDbio3(NDbio3d) )\n' + ' '*nspace + 'END IF\n')
   f.write('!\n!  Allocate FABM index vectors idfabmD2 and idfabmD3.\n')
   f.write('!  These are the indices of the compiled diagnostic subsets in the complete FABM lists\n!\n')
   f.write(' '*nspace + 'IF (.not.allocated(idfabmD2)) THEN\n' + ' '*nspace +
           '  allocate ( idfabmD2(NDbio2d) )\n' + ' '*nspace + 'END IF\n\n')
   f.write(' '*nspace + 'IF (.not.allocated(idfabmD3)) THEN\n' + ' '*nspace +
           '  allocate ( idfabmD3(NDbio3d) )\n' + ' '*nspace + 'END IF\n')

   f.write('#endif\n\n')
   f.write('#ifdef DEBUG\n' + ' '*nspace + 'write(*,*) "Done all allocation in ' + outfile1 + '"\n#endif\n')

   #Initialize pelagic biological tracer index vector
   f.write('!\n!  Initialize pelagic tracer identification indices.\n!\n')
   f.write(' '*nspace + 'ic=NAT+NPT+NCS+NNS\n')
   f.write(' '*nspace + 'DO i=1,NBT\n')
   f.write(' '*nspace + '  idbio(i)=ic+i\n')
   f.write(' '*nspace + 'END DO\n\n')

   #Initialize individual pelagic tracer indices
   f.write('!  Here we set the indices in the tracer array "t"\n')
   for i,variable in enumerate(model.interior_state_variables,1):
      l1 = len(variable.output_name)
      if l1>lmax:
         vstr = 'BS' + str(i)
         l1 = len(vstr)
      else:
         vstr = variable.output_name

      f.write(' '*nspace + 'i' + vstr + '=ic+' + str(i) + '\n')

   if NSAT>0:
      #Initialize individual surface-attached state variable indices
      f.write('\n!  Here we set the indices in the tracer array "state_sf"\n')
      for i,variable in enumerate(model.surface_state_variables,1):
         l1 = len(variable.output_name)
         if l1>lmax:
            vstr = 'SA' + str(i)
            l1 = len(vstr)
         else:
            vstr = variable.output_name

         f.write(' '*nspace + 'i' + vstr + '=' + str(i) + '\n')

   if NBAT>0:
      #Initialize individual benthic state variable indices
      f.write('\n!  Here we set the indices in the tracer array "state_bt"\n')
      for i,variable in enumerate(model.bottom_state_variables,1):
         l1 = len(variable.output_name)
         if l1>lmax:
            vstr = 'BA' + str(i)
            l1 = len(vstr)
         else:
            vstr = variable.output_name

         f.write(' '*nspace + 'i' + vstr + '=' + str(i) + '\n')

   f.write('\n#if defined DIAGNOSTICS && defined DIAGNOSTICS_BIO\n')
   f.write('!\n! Set 2D diagnostic indices within the complete FABM list.\n!\n')
   for i,ind in enumerate(idfabmD2,1):
      f.write(' '*nspace + 'idfabmD2(' + str(i) + ')=' + str(ind) + '\n')
   f.write('!\n! Set 3D diagnostic indices within the complete FABM list.\n!\n')
   for i,ind in enumerate(idfabmD3,1):
      f.write(' '*nspace + 'idfabmD3(' + str(i) + ')=' + str(ind) + '\n')
   f.write('#endif\n\n')

   f.write('\n' + ' '*nspace + 'RETURN\n')
   f.write(' '*nspace + 'END SUBROUTINE initialize_biology')




## outname_def.h
#  This code defines internal biological parameter names in output NetCDF files.
#  With FABM coupling this parameter set is EMPTY, so this file is FIXED.




## outname_inp.h
#  This code assigns values to the ROMS internal biological parameters set in rfabm.in.
#  NOTE: FABM model parameter values are set in the .yaml file
#        rfabm.in is only used to set transport parameters and output / diagnostic options
if roms_branch==0:
   base_inp = 'base_v37_KATE/base_rfabm_inp.h'
   line_split = 471
if roms_branch==1:
   base_inp = 'base_Arango/base_rfabm_inp.h'
   line_split = 552 #500 with 26/07/2019 code
with open(base_inp,'r') as f:
   data = f.readlines()
#Here we split the base file into two files at a suitable point
   data1 = data[0:line_split]
   data2 = data[line_split+1:]

nspace = 12
with open(outfile2,'w') as f:
   f.writelines(data1)

   f.write('\n\n!!! PWA: HERE IS WHERE THE FABM-DEPENDENT INFO IS INSERTED !!!!!!\n')

   f.write('#ifdef DIAGNOSTICS_BIO\n')
   f.write('!\n!  FABM horizontal diagnostic variables\n!\n\n')
   for i,variable in enumerate(DiaBio_2d,1):
      l1 = len(variable.output_name)
      if l1>lmax:
         vstr = 'HD' + str(i)
         l1 = len(vstr)
      else:
         vstr = variable.output_name

      f.write(' '*nspace + 'CASE(' + ' '*55 + '&\n')
      f.write('     & \'Dout(i' + vstr + ')\')\n')
      f.write(' '*nspace + '  IF (iDbio2(' + ' '*47 + '&\n')
      f.write('     & i' + vstr + ').eq.0) THEN\n')
      f.write(' '*nspace + '    IF (Master) WRITE (out,40)' + ' '*30 + '&\n')
      f.write('     & \'iDbio2(i' + vstr + ')\'\n')
      f.write(' '*nspace + '    exit_flag=5\n')
      f.write(' '*nspace + '    RETURN\n')
      f.write(' '*nspace + '  END IF\n')
      f.write(' '*nspace + '  Npts=load_l(Nval, Cval, Ngrids, Lbio)\n')
      f.write(' '*nspace + '  i=iDbio2(' + ' '*49 + '&\n')
      f.write('     & i' + vstr + ')\n')
      f.write(' '*nspace + '  DO ng=1,Ngrids\n')
      f.write(' '*nspace + '    Dout(i,ng)=Lbio(ng)\n')
      f.write(' '*nspace + '  END DO\n')

   f.write('\n!\n!  FABM interior diagnostic variables\n!\n\n')
   for i,variable in enumerate(DiaBio_3d,1):
      l1 = len(variable.output_name)
      if l1>lmax:
         vstr = 'BD' + str(i)
         l1 = len(vstr)
      else:
         vstr = variable.output_name

      f.write(' '*nspace + 'CASE(' + ' '*55 + '&\n')
      f.write('     & \'Dout(i' + vstr + ')\')\n')
      f.write(' '*nspace + '  IF (iDbio3(' + ' '*47 + '&\n')
      f.write('     & i' + vstr + ').eq.0) THEN\n')
      f.write(' '*nspace + '    IF (Master) WRITE (out,40)' + ' '*30 + '&\n')
      f.write('     & \'iDbio3(i' + vstr + ')\'\n')
      f.write(' '*nspace + '    exit_flag=5\n')
      f.write(' '*nspace + '    RETURN\n')
      f.write(' '*nspace + '  END IF\n')
      f.write(' '*nspace + '  Npts=load_l(Nval, Cval, Ngrids, Lbio)\n')
      f.write(' '*nspace + '  i=iDbio3(' + ' '*49 + '&\n')
      f.write('     & i' + vstr + ')\n')
      f.write(' '*nspace + '  DO ng=1,Ngrids\n')
      f.write(' '*nspace + '    Dout(i,ng)=Lbio(ng)\n')
      f.write(' '*nspace + '  END DO\n')

   f.write('# ifdef DEBUG\n' + ' '*nspace + 'write(*,*) "Done input of diagnostic output options Dout"\n# endif\n')
   f.write('#endif\n')
   f.write('\n!!!PWA FABM-DEPENDENT INFO INSERTED ABOVE !!!!!!\n\n\n')
   f.writelines(data2)



## outname_var.h
#  This code defines the metadata indices e.g. idTvar(iOXY_) that are used for NetCDF IO.
#  Metadata info is read from varinfo.dat
#  For example, in mod_ncparam.F/initialize_ncparam, varinfo.dat is read and
#     internal ID variables corresponding to the strings in Vinfo(6) are assigned an integer ID (varid).
#     The biological indices are grouped into a vector idTvar(). To assign named indices
#     within this vector e.g. idTvar(iDOM_), we must have model-dependent code (include file)
#     because iDOM_ is model-specific. Hence we have include files e.g. oxydep_var.h.
#  With FABM coupling this file must therefore be adapted to the number of model variables and their index names.
nspace2 = 14
with open(outfile3,'w') as f:

   f.write('/*\n**  ROMS input file generated by fabm_roms_a20_v3_npzd_Franks.py                        **\n')
   f.write('**  Assigns metadata indices for the FABM model variables that are used in NetCDF IO  **\n')
   f.write('**  The metadata information is read from "varinfo.dat"                               **\n')
   f.write('**  This file is included in file "mod_ncparam.F", routine "initialize_ncparm"        **\n')
   f.write('**  Phil Wallhead, NIVA (2015)  **\n*/\n\n')

   f.write('/*\n**  Model state biological tracers.\n*/\n')
   f.write('#ifdef DEBUG\n' + ' '*nspace + 'write(*,*) "Assigning idTvar in rfabm_var.h"\n#endif\n')
   for i,variable in enumerate(model.interior_state_variables,1):
      l1 = len(variable.output_name)
      if l1>lmax:
         vstr = 'BS' + str(i)
         l1 = len(vstr)
      else:
         vstr = variable.output_name

      f.write(' '*nspace2 + 'CASE(\'idTvar(i' + vstr + ')\')\n')
      f.write(' '*nspace2 + '  idTvar(i' + vstr + ') = varid\n')

   if NSAT>0:
      f.write('\n/*\n**  Surface-attached variables.\n*/\n')
      f.write('#ifdef DEBUG\n' + ' '*nspace + 'write(*,*) "Assigning idSAvar in rfabm_var.h"\n#endif\n')
      for i,variable in enumerate(model.surface_state_variables,1):
         l1 = len(variable.output_name)
         if l1>lmax:
            vstr = 'SA' + str(i)
            l1 = len(vstr)
         else:
            vstr = variable.output_name

         f.write(' '*nspace2 + 'CASE(\'idSAvar(i' + vstr + ')\')\n')
         f.write(' '*nspace2 + '  idSAvar(i' + vstr + ') = varid\n')

   if NBAT>0:
      f.write('\n/*\n**  Bottom-attached variables.\n*/\n')
      f.write('#ifdef DEBUG\n' + ' '*nspace + 'write(*,*) "Assigning idBAvar in rfabm_var.h"\n#endif\n')
      for i,variable in enumerate(model.bottom_state_variables,1):
         l1 = len(variable.output_name)
         if l1>lmax:
            vstr = 'BA' + str(i)
            l1 = len(vstr)
         else:
            vstr = variable.output_name

         f.write(' '*nspace2 + 'CASE(\'idBAvar(i' + vstr + ')\')\n')
         f.write(' '*nspace2 + '  idBAvar(i' + vstr + ') = varid\n')

   f.write('\n/*\n**  Adjoint sensitivity state biological tracers.\n*/\n\n')
   f.write('#if defined AD_SENSITIVITY   || defined IS4DVAR_SENSITIVITY || \\\n')
   f.write('    defined OPT_OBSERVATIONS || defined SENSITIVITY_4DVAR   || \\\n')
   f.write('    defined SO_SEMI\n')
   f.write('\n/*\n**  Adjoint sensitivity.\n*/\n')
   f.write('#ifdef DEBUG\n' + ' '*nspace + 'write(*,*) "Assigning idTads in rfabm_var.h"\n#endif\n')
   for i,variable in enumerate(model.interior_state_variables,1):
      l1 = len(variable.output_name)
      if l1>lmax:
         vstr = 'BS' + str(i)
         l1 = len(vstr)
      else:
         vstr = variable.output_name

      f.write(' '*nspace2 + 'CASE(\'idTads(i' + vstr + ')\')\n')
      f.write(' '*nspace2 + '  idTads(i' + vstr + ') = varid\n')
   f.write('#endif')

   f.write('\n/*\n**  Biological tracers open boundary conditions.\n*/\n')
   f.write('#ifdef DEBUG\n' + ' '*nspace + 'write(*,*) "Assigning idTbry in rfabm_var.h"\n#endif\n')
   for i,variable in enumerate(model.interior_state_variables,1):
      l1 = len(variable.output_name)
      if l1>lmax:
         vstr = 'BS' + str(i)
         l1 = len(vstr)
      else:
         vstr = variable.output_name

      f.write(' '*nspace2 + 'CASE(\'idTbry(iwest,i' + vstr + ')\')\n')
      f.write(' '*nspace2 + '  idTbry(iwest,i' + vstr + ') = varid\n')
      f.write(' '*nspace2 + 'CASE(\'idTbry(ieast,i' + vstr + ')\')\n')
      f.write(' '*nspace2 + '  idTbry(ieast,i' + vstr + ') = varid\n')
      f.write(' '*nspace2 + 'CASE(\'idTbry(isouth,i' + vstr + ')\')\n')
      f.write(' '*nspace2 + '  idTbry(isouth,i' + vstr + ') = varid\n')
      f.write(' '*nspace2 + 'CASE(\'idTbry(inorth,i' + vstr + ')\')\n')
      f.write(' '*nspace2 + '  idTbry(inorth,i' + vstr + ') = varid\n\n')

   f.write('\n/*\n**  Biological tracers point Source/Sinks (river runoff)..\n*/\n')
   f.write('#ifdef DEBUG\n' + ' '*nspace + 'write(*,*) "Assigning idRtrc in rfabm_var.h"\n#endif\n')
   for i,variable in enumerate(model.interior_state_variables,1):
      l1 = len(variable.output_name)
      if l1>lmax:
         vstr = 'BS' + str(i)
         l1 = len(vstr)
      else:
         vstr = variable.output_name

      f.write(' '*nspace2 + 'CASE(\'idRtrc(i' + vstr + ')\')\n')
      f.write(' '*nspace2 + '  idRtrc(i' + vstr + ') = varid\n')

   f.write('\n#ifdef DIAGNOSTICS_BIO\n')
   f.write('!\n!  FABM horizontal diagnostic variables\n!\n\n')
   f.write('#ifdef DEBUG\n' + ' '*nspace + 'write(*,*) "Assigning iDbio2 in rfabm_var.h"\n#endif\n')
   for i,variable in enumerate(DiaBio_2d,1):
      l1 = len(variable.output_name)
      if l1>lmax:
         vstr = 'HD' + str(i)
         l1 = len(vstr)
      else:
         vstr = variable.output_name

      f.write(' '*nspace + 'CASE(' + ' '*55 + '&\n')
      f.write('     & \'iDbio2(i' + vstr + ')\')\n')
      f.write(' '*nspace + '  iDbio2(' + ' '*51 + '&\n')
      f.write('     & i' + vstr + ')=varid\n')

   f.write('\n!\n!  FABM interior diagnostic variables\n!\n\n')
   f.write('#ifdef DEBUG\n' + ' '*nspace + 'write(*,*) "Assigning iDbio3 in rfabm_var.h"\n#endif\n')
   for i,variable in enumerate(DiaBio_3d,1):
      l1 = len(variable.output_name)
      if l1>lmax:
         vstr = 'BD' + str(i)
         l1 = len(vstr)
      else:
         vstr = variable.output_name

      f.write(' '*nspace + 'CASE(' + ' '*55 + '&\n')
      f.write('     & \'iDbio3(i' + vstr + ')\')\n')
      f.write(' '*nspace + '  iDbio3(' + ' '*51 + '&\n')
      f.write('     & i' + vstr + ')=varid\n')

   f.write('#endif\n')
   f.write('#ifdef DEBUG\n' + ' '*nspace + 'write(*,*) "Done assigning ids in rfabm_var.h"\n#endif\n')



## varinfo.dat
if roms_branch==0:
   base_varinfo = 'base_v37_KATE/' + base_varinfo_name
   line_split = 3871 #NOTE: This must be adapted to the user's varinfo.dat file
if roms_branch==1:
   base_varinfo = 'base_Arango/' + base_varinfo_name
   line_split = 2862 #NOTE: This must be adapted to the user's varinfo.dat file
with open(base_varinfo,'r') as f: 
   data = f.readlines()

#Here we split the base file into two files at a suitable point
#that will avoid SIGSEGV errors generated by tracer2 if AVERAGES is defined
#One such point is after "wet_dry_masking" and before Fennel model variables
   data1 = data[0:line_split]
   data2 = data[line_split+1:]

nspace = 45
with open(outfile4,'w') as f:
   f.writelines(data1)

   f.write('\n\n!!! PWA: HERE IS WHERE THE FABM-DEPENDENT INFO IS INSERTED !!!!!!\n')
   f.write('\n!\n!  FABM model boundary variables\n!\n\n')
   for i,variable in enumerate(model.interior_state_variables,1):
      l1 = len(variable.output_name)
      if l1>lmax:
         vstr = 'BS' + str(i)
         l1 = len(vstr)
      else:
         vstr = variable.output_name

      f.write('\'' + vstr + '_east\'' + ' '*nspace +'! Input\n')
      f.write('  \'' + variable.long_path + ' eastern boundary condition\'\n')
      f.write('  \'' + variable.units + '\'\n')
      f.write('  \'' + vstr + '_east, scalar, series\'\n')
      f.write('  \''+brystr+'\'\n')
      f.write('  \'idTbry(ieast,i' + vstr + ')\'\n')
      f.write('  \'nulvar\'\n')
      f.write('  1.0d0\n\n')

      f.write('\'' + vstr + '_west\'' + ' '*nspace +'! Input\n')
      f.write('  \'' + variable.long_path + ' western boundary condition\'\n')
      f.write('  \'' + variable.units + '\'\n')
      f.write('  \'' + vstr + '_west, scalar, series\'\n')
      f.write('  \''+brystr+'\'\n')
      f.write('  \'idTbry(iwest,i' + vstr + ')\'\n')
      f.write('  \'nulvar\'\n')
      f.write('  1.0d0\n\n')

      f.write('\'' + vstr + '_south\'' + ' '*nspace +'! Input\n')
      f.write('  \'' + variable.long_path + ' southern boundary condition\'\n')
      f.write('  \'' + variable.units + '\'\n')
      f.write('  \'' + vstr + '_south, scalar, series\'\n')
      f.write('  \''+brystr+'\'\n')
      f.write('  \'idTbry(isouth,i' + vstr + ')\'\n')
      f.write('  \'nulvar\'\n')
      f.write('  1.0d0\n\n')

      f.write('\'' + vstr + '_north\'' + ' '*nspace +'! Input\n')
      f.write('  \'' + variable.long_path + ' northern boundary condition\'\n')
      f.write('  \'' + variable.units + '\'\n')
      f.write('  \'' + vstr + '_north, scalar, series\'\n')
      f.write('  \''+brystr+'\'\n')
      f.write('  \'idTbry(inorth,i' + vstr + ')\'\n')
      f.write('  \'nulvar\'\n')
      f.write('  1.0d0\n\n')

   f.write('\n!\n!  FABM interior model variables\n!\n')
   f.write('!!!NOTE PWA: These MUST appear BEFORE tracer2 if AVERAGES is defined\n\n')
   for i,variable in enumerate(model.interior_state_variables,1):
      l1 = len(variable.output_name)
      if l1>lmax:
         vstr = 'BS' + str(i)
         l1 = len(vstr)
      else:
         vstr = variable.output_name

      f.write('\'' + vstr + '\'' + ' '*nspace + '! Input/Output\n')
      f.write('  \'' + variable.long_path + '\'\n')
      f.write('  \'' + variable.units + '\'\n')
      f.write('  \'' + vstr + ', scalar, series\'\n')
      f.write('  \'ocean_time\'\n')
      f.write('  \'idTvar(i' + vstr + ')\'\n')
      f.write('  \'r3dvar\'\n')
      f.write('  1.0d0\n\n')

   if NSAT>0:
      f.write('\n!\n!  FABM surface attached variables\n!\n')
      for i,variable in enumerate(model.surface_state_variables,1):
         l1 = len(variable.output_name)
         if l1>lmax:
            vstr = 'SA' + str(i)
            l1 = len(vstr)
         else:
            vstr = variable.output_name

         f.write('\'' + vstr + '\'' + ' '*nspace + '! Input/Output\n')
         f.write('  \'' + variable.long_path + '\'\n')
         f.write('  \'' + variable.units + '\'\n')
         f.write('  \'' + vstr + ', scalar, series\'\n')
         f.write('  \'ocean_time\'\n')
         f.write('  \'idSAvar(i' + vstr + ')\'\n')
         f.write('  \'r2dvar\'\n')
         f.write('  1.0d0\n\n')

   if NBAT>0:
      f.write('\n!\n!  FABM bottom attached variables\n!\n')
      for i,variable in enumerate(model.bottom_state_variables,1):
         l1 = len(variable.output_name)
         if l1>lmax:
            vstr = 'BA' + str(i)
            l1 = len(vstr)
         else:
            vstr = variable.output_name

         f.write('\'' + vstr + '\'' + ' '*nspace + '! Input/Output\n')
         f.write('  \'' + variable.long_path + '\'\n')
         f.write('  \'' + variable.units + '\'\n')
         f.write('  \'' + vstr + ', scalar, series\'\n')
         f.write('  \'ocean_time\'\n')
         f.write('  \'idBAvar(i' + vstr + ')\'\n')
         f.write('  \'r2dvar\'\n')
         f.write('  1.0d0\n\n')

   f.write('\n!\n!  FABM model point sources/sinks (river runoff)\n!\n\n')
   for i,variable in enumerate(model.interior_state_variables,1):
      l1 = len(variable.output_name)
      if l1>lmax:
         vstr = 'BS' + str(i)
         l1 = len(vstr)
      else:
         vstr = variable.output_name

      f.write('\'river_' + vstr + '\'' + ' '*nspace + '! Input\n')
      f.write('  \'river runoff ' + variable.long_path + '\'\n')
      f.write('  \'' + variable.units + '\'\n')
      f.write('  \'' + vstr + ', scalar, series\'\n')
      f.write('  \'river_time\'\n')
      f.write('  \'idRtrc(i' + vstr + ')\'\n')
      f.write('  \'nulvar\'\n')
      f.write('  1.0d0\n\n')

   f.write('\n!\n!  FABM model forcing variables\n!\n\n')

   f.write('\'' + 'swradWm2' + '\'' + ' '*nspace + '! Output\n')
   f.write('  \'' + 'solar shortwave radiation flux (supplied to FABM)' + '\'\n')
   f.write('  \'' + 'watt meter-2' + '\'\n')
   f.write('  \'' + 'shortwave radiation' + ', scalar, series\'\n')
   f.write('  \'ocean_time\'\n')
   f.write('  \'idSradWm2\'\n')
   f.write('  \'r2dvar\'\n')
   f.write('  1.0d0\n\n')

   f.write('\'' + 'Asilt' + '\'' + ' '*nspace + '! Input/Output\n')
   f.write('  \'' + 'absorption due to silt' + '\'\n')
   f.write('  \'' + 'm-1' + '\'\n')
   f.write('  \'' + 'Asilt_time' + ', scalar, series\'\n')
   f.write('  \'Asilt_time\'\n')
   f.write('  \'idAsilt\'\n')
   f.write('  \'r3dvar\'\n')
   f.write('  1.0d0\n\n')

   f.write('\'' + 'xCO2atm' + '\'' + ' '*nspace + '! Input/Output\n')
   f.write('  \'' + 'dry air mole fraction of CO2 in air' + '\'\n')
   f.write('  \'' + 'ppm' + '\'\n')
   f.write('  \'' + 'xCO2atm' + ', scalar, series\'\n')
   f.write('  \'xCO2atm_time\'\n')
   f.write('  \'idxCO2atm\'\n')
   f.write('  \'r2dvar\'\n')
   f.write('  1.0d0\n\n')

   f.write('\'' + 'pCO2atm' + '\'' + ' '*nspace + '! Input/Output\n')
   f.write('  \'' + 'partial pressure of CO2 in air' + '\'\n')
   f.write('  \'' + 'uatm' + '\'\n')
   f.write('  \'' + 'pCO2atm' + ', scalar, series\'\n')
   f.write('  \'xCO2atm_time\'\n')
   f.write('  \'idpCO2atm\'\n')
   f.write('  \'r2dvar\'\n')
   f.write('  1.0d0\n\n')

   f.write('\'' + 'ADY_0' + '\'' + ' '*nspace + '! Input/Output\n')
   f.write('  \'' + 'absorption due to dissolved+detrital matter' + '\'\n')
   f.write('  \'' + 'm-1' + '\'\n')
   f.write('  \'' + 'ADY_0' + ', scalar, series\'\n')
   f.write('  \'ADY_0_time\'\n')
   f.write('  \'idADY_0\'\n')
   f.write('  \'r2dvar\'\n')
   f.write('  1.0d0\n\n')

   f.write('\'' + 'N3atmd' + '\'' + ' '*nspace + '! Input/Output\n')
   f.write('  \'' + 'atmospheric deposition of oxidized nitrogen' + '\'\n')
   f.write('  \'' + 'mmolN/m2/s' + '\'\n')
   f.write('  \'' + 'N3atmd' + ', scalar, series\'\n')
   f.write('  \'atmd_time\'\n')
   f.write('  \'idN3atmd\'\n')
   f.write('  \'r2dvar\'\n')
   f.write('  1.0d0\n\n')

   f.write('\'' + 'N4atmd' + '\'' + ' '*nspace + '! Input/Output\n')
   f.write('  \'' + 'atmospheric deposition of reduced nitrogen' + '\'\n')
   f.write('  \'' + 'mmolN/m2/s' + '\'\n')
   f.write('  \'' + 'N4atmd' + ', scalar, series\'\n')
   f.write('  \'atmd_time\'\n')
   f.write('  \'idN4atmd\'\n')
   f.write('  \'r2dvar\'\n')
   f.write('  1.0d0\n\n')

   #Interior sources (e.g. fish farms)
   for i,variable in enumerate(model.interior_state_variables,1):
      l1 = len(variable.output_name)
      if l1>lmax:
         vstr = 'BS' + str(i)
         l1 = len(vstr)
      else:
         vstr = variable.output_name
      f.write('\'isource_' + vstr + '\'' + ' '*nspace + '! Input\n')
      f.write('  \'interior source ' + variable.long_path + '\'\n')
      f.write('  \'' + variable.units + '/s\'\n')
      f.write('  \'' + vstr + ', scalar, series\'\n')
      f.write('  \'isource_time\'\n')
      f.write('  \'idIsources(i' + vstr + ')\'\n')
      f.write('  \'nulvar\'\n')
      f.write('  1.0d0\n\n')

   f.write('\n!\n!  FABM horizontal diagnostic variables\n!\n\n')
   for i,variable in enumerate(DiaBio_2d,1):
      l1 = len(variable.output_name)
      if l1>lmax:
         vstr = 'HD' + str(i)
         l1 = len(vstr)
      else:
         vstr = variable.output_name

      f.write('\'' + vstr + '\'' + ' '*nspace + '! Input/Output\n')
      f.write('  \'' + variable.long_path + '\'\n')
      f.write('  \'' + variable.units + '\'\n')
      f.write('  \'' + vstr + ', scalar, series\'\n')
      f.write('  \'ocean_time\'\n')
      f.write('  \'iDbio2(i' + vstr + ')\'\n')
      f.write('  \'r2dvar\'\n')
      f.write('  1.0d0\n\n')

   f.write('\n!\n!  FABM interior diagnostic variables\n!\n\n')
   for i,variable in enumerate(DiaBio_3d,1):
      l1 = len(variable.output_name)
      if l1>lmax:
         vstr = 'BD' + str(i)
         l1 = len(vstr)
      else:
         vstr = variable.output_name

      f.write('\'' + vstr + '\'' + ' '*nspace + '! Input/Output\n')
      f.write('  \'' + variable.long_path + '\'\n')
      f.write('  \'' + variable.units + '\'\n')
      f.write('  \'' + vstr + ', scalar, series\'\n')
      f.write('  \'ocean_time\'\n')
      f.write('  \'iDbio3(i' + vstr + ')\'\n')
      f.write('  \'r3dvar\'\n')
      f.write('  1.0d0\n\n')

   f.write('!!!PWA FABM-DEPENDENT INFO INSERTED ABOVE !!!!!!\n\n\n')
   f.writelines(data2)




## outname_wrt.h
#  This code writes the internal biological parameter values to NetCDF output files.
#  With FABM coupling this parameter set is EMPTY, so this file is FIXED.




## rfabm.in
#The code below writes a "default" rfabm.in file -- this can be modified at run time
if roms_branch==0:
   base_in = 'base_v37_KATE/base_rfabm.in'
   line_split1 = 39
   line_split2 = 74
if roms_branch==1:
   base_in = 'base_Arango/base_rfabm.in'
   line_split1 = 40
   line_split2 = 61
   line_split3 = 96

with open(base_in,'r') as f:    #Adjust if varinfo base file is updated
   data = f.readlines()
   data1 = data[0:line_split1] # File introductory comments
   
   if roms_branch==0:
      data2 = data[line_split1+1:line_split2] # Text explaining LBC notation
      data3 = data[line_split2+1:] # Glossary
   if roms_branch==1:
      data2 = data[line_split1+1:line_split2] # Text explaining advection scheme switches
      data3 = data[line_split2+1:line_split3] # Text explaining LBC notation
      data4 = data[line_split3+1:] # Glossary

nspace = 20

with open(outfile5,'w') as f:
   f.writelines(data1)

   f.write('\n! Switch to control the computation of biology within nested and/or multiple\n')
   f.write('! connected grids.\n\n')
   f.write('    Lbiology == T\n\n')
   
   f.write('! Number of initial iterations for which FABM-supplied SMS is checked\n\n')
   f.write('   icheckmax == 0\n\n') # Modify to activate initial FABM output checks

   f.write('! Maximum acceptable absolute value of the FABM-supplied SMS, if checked\n\n')
   f.write('    dBdt1max == 1.0d3\n\n') # Modify for initial FABM output checks

   f.write('! Harmonic/biharmonic horizontal diffusion of biological tracer for\n')
   f.write('! nonlinear model and adjoint-based algorithms: [1:NBT,Ngrids].\n\n')
   f.write('        TNU2 == ' + str(NBT) + '*' + TNU2str + '                        ! m2/s\n')
   f.write('        TNU4 == ' + str(NBT) + '*' + TNU4str + '                        ! m4/s\n\n')
   f.write('     ad_TNU2 == ' + str(NBT) + '*' + TNU2str + '                        ! m2/s\n')
   f.write('     ad_TNU4 == ' + str(NBT) + '*' + TNU4str + '                        ! m4/s\n\n')

   f.write('! Logical switches (TRUE/FALSE) to increase/decrease horizontal diffusivity\n')
   f.write('! in specific areas of the application domain (like sponge areas) for the\n')
   f.write('! desired grid: [1:NBT,Ngrids] ![Ngrids]\n\n')

   if default_LtracerSponge==1:
      f.write('  LtracerSponge == ' + LtracerSponge + '\n\n')
   else:
      f.write('  LtracerSponge == ' + str(NBT) + '*F\n\n')

   f.write('! Vertical mixing coefficients for biological tracers for nonlinear\n')
   f.write('! model and basic state scale factor in adjoint-based algorithms:\n')
   f.write('! [1:NBT,Ngrids].\n\n')
   f.write('     AKT_BAK == ' + str(NBT) + '*' + AKT_BAKstr + '                       ! m2/s\n')
   f.write('  ad_AKT_fac == ' + str(NBT) + '*1.0d0                        ! nondimensional\n\n')

   f.write('! Nudging/relaxation time scales, inverse scales will be computed\n')
   f.write('! internally: [1:NBT,Ngrids].\n\n')
   f.write('       TNUDG == ' + str(NBT) + '*' + TNUDGstr + '                       ! days\n\n')

   if roms_branch==0:
      f.writelines(data2)
   if roms_branch==1:
      f.writelines(data2)
      f.write('\n   Hadvection == ' + Hadvectionstr + '\n')
      f.write('   Vadvection == ' + Vadvectionstr + '\n\n')
      f.writelines(data3)

   for i,variable in enumerate(model.interior_state_variables,1):
      if i==NBT:
         s1='  '
      else:
         s1='\ '

      l1 = len(variable.output_name)
      if l1>lmax:
         vstr = 'BS' + str(i)
         l1 = len(vstr)
      else:
         vstr = variable.output_name

      if i==1:
         if isforced_bc[i-1]==1:
            f.write('\n     LBC(isTvar) == ' + bc_forced + s1 + '   ! idbio( 1), ' + vstr + '\n')
         else:
            f.write('\n     LBC(isTvar) == ' + bc_unforced + s1 + '   ! idbio( 1), ' + vstr + '\n')
      else:
         if isforced_bc[i-1]==1:
            f.write(' '*20 + bc_forced + s1 + '   ! idbio( ' + str(i) + '), ' + vstr + '\n')
         else:
            f.write(' '*20 + bc_unforced + s1 + '   ! idbio( ' + str(i) + '), ' + vstr + '\n')

   f.write('\n! Adjoint-based algorithms can have different lateral boundary\n')
   f.write('! conditions keywords.\n\n')
   f.write('ad_LBC(isTvar) ==   Gra     Clo     Gra     Gra\n\n')

   f.write('! Logical switches (TRUE/FALSE) to activate biological tracers point\n')
   f.write('! Sources/Sinks (like river runoff) and to specify which tracer variable\n')
   f.write('! to consider: [NBT,Ngrids] values are expected. See glossary below for\n')
   f.write('! details.\n\n')

   if default_LtracerSrc==1:
      f.write('  LtracerSrc == ' + LtracerSrc + '\n\n')
   else:
      f.write('  LtracerSrc == ' + str(NBT) + '*F\n\n')

   f.write('! Logical switches (TRUE/FALSE) to read and process biological tracer\n')
   f.write('! climatology fields: [NBT,Ngrids] values are expected. See glossary below\n')
   f.write('! for details.\n\n')

   f.write('  LtracerCLM == ' + str(NBT) + '*F\n\n')

   f.write('! Logical switches (TRUE/FALSE) to nudge the desired biological tracer\n')
   f.write('! climatology field. If not analytical climatology fields, users need to\n')
   f.write('! turn on the logical switches above to process the fields from the\n')
   f.write('! climatology NetCDF file that are needed for nudging; [NBT,Ngrids]\n')
   f.write('! values are expected. See glossary below for details.\n\n')

   f.write('  LnudgeTCLM == ' + str(NBT) + '*F\n\n')

   f.write('! Logical switches (TRUE/FALSE) to activate writing of biological fields\n')
   f.write('! into HISTORY output file: [1:NBT,Ngrids].\n\n')

   f.write('Hout(idTvar) == ' + str(NBT) + '*T    ! ..., NO3, ...           biological tracer\n')
   if NSAT>0:
      f.write('Hout(idSAvar) == ' + str(NSAT) + '*T   !                         surface-attached variables\n')
   if NBAT>0:
      f.write('Hout(idBAvar) == ' + str(NBAT) + '*T   !                         bottom-attached variables\n')

   if roms_branch==1:
      f.write('\n! Logical switches (TRUE/FALSE) to activate writing of biological fields\n')
      f.write('! into QUICKSAVE output file: [1:NBT,Ngrids].\n\n')

      f.write('Qout(idTvar) == ' + str(NBT) + '*F    ! ..., NO3, ...           biological tracer\n')
      if NSAT>0:
         f.write('Qout(idSAvar) == ' + str(NSAT) + '*F   !                         surface-attached variables\n')
      if NBAT>0:
         f.write('Qout(idBAvar) == ' + str(NBAT) + '*F   !                         bottom-attached variables\n')

   f.write('\n! Logical switches (TRUE/FALSE) to activate writing of time-averaged fields\n')
   f.write('! into AVERAGE output file: [1:NBT,Ngrids].\n\n')

   f.write('Aout(idTvar) == ' + str(NBT) + '*T    ! ..., NO3, ...           biological tracer\n')
   if NSAT>0:
      f.write('Aout(idSAvar) == ' + str(NSAT) + '*T   !                         surface-attached variables\n')
   if NBAT>0:
      f.write('Aout(idBAvar) == ' + str(NBAT) + '*T   !                         bottom-attached variables\n')

   f.write('\n! Logical switches (TRUE/FALSE) to activate writing of time-averaged,\n')
   f.write('! biological tracer diagnostic terms into DIAGNOSTIC output file:\n')
   f.write('! [1:NBT,Ngrids].\n\n')

   f.write('! Logical switches (TRUE/FALSE) to activate writing of time-averaged,\n')
   f.write('! biological processes diagnostics terms into DIAGNOSTIC output file [Ngrids].\n\n')

   if default_Dout==1:
      #Define default set of diagnostics to output
      #(application specific, can be modified at run-time)

      f.write('! O2 related\n')
      for i,v in enumerate(defDoutO2_2d,1):
         if len(v)>lmax:
            vstr = 'HD' + str(i)
         else:
            vstr = v
         f.write('Dout(i' + vstr + ') == T\n')
      for i,v in enumerate(defDoutO2_3d,1):
         if len(v)>lmax:
            vstr = 'BD' + str(i)
         else:
            vstr = v
         f.write('Dout(i' + vstr + ') == T\n')

      f.write('\n! DIC related\n')
      for i,v in enumerate(defDoutDIC_2d,1):
         if len(v)>lmax:
            vstr = 'HD' + str(i)
         else:
            vstr = v
         f.write('Dout(i' + vstr + ') == T\n')
      for i,v in enumerate(defDoutDIC_3d,1):
         if len(v)>lmax:
            vstr = 'BD' + str(i)
         else:
            vstr = v
         f.write('Dout(i' + vstr + ') == T\n')

      f.write('\n! TA related\n')
      for i,v in enumerate(defDoutTA_2d,1):
         if len(v)>lmax:
            vstr = 'HD' + str(i)
         else:
            vstr = v
         f.write('Dout(i' + vstr + ') == T\n')
      for i,v in enumerate(defDoutTA_3d,1):
         if len(v)>lmax:
            vstr = 'BD' + str(i)
         else:
            vstr = v
         f.write('Dout(i' + vstr + ') == T\n')

      f.write('\n! DOM related\n')
      for i,v in enumerate(defDoutDOM_2d,1):
         if len(v)>lmax:
            vstr = 'HD' + str(i)
         else:
            vstr = v
         f.write('Dout(i' + vstr + ') == T\n')
      for i,v in enumerate(defDoutDOM_3d,1):
         if len(v)>lmax:
            vstr = 'BD' + str(i)
         else:
            vstr = v
         f.write('Dout(i' + vstr + ') == T\n')

      f.write('\n! Aggregate variables\n')
      for i,v in enumerate(defDoutAgg_2d,1):
         if len(v)>lmax:
            vstr = 'HD' + str(i)
         else:
            vstr = v
         f.write('Dout(i' + vstr + ') == T\n')
      for i,v in enumerate(defDoutAgg_3d,1):
         if len(v)>lmax:
            vstr = 'BD' + str(i)
         else:
            vstr = v
         f.write('Dout(i' + vstr + ') == T\n')

      f.write('\n! Light related\n')
      for i,v in enumerate(defDoutLight_2d,1):
         if len(v)>lmax:
            vstr = 'HD' + str(i)
         else:
            vstr = v
         f.write('Dout(i' + vstr + ') == T\n')
      for i,v in enumerate(defDoutLight_3d,1):
         if len(v)>lmax:
            vstr = 'BD' + str(i)
         else:
            vstr = v
         f.write('Dout(i' + vstr + ') == T\n')

      f.write('\n! Other diagnostics\n')
      for i,v in enumerate(defDoutOther_2d,1):
         if len(v)>lmax:
            vstr = 'HD' + str(i)
         else:
            vstr = v
         f.write('Dout(i' + vstr + ') == T\n')
      for i,v in enumerate(defDoutOther_3d,1):
         if len(v)>lmax:
            vstr = 'BD' + str(i)
         else:
            vstr = v
         f.write('Dout(i' + vstr + ') == T\n')

   if roms_branch==0:
      f.writelines(data3)
   if roms_branch==1:
      f.writelines(data4)



## fabm.cdl
   os.system('python fabm_describe_model_pwa.py ' + yamlfile + ' > ' + outfile6)



# Useful python commands:
#
#model.printInformation() #To output all model information
#print dir(model) #To output all model attributes
#print model.findParameter('light/EPS0r',case_insensitive=True).long_name #Get the long name
