#!/bin/bash
#
# svn $Id: build.sh 139 2008-01-10 00:17:29Z arango $
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::: John Wilkin :::
# Copyright (c) 2002-2008 The ROMS/TOMS Group                           :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.txt                                                :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::: Hernan G. Arango :::
#                                                                       :::
#                                                                       :::
# Script to compile an user application where the application-specific  :::
# files are kept separate from the ROMS source code.                    :::
#                                                                       :::
# Q: How/why does this script work?                                     :::
#                                                                       :::
# A: The ROMS makefile configures user-defined options with a set of    :::
#    flags such as ROMS_APPLICATION. Browse the makefile to see these.  :::
#    If an option in the makefile uses the syntax ?= in setting the     :::
#    default, this means that make will check whether an environment    :::
#    variable by that name is set in the shell that calls make. If so   :::
#    the environment variable value overrides the default (and the      :::
#    user need not maintain separate makefiles, or frequently edit      :::
#    the makefile, to run separate applications).                       :::
#                                                                       :::
# Usage:                                                                :::
#                                                                       :::
#    ./build.sh [options]                                               :::
#                                                                       :::
# Options:                                                              :::
#                                                                       :::
#    -j [N]      Compile in parallel using N CPUs                       :::
#                  omit argument for all available CPUs                 :::
#    -noclean    Do not clean already compiled objects                  :::
#                                                                       :::
# Notice that sometimes the parallel compilation fail to find MPI       :::
# include file "mpif.h".                                                :::
#                                                                       :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


parallel=1
clean=1

while [ $# -gt 0 ]
do
  case "$1" in
    -j )
      shift
      parallel=1
      test=`echo $1 | grep -P '^\d+$'`
      if [ "$test" != "" ]; then
        NCPUS="-j $1"
        shift
      else
        NCPUS="-j"
      fi
      ;;

    -noclean )
      shift
      clean=0
      ;;

    * )
      echo ""
      echo "$0 : Unknown option [ $1 ]"
      echo ""
      echo "Available Options:"
      echo ""
      echo "-j [N]      Compile in parallel using N CPUs"
      echo "              omit argument for all avaliable CPUs"
      echo "-noclean    Do not clean already compiled objects"
      echo ""
      exit 1
      ;;
  esac
done

# Set the CPP option defining the particular application. This will
# determine the name of the ".h" header file with the application
# CPP definitions.

export ROMS_APPLICATION=a20_v3_fabm_npzd_franks

# Set a local environmental variable to define the path to the directories
# where all this project's files are kept.

export MY_ROOT_DIR=/cluster/projects/nn9630k/A20/

export MY_PROJECT_DIR=${MY_ROOT_DIR}/Run3
export MY_COMPILE_DIR=${MY_PROJECT_DIR}
echo "Compiling ",$ROMS_APPLICATION
# The path to the user's local current ROMS source code.
#
# If using svn locally, this would be the user's Working Copy Path (WCPATH).
# Note that one advantage of maintaining your source code locally with svn
# is that when working simultaneously on multiple machines (e.g. a local
# workstation, a local cluster and a remote supercomputer) you can checkout
# the latest release and always get an up-to-date customized source on each
# machine. This script is designed to more easily allow for differing paths
# to the code and inputs on differing machines.

#
export MY_ROMS_SRC=${MY_ROOT_DIR}/ROMS_v37_KATE

###Inserted PWA 01/07/2014
NEW_SRC_DIR=${MY_ROMS_SRC}/ROMS/Nonlinear/Biology


# # KHC - 20110209
# # Check if we have any modified source files

if [ -s Src_modify ]; then
    cd Src_modify/
#    gotModifiedSource=`ls *.F *.h *.mk`
    gotModifiedSource=`ls *.F *.h *.mk makefile`
###Modified PWA 04/03/2015
    cd ..
fi

# Replace the original files with the modifications
if [ "$gotModifiedSource" != "" ]; then

    # Copy locally modified source to main ROMS directory
    for ModSrc_modify in $gotModifiedSource; do

        # Check where original resides
###PWA Modified 04/03/2015
        #	origFile=`find $MY_ROMS_SRC -name $ModSrc_modify
        if [ "$ModSrc_modify" == "makefile" ]; then
        # For the makefile, consider only the makefile in $MY_ROMS_SRC (PWA)
            echo "Copying makefile"
            echo $ModSrc_modify
            origFile=$MY_ROMS_SRC/makefile
        else
            origFile=`find $MY_ROMS_SRC -name $ModSrc_modify`
            if [ "$origFile" != "" ]; then
                echo "Found original file:"
                echo $origFile
            fi
        fi
###PWA Modified 04/03/2015

        if [ -f "$origFile" ]; then

            # Moving original and copying user-modifed source code
            # first checking if the original already exists with
            # the .orig extension
            if [ ! -f "$origFile.orig" ]; then
                mv $origFile $origFile.orig
                echo "Moving $origFile to $origFile.orig"
            fi

            # Copying from local source directory to repository
            cp Src_modify/$ModSrc_modify $origFile
            echo "Copying Src_modify/$ModSrc_modify to $origFile"

            if [ ! -f USER_MODIFIED_CODE_IN_REPO ]; then

                # Touch file to notify that user modified code has been
                # placed in the repository
                touch USER_MODIFIED_CODE_IN_REPO

            fi
        else

	    ##PWA changes

        # No such file in repository, proceed to copy in new file (PWA)
            cp Src_modify/$ModSrc_modify $NEW_SRC_DIR/.
            echo "Copying Src_modify/$ModSrc_modify to $NEW_SRC_DIR"

            if [ ! -f USER_MODIFIED_CODE_IN_REPO ]; then

                # Touch file to notify that user modified code has been
                # placed in the repository
                touch USER_MODIFIED_CODE_IN_REPO

            fi
	    ##End of PWA changes
        fi
    done
fi
# Removing user modified source code in repository
# KHC - 20110209
# NMK - 2013
rollback() {
    cd $MY_COMPILE_DIR

    if [ -f USER_MODIFIED_CODE_IN_REPO ]; then

    # Find source code files with ".orig"-ending and
    # remove ending
        filelist=`find "$MY_ROMS_SRC" -name *.orig`

        if [ "$filelist" != "" ]; then

            for oldFileName in $filelist; do

            # extract basename
            newFileName=`basename $oldFileName .orig`
            fileDirectory=`dirname $oldFileName`
            mv $oldFileName  $fileDirectory/$newFileName

            echo "Moved $oldFileName  to $fileDirectory/$newFileName"

            done

        else # Empty filelist, no such files in repository

            echo "Did not find any .orig-files in the repository, empty file deleted"

        fi

        # Remove empty file

        rm -f USER_MODIFIED_CODE_IN_REPO

    fi
}
trap 'rollback; exit 99' 0


# Set tunable CPP options.
#
# Sometimes it is desirable to activate one or more CPP options to run
# different variants of the same application without modifying its header
# file. If this is the case, specify each options here using the -D syntax.
# Notice also that you need to use shell's quoting syntax to enclose the
# definition.  Both single or double quotes works. For example, to write
# time-averaged fields set:
#
#    export MY_CPP_FLAGS="-DAVERAGES"

#export MY_CPP_FLAGS="-D"

# Other user defined environmental variables. See the ROMS makefile for
# details on other options the user might want to set here.

 export USE_MPI=on
 export USE_MPIF90=on
 export FORT=ifort

#export USE_OpenMP=on

# export USE_DEBUG=on
 export USE_LARGE=on
 export USE_NETCDF4=on

# Activate Data Access Protocol (like OPeNDAP) support for input
# NetCDF files.  This is only possible for NetCDF library version
# 4.1.1 or higher. Also, add the path of the "curl-config" script
# that includes all linking libraries for DAP support.

if [ -n "${USE_NETCDF4:+1}" ]; then
    export USE_DAP=on
    export PATH=/usr/bin:$PATH
fi

# There are several MPI libraries out there. The user can select here the
# appropriate "mpif90" script to compile, provided that the makefile
# macro file (say, Linux-pgi.mk) in the Compilers directory has:
#
#         FC := mpif90
#
# "mpif90" defined without any path. Recall that you still need to use the
# appropriate "mpirun" to execute. Also notice that the path where the
# MPI library is installed is computer dependent.

if [ -n "${USE_MPIF90:+1}" ]; then
  case "$FORT" in
    ifort )
	  echo "exporting"
#      export PATH=/software/intel/impi/4.0.0.027/bin:$PATH
#      export PATH=/software/intel/impi/3.2.0.011/bin:$PATH
#      export PATH=/opt/intelsoft/mpich/bin:$PATH
#     export PATH=/opt/intelsoft/mpich2/bin:$PATH
#     export PATH=/opt/intelsoft/openmpi/bin:$PATH
      ;;

    pgi )
      export PATH=/opt/pgisoft/mpich/bin:$PATH
#     export PATH=/opt/pgisoft/mpich2/bin:$PATH
#     export PATH=/opt/pgisoft/openmpi/bin:$PATH
      ;;

    gfortran )
#     export PATH=/opt/gfortransoft/mpich2/bin:$PATH
      #export PATH=/disk1/openmpi/bin:$PATH
      ;;

    g95 )
#     export PATH=/opt/g95soft/mpich2/bin:$PATH
      export PATH=/opt/g95soft/openmpi/bin:$PATH
      ;;

  esac
fi

# The path of the libraries required by ROMS can be set here using
# environmental variables which take precedence to the values
# specified in the makefile macro definitions file (Compilers/*.mk).
# If so desired, uncomment the local USE_MY_LIBS definition below
# and edit the paths to your values. For most applications, only
# the location of the NetCDF library (NETCDF_LIBDIR) and include
# directorry (NETCDF_INCDIR) are needed!
#
# Notice that when the USE_NETCDF4 macro is activated, we need a
# serial and parallel version of the NetCDF-4/HDF5 library. The
# parallel library uses parallel I/O through MPI-I/O so we need
# compile also with the MPI library. This is fine in ROMS
# distributed-memory applications.  However, in serial or
# shared-memory ROMS applications we need to use the serial
# version of the NetCDF-4/HDF5 to avoid conflicts with the
# compiler. Recall also that the MPI library comes in several
# flavors: MPICH, MPICH2, OpenMPI.

#export            USE_MY_LIBS=on

if [ -n "${USE_MY_LIBS:+1}" ]; then
  case "$FORT" in
    ifort )
      export           ESMF_DIR=/opt/intelsoft/esmf-3.1.0
      export            ESMF_OS=Linux
      export      ESMF_COMPILER=ifort
      export          ESMF_BOPT=O
      export           ESMF_ABI=64
      export          ESMF_COMM=mpich
      export          ESMF_SITE=default
      export         MCT_INCDIR=/opt/intelsoft/mct/include
      export         MCT_LIBDIR=/opt/intelsoft/mct/lib

      export      ARPACK_LIBDIR=/opt/intelsoft/serial/ARPACK
      if [ -n "${USE_MPI:+1}" ]; then
        export   PARPACK_LIBDIR=/opt/intelsoft/mpich/PARPACK
#       export   PARPACK_LIBDIR=/opt/intelsoft/mpich2/PARPACK
#       export   PARPACK_LIBDIR=/opt/intelsoft/openmpi/PARPACK
      fi

      if [ -n "${USE_NETCDF4:+1}" ]; then
        if [ -n "${USE_MPI:+1}" ]; then
          export  NETCDF_INCDIR=/opt/intelsoft/mpich/netcdf4/include
          export  NETCDF_LIBDIR=/opt/intelsoft/mpich/netcdf4/lib
          export    HDF5_LIBDIR=/opt/intelsoft/mpich/hdf5/lib

#         export  NETCDF_INCDIR=/opt/intelsoft/mpich2/netcdf4/include
#         export  NETCDF_LIBDIR=/opt/intelsoft/mpich2/netcdf4/lib
#         export    HDF5_LIBDIR=/opt/intelsoft/mpich2/hdf5/lib

#         export  NETCDF_INCDIR=/opt/intelsoft/openmpi/netcdf4/include
#         export  NETCDF_LIBDIR=/opt/intelsoft/openmpi/netcdf4/lib
#         export    HDF5_LIBDIR=/opt/intelsoft/openmpi/hdf5/lib
        else
          export  NETCDF_INCDIR=/opt/intelsoft/serial/netcdf4/include
          export  NETCDF_LIBDIR=/opt/intelsoft/serial/netcdf4/lib
          export    HDF5_LIBDIR=/opt/intelsoft/serial/hdf5/lib
        fi
      else
        export    NETCDF_INCDIR=/opt/intelsoft/serial/netcdf3/include
        export    NETCDF_LIBDIR=/opt/intelsoft/serial/netcdf3/lib
      fi
      ;;

    pgi )
      export           ESMF_DIR=/opt/pgisoft/esmf-3.1.0
      export            ESMF_OS=Linux
      export      ESMF_COMPILER=pgi
      export          ESMF_BOPT=O
      export           ESMF_ABI=64
      export          ESMF_COMM=mpich
      export          ESMF_SITE=default
      export         MCT_INCDIR=/opt/pgisoft/mct/include
      export         MCT_LIBDIR=/opt/pgisoft/mct/lib

      export      ARPACK_LIBDIR=/opt/pgisoft/serial/ARPACK
      if [ -n "${USE_MPI:+1}" ]; then
        export   PARPACK_LIBDIR=/opt/pgisoft/mpich/PARPACK
#       export   PARPACK_LIBDIR=/opt/pgisoft/mpich2/PARPACK
#       export   PARPACK_LIBDIR=/opt/pgisoft/openmpi/PARPACK
      fi

      if [ -n "${USE_NETCDF4:+1}" ]; then
        if [ -n "${USE_MPI:+1}" ]; then
          export  NETCDF_INCDIR=/opt/pgisoft/mpich/netcdf4/include
          export  NETCDF_LIBDIR=/opt/pgisoft/mpich/netcdf4/lib
          export    HDF5_LIBDIR=/opt/pgisoft/mpich/hdf5/lib

#         export  NETCDF_INCDIR=/opt/pgisoft/mpich2/netcdf4/include
#         export  NETCDF_LIBDIR=/opt/pgisoft/mpich2/netcdf4/lib
#         export    HDF5_LIBDIR=/opt/pgisoft/mpich2/hdf5/lib

#         export  NETCDF_INCDIR=/opt/pgisoft/openmpi/netcdf4/include
#         export  NETCDF_LIBDIR=/opt/pgisoft/openmpi/netcdf4/lib
#         export    HDF5_LIBDIR=/opt/pgisoft/openmpi/hdf5/lib
        else
          export  NETCDF_INCDIR=/opt/pgisoft/serial/netcdf4/include
          export  NETCDF_LIBDIR=/opt/pgisoft/serial/netcdf4/lib
          export    HDF5_LIBDIR=/opt/pgisoft/serial/hdf5/lib
        fi
      else
        export    NETCDF_INCDIR=/opt/pgisoft/serial/netcdf3/include
        export    NETCDF_LIBDIR=/opt/pgisoft/serial/netcdf3/lib
      fi
      ;;

    gfortran )
      export         MCT_INCDIR=/opt/gfortransoft/mct/include
      export         MCT_LIBDIR=/opt/gfortransoft/mct/lib

      export      ARPACK_LIBDIR=/disk1/ROMS/Lib/ARPACK
      if [ -n "${USE_MPI:+1}" ]; then
#       export   PARPACK_LIBDIR=/opt/gfortransoft/mpich/PARPACK
#       export   PARPACK_LIBDIR=/opt/gfortransoft/mpich2/PARPACK
        export   PARPACK_LIBDIR=/opt/gfortransoft/openmpi/PARPACK
      fi

      if [ -n "${USE_NETCDF4:+1}" ]; then
        if [ -n "${USE_MPI:+1}" ]; then
#         export  NETCDF_INCDIR=/opt/gfortransoft/mpich/netcdf4/include
#         export  NETCDF_LIBDIR=/opt/gfortransoft/mpich/netcdf4/lib
#         export    HDF5_LIBDIR=/opt/gfortransoft/mpich/hdf5/lib

#         export  NETCDF_INCDIR=/opt/gfortransoft/mpich2/netcdf4/include
#         export  NETCDF_LIBDIR=/opt/gfortransoft/mpich2/netcdf4/lib
#         export    HDF5_LIBDIR=/opt/gfortransoft/mpich2/hdf5/lib

          export  NETCDF_INCDIR=//include
          export  NETCDF_LIBDIR=/opt/gfortransoft/openmpi/netcdf4/lib
          export    HDF5_LIBDIR=/opt/gfortransoft/openmpi/hdf5/lib
        else
          export  NETCDF_INCDIR=/opt/gfortransoft/serial/netcdf4/include
          export  NETCDF_LIBDIR=/opt/gfortransoft/serial/netcdf4/lib
          export    HDF5_LIBDIR=/opt/gfortransoft/serial/hdf5/lib
        fi
      else
          export  NETCDF_INCDIR=/opt/gfortransoft/serial/netcdf3/include
          export  NETCDF_LIBDIR=/opt/gfortransoft/serial/netcdf3/lib
      fi
      ;;

    g95 )
      export         MCT_INCDIR=/opt/g95soft/mct/include
      export         MCT_LIBDIR=/opt/g95soft/mct/lib

      export      ARPACK_LIBDIR=/opt/g95soft/serial/ARPACK
      if [ -n "${USE_MPI:+1}" ]; then
#       export   PARPACK_LIBDIR=/opt/g95soft/mpich/PARPACK
#       export   PARPACK_LIBDIR=/opt/g95soft/mpich2/PARPACK
        export   PARPACK_LIBDIR=/opt/g95soft/openmpi/PARPACK
      fi

      if [ -n "${USE_NETCDF4:+1}" ]; then
        if [ -n "${USE_MPI:+1}" ]; then
#         export  NETCDF_INCDIR=/opt/g95soft/mpich2/netcdf4/include
#         export  NETCDF_LIBDIR=/opt/g95soft/mpich2/netcdf4/lib
#         export    HDF5_LIBDIR=/opt/g95soft/mpich2/hdf5/lib

          export  NETCDF_INCDIR=/opt/g95soft/openmpi/netcdf4/include
          export  NETCDF_LIBDIR=/opt/g95soft/openmpi/netcdf4/lib
          export    HDF5_LIBDIR=/opt/g95soft/openmpi/hdf5/lib
        else
          export  NETCDF_INCDIR=/opt/g95soft/serial/netcdf4/include
          export  NETCDF_LIBDIR=/opt/g95soft/serial/netcdf4/lib
          export    HDF5_LIBDIR=/opt/g95soft/serial/hdf5/lib
        fi
      else
        export    NETCDF_INCDIR=/opt/g95soft/serial/netcdf3/include
        export    NETCDF_LIBDIR=/opt/g95soft/serial/netcdf3/lib
      fi
      ;;

  esac
fi

# The rest of this script sets the path to the users header file and
# analytical source files, if any. See the templates in User/Functionals.

 export MY_HEADER_DIR=${MY_PROJECT_DIR}/Include

 export MY_ANALYTICAL_DIR=${MY_PROJECT_DIR}

# Put the binary to execute in the following directory.

 export BINDIR=${MY_PROJECT_DIR}

# Put the f90 files in a project specific Build directory to avoid conflict
# with other projects.

 export SCRATCH_DIR=${MY_PROJECT_DIR}/Build


###PWA Inserted 04/03/2015
 export FABM_INCDIR=/cluster/home/pwa/local/fabm/roms/include
 export FABM_LIBDIR=/cluster/home/pwa/local/fabm/roms/lib
#This defines the location of the FABM *.mod files and library,
#for use in the makefile if needed.
###PWA Inserted 04/03/2015


# Go to the users source directory to compile. The options set above will
# pick up the application-specific code from the appropriate place.

 cd ${MY_ROMS_SRC}

# Stop if activating both MPI and OpenMP at the same time.

#if [ ${?USE_MPI} & ${?USE_OpenMP} ]; then
#  echo "You cannot activate USE_MPI and USE_OpenMP at the same time!"
#  exit 1
#fi

# Remove build directory.

if [ $clean -eq 1 ]; then
  make clean
fi

# Compile (the binary will go to BINDIR set above).

if [ $parallel -eq 1 ]; then
  make $NCPUS
else
  make
fi
