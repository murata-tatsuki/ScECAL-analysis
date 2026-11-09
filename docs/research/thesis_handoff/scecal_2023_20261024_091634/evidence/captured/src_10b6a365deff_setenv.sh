# Source this script to set up
#  GCC14.2.1 (gcc-toolset-14)
#  Python3.12.12 (python3.12 via dnf)
#  ROOT 6-36-08
#  GenFit 
#  Geant4 11.01.p03
#


drop_from_path_my()
{
   # Assert that we got enough arguments
   if test $# -ne 2 ; then
      echo "drop_from_path_my: needs 2 arguments"
      return 1
   fi

   p=$1
   drop=$2

   newpath=`echo $p | sed -e "s;:${drop}:;:;g" \
                          -e "s;:${drop};;g"   \
                          -e "s;${drop}:;;g"   \
                          -e "s;${drop};;g"`
}

export PATH=${PATH}:/meg/home/meg/local/bin
export MEG2LIB=/meg/home/meg/offline/meg2lib

# Astyle
if [ -n "${ASTYLE_PATH}" ] ; then
   if [ -n "${PATH}" ]; then
       drop_from_path_my $PATH $ASTYLE_PATH
       PATH=$newpath
   fi
fi
export ASTYLE_PATH="/meg/home/meg/cern/astyle3.1"
export PATH=${PATH}:${ASTYLE_PATH}


# ONNX Runtime
if [ -n "${ONNXSYS}" ] ; then
    old_onnxsys=${ONNXSYS}
    if [ -n "${LD_LIBRARY_PATH}" ]; then
        drop_from_path_my $LD_LIBRARY_PATH ${old_onnxsys}/lib64
        LD_LIBRARY_PATH=$newpath
    fi
    if [ -n "${ROOT_INCLUDE_PATH}" ]; then
	drop_from_path_my $ROOT_INCLUDE_PATH ${old_onnxsys}/include/onnxruntime/core/session
	ROOT_INCLUDE_PATH=$newpath
    fi
fi
export ONNXSYS=/meg/home/meg/cern/onnxruntime_v1.19.2_gcc1421_release
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ONNXSYS/lib64
if [ -n "$ROOT_INCLUDE_PATH" ] ; then
    export ROOT_INCLUDE_PATH=${ROOT_INCLUDE_PATH}:${ONNXSYS}/include/onnxruntime/core/session
else   
    export ROOT_INCLUDE_PATH=${ONNXSYS}/include/onnxruntime/core/session
fi


# GEANT4
if [ -n "${GEANT4SYS}" ] ; then
   if [ -n "${PATH}" ]; then
       drop_from_path_my $PATH $GEANT4SYS/bin
       PATH=$newpath
   fi
   if [ -n "${LD_LIBRARY_PATH}" ]; then
      drop_from_path_my $LD_LIBRARY_PATH $GEANT4SYS/lib64
      LD_LIBRARY_PATH=$newpath
   fi
fi
if [ -n "${G4WORKDIR}" ] ; then
   old_g4workdir=${G4WORKDIR}
   if [ -n "${PATH}" ]; then
      drop_from_path_my $PATH ${old_g4workdir}/bin/${G4SYSTEM}
      PATH=$newpath
   fi
fi
export GEANT4SYS=/meg/home/meg/cern/geant4.11.1.3_gcc1421_release
CURRENT_DIR=`pwd`
cd ${GEANT4SYS}/share/Geant4/geant4make
source ${GEANT4SYS}/share/Geant4/geant4make/geant4make.sh
cd ${CURRENT_DIR}
drop_from_path_my $LD_LIBRARY_PATH /usr/lib64
LD_LIBRARY_PATH=$newpath

# ROOT
source /meg/home/meg/cern/root-v6.36.08_gcc1421_python312_release/bin/thisroot.sh

if [ -n "${GENFIT2SYS}" ] ; then
    old_genfit2sys=${GENFIT2SYS}
    if [ -n "${LD_LIBRARY_PATH}" ]; then
	drop_from_path_my $LD_LIBRARY_PATH ${old_genfit2sys}/lib64
	drop_from_path_my $newpath ${old_genfit2sys}/lib
	LD_LIBRARY_PATH=$newpath
    fi
fi
export GENFIT2SYS=/meg/home/meg/cern/GenFit2.0.4_gcc1421_root6.36.08_release
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${GENFIT2SYS}/lib64

export SQLITESYS=/meg/home/meg/cern/sqlite3.7.17
export PATH=$SQLITESYS/bin:$PATH
#export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:$SQLITESYS/lib

source /opt/rh/gcc-toolset-14/enable
#source /opt/rh/rh-python38/enable
#export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH
