#!/bin/bash

SLC_PFX=/usr/local/Slicer3/
FSL_PFX=/usr/local/fsl/4.1.9/
#SLC_PFX=/ISIS/proc1/aglatz/Slicer3/
#FSL_PFX=/usr/local/fsl/
FLIRT_PASS1="-dof 6"
FLIRT_PASS2=""
#FLIRT_PASS2="-dof 6"
AR_IDX="1"
WM_IDX="3"
GM_IDX="4"
CS_IDX="1"
#AR_IDX="4"
#WM_IDX="2"
#GM_IDX="3"
#CS_IDX="4"
BET_GRE="-f 0.4"
BET_T1W="-f 0.25 -B"
BET_T2W="-f 0.25 -B"
BET_FLAIR="-f 0.25 -B"

# Return environment variable name
eval "echo \$$1"
exit 0
