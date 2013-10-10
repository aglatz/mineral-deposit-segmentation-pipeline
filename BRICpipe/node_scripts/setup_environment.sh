#!/bin/bash

SLC_PFX=/media/LP3TBdisk/Andreas_PhD/Slicer3
FSL_PFX=/usr/local/fsl
FLIRT_PASS1="-dof 6 -cost normmi -searchcost normmi"
FLIRT_PASS2=""
AR_IDX="1"
WM_IDX="4"
GM_IDX="3"
CS_IDX="2"
BET_GRE="-f 0.3"
BET_T1W="-f 0.3 -B"
BET_T2W="-f 0.3 -B"
BET_FLAIR="-f 0.25 -B"

# Return environment variable name
eval "echo \$$1"
exit 0
