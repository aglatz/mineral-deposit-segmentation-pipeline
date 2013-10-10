#!/bin/bash

SLC_PFX=/usr/local/Slicer3/
FSL_PFX=/usr/local/fsl/4.1.9/
FLIRT_PASS1="-dof 6"
FLIRT_PASS2=""
AR_IDX="1"
WM_IDX="3"
GM_IDX="4"
CS_IDX="1"
BET_GRE="-f 0.4"
BET_T1W="-f 0.25 -B"
BET_T2W="-f 0.25 -B"
BET_FLAIR="-f 0.25 -B"

# Return environment variable name
eval "echo \$$1"
exit 0
