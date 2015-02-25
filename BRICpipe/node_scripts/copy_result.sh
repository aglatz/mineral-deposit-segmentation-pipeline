#!/bin/bash

DST_PATH=$1; shift
OUT_FILES="$@"

rsync -azi ${OUT_FILES} ${DST_PATH} > /dev/zero

exit 0
