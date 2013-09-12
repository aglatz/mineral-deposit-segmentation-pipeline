#!/bin/bash

OUT_FILES=$1
DST_PATH=$2

rsync -azi ${OUT_FILES} ${DST_PATH} > /dev/zero

exit 0
