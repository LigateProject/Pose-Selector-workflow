#!/bin/bash
set -e

DATA=

echo "Poses that could not be completed:"
poses=$(find ${DATA} -name hq-work)

for pose in ${poses}; do
  echo ${pose}
  # clean up
  rm -r ${pose}
done
