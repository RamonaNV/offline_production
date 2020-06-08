#!/bin/bash

# Add decreasing priorities to an existing DAG, ensuring
# that jobs are scheduled in the order they appear.

if [[ -z "$1" ]]; then
  echo "Usage: $0 foo.dag >> foo.dag"
fi

i=$(((1<<32 - 1)))

while read line; do
  fields=( $line )
  if [[ ${fields[0]} == "JOB" ]]; then
    i=$(($i-1))
    echo PRIORITY ${fields[1]} $i
  fi
done < "$1"
