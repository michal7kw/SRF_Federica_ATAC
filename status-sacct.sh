#!/bin/bash
jobid="$1"

if [[ "$jobid" == "" ]]; then
    echo "No jobid given" >&2
    exit 1
fi

output=$(sacct -j "$jobid" --format=State --noheader | head -n 1 | awk '{print $1}')

if [[ "$output" == "COMPLETED" ]]; then
    echo "success"
elif [[ "$output" == "RUNNING" ]]; then
    echo "running"
elif [[ "$output" == "PENDING" ]]; then
    echo "running"
else
    echo "failed"
fi
