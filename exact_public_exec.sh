#!/bin/bash

for file in exact-public-instances/*
do
    echo "------"
    echo "running on $file"
    timeout 60 /usr/bin/time -f "exec time: %E (h:m:s)" ./target/release/sherby-pace-2024 $file
    exit_status=$?
    if [[ $exit_status -eq 124 ]]; then
        echo "TIME OUT"
    fi
done
