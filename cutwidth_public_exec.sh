#!/bin/bash

arr=()

for i in {1..125}; 
do
    echo "------"
    echo "running on file $i"
    systemd-run --send-sighup --scope -p MemoryLimit=8000M timeout 1800 /usr/bin/time -f "exec time: %E (h:m:s)" ./target/release/sherby-cutwidth < cutwidth-public//$i.gr > cutwidth-public//$i.sol
    exit_status=$?
    if [[ $exit_status -eq 124 ]]; then
        echo "TIME OUT"
        arr+=('\033[0;31m#\033[0m]')
    fi
    if [[ $exit_status -eq 143 ]]; then
        echo "OUT OF MEMORY"
        arr+=('\033[0;31m#\033[0m]')
    fi
    if [[ $exit_status -eq 0 ]]; then
        echo "SUCCESS"
        arr+=('\032[0;31m#\033[0m]')
    fi
done

for c in "${arr[@]}"
do
    echo $c
done
