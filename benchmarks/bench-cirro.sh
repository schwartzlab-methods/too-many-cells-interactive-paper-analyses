#! /usr/bin/env bash

# Script for benchmarking performance of cirrocumulus again test dataset.
# NOTE: prerequeisite is image tagged `cirro-runner` built from Dockerfile.cirrocumulus

set -euo pipefail

start_time=$(date +"%T.%N")
cgroup_mem_root=$(grep cgroup /proc/mounts | grep memory | awk '{print $2}')/docker

results_file=_results-cirro-tmci-full.csv

container_name=cirro-runner

echo "type,dataset,max_mem,start,finish" > $results_file

get_max_mem () {

    max_mem="0"

    while [[ -n $(docker ps -q --filter name="${1}") ]]; do

        container_id=$(docker ps --no-trunc -q --filter name="${1}")

        sleep .1

        if [[ -f "$cgroup_mem_root"/"$container_id"/memory.max_usage_in_bytes ]]; then
            mem=$(cat "$cgroup_mem_root"/"$container_id"/memory.max_usage_in_bytes)
        fi


    done

    echo $mem
}

docker run --rm -i -d -p 5004:5004 --name ${container_name} --log-driver=journald -v $PWD/data/tm-full:/code/matrices:ro cirro-runner

# we'll kill this container manually
max_mem=$(get_max_mem ${container_name})

echo "cirro,$max_mem,$start_time," >> $results_file
