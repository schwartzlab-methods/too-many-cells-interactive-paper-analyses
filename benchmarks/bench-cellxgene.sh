#! /usr/bin/env bash

set -euo pipefail

# This script assists in comparing performance of cellxgene and tmci.
# Both processes begin with a MEX matrix directory and end with a running webapp w/ visualizations and gene expression overlay
# NOTE: the input for tmci includes cluster_tree.json and labels.csv, which have been generated for the test matrix prior to running this program.
# NOTE: manually update variables below to change dataset or run type
# NOTE: should be run in container w/ image tagged `cellxgene-runner` built w/ Dockerfile.cellxgene

start_time=$(date +"%T.%N")
cgroup_mem_root=$(grep cgroup /proc/mounts | grep memory | awk '{print $2}')/docker

results_file=_results-cellxgene.csv
tmci_db_container_name=postgres
tmci_container_name=tmci

run_type=cellxgene

cellxgene_container_name=cellxgene-runner

dataset=full

tmci_webserver_container_name=tmci_node_server

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

get_max_mem_dual () {

    max_mem="0"

    while [[ -n $(docker ps -q --filter name="${1}") ]] && [[ -n $(docker ps -q --filter name="${2}") ]]; do

        container_ids=$(docker ps --no-trunc -q --filter name="${1}" --filter name="${2}" )

        total_mem="0"

        while IFS= read -r line; do
            if [[ -f "$cgroup_mem_root"/"$line"/memory.max_usage_in_bytes ]]; then
                cmem=$(cat "$cgroup_mem_root"/"$line"/memory.max_usage_in_bytes)
                total_mem=$(($total_mem + $cmem))
            fi
        done <<< "$container_ids"

        if [[ $total_mem -gt $max_mem ]]; then
            max_mem=$total_mem
        fi

        sleep .1

    done

    echo $max_mem
}

# This benchmarking process is pretty manual. Because each benchmarked code path starts a server in the foreground, there is no exit time that can be easily recorded by the script.
# Instead, we need to inspect journald to get the time stamp of the "server started" message after manually killing the process. That's why there's no loop here.
# e.g., journalctl -f CONTAINER_NAME=cellxgene-runner

if [[ $run_type == 'tmci' ]]; then

    max_mem="0"

    # for tmci, we are just benchmarking the run-time of the start-and-load script with the matrix and pre-generated labels and tree.

    # start postgres in the background
    docker-compose -f docker-compose.yaml run --rm --name ${tmci_db_container_name} -d postgres

    #load data
    docker-compose -f docker-compose.yaml run --rm -d --name ${tmci_container_name} \
        -v "$(readlink -f ./data/tm-full):/usr/data/matrices:ro" \
        node --init

    max_mem_load=$(get_max_mem_dual ${tmci_db_container_name} ${tmci_container_name})


    docker-compose -f docker-compose.yaml run -d --name ${tmci_webserver_container_name} --rm  \
    -v "$(readlink -f ./data/prior/tm-full/cluster_tree.json):/usr/app/static/files/cluster_tree.json:ro" \
    -v "$(readlink -f ./data/prior/tm-full/labels.csv):/usr/app/static/files/labels.csv:ro" \
    node --prod

    # we'll kill this container manually
    max_mem_run=$(get_max_mem_dual ${tmci_db_container_name} ${tmci_webserver_container_name})

    if [ $max_mem_run -gt  $max_mem_load ]; then
        max_mem=$max_mem_run
    else
        max_mem=$max_mem_load
    fi

    echo "tmci,$dataset,$max_mem,$start_time," >> $results_file

    docker stop ${tmci_db_container_name}
    docker rm ${tmci_db_container_name}
    docker rm ${tmci_container_name}
    docker rm ${tmci_webserver_container_name}

else
    docker run -d -i --name ${cellxgene_container_name} --rm --log-driver=journald -v "${PWD}"/data/tm-full:/code/matrices -p 5004:5004 cellxgene-runner:latest

    #we'll kill this container manually
    max_mem=$(get_max_mem ${cellxgene_container_name})

    echo "cellxgene,$dataset,$max_mem,$start_time," >> $results_file

fi



