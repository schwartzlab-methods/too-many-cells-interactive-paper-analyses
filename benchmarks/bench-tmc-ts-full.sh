#! /usr/bin/env bash

set -euo pipefail

result_file=_results-tmc-th.csv
tmc_container_name=tmc
tmci_container_name=tmci_node
tmci_db_container_name=postgres
prior=tabula-sapiens

echo "type,start,finish,RAM(b)," > $result_file

five_features=("Cd4" "Apoe" "Ahr" "Myc" "Serp1")

cgroup_mem_root=$(grep cgroup /proc/mounts | grep memory | awk '{print $2}')/docker

get_max_mem () {

    max_mem="0"
    mem="0"

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
                cmem=$(cat "$cgroup_mem_root"/"$line"/memory.max_usage_in_bytes || "0")
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

# for i in {1..10}; do
#     # we run the process in the background and watch for exit on mem usage loop and remove manually after benchmarking time
#     docker run -d -v `pwd`:`pwd` -w `pwd` \
#         --user 1000 --name ${tmc_container_name} \
#         gregoryschwartz/too-many-cells:3.0.1.0 make-tree \
#         --prior ./$prior/ \
#         --labels-file ./$prior/labels.csv \
#         --dendrogram-output newtree.svg \
#         --output ./data

#     max_mem=$(get_max_mem ${tmc_container_name})

#     start_finish=$(docker container inspect ${tmc_container_name} | jq .[0].State.StartedAt,.[0].State.FinishedAt | tr '\n' ',')

#     docker rm ${tmc_container_name}

#     echo -e "tmc,basic,${start_finish}${max_mem}," >> $result_file

# done

# for i in {1..10}; do
#     docker run -d -v `pwd`:`pwd` -w `pwd` \
#         --user 1000 --name ${tmc_container_name} \
#         gregoryschwartz/too-many-cells:3.0.1.0 make-tree \
#         --prior ./$prior/ \
#         --min-size 1000 \
#         --labels-file ./$prior/labels.csv \
#         --dendrogram-output newtree.svg \
#         --output ./data

#     max_mem=$(get_max_mem ${tmc_container_name})

#     start_finish=$(docker container inspect ${tmc_container_name} | jq .[0].State.StartedAt,.[0].State.FinishedAt | tr '\n' ',')

#     docker rm ${tmc_container_name}

#     echo -e "tmc,prune,${start_finish}${max_mem}," >> $result_file

# done


# these will run out of memory
# for i in {1..10}; do
#     docker run -d -v `pwd`:`pwd` -w `pwd` \
#         --user 1000 --name ${tmc_container_name} \
#         gregoryschwartz/too-many-cells:3.0.1.0 make-tree \
#         --matrix-path ./matrix-10x \
#         --labels-file ./$prior/labels.csv \
#         --prior ./$prior \
#         --dendrogram-output newtree.svg \
#         --output ./data

#     max_mem=$(get_max_mem ${tmc_container_name})

#     start_finish=$(docker container inspect ${tmc_container_name} | jq .[0].State.StartedAt,.[0].State.FinishedAt | tr '\n' ',')

#     docker rm ${tmc_container_name}

#     echo -e "tmc,with-matrix,${start_finish}${max_mem}," >> $result_file

# done


# for i in {1..10}; do

#     docker compose -f docker-compose.yaml \
#         run -d --entrypoint="" \
#         --name ${tmci_container_name} \
#         --no-deps \
#         -v "$(readlink -f ./$prior/labels.csv):/tmp/labels.csv:ro" \
#         -v "$(readlink -f ./$prior/cluster_tree.json):/tmp/cluster_tree.json:ro" \
#         -v "$(readlink -f ./blank-config.json):/tmp/config.json:ro" \
#         -v "$(readlink -f /tmp):/tmp/results" \
#         node node dist/exportTree.js \
#         --labelPath /tmp/labels.csv \
#         --treePath /tmp/cluster_tree.json \
#         --outPath /tmp/results/export667.svg \
#         --configPath /tmp/config.json

#     max_mem=$(get_max_mem ${tmci_container_name})

#     docker stop ${tmci_container_name}
#     start_finish=$(docker container inspect ${tmci_container_name} | jq .[0].State.StartedAt,.[0].State.FinishedAt | tr '\n' ',')
#     docker rm $tmci_container_name

#     echo -e "tmci,basic,${start_finish}${max_mem}," >> $result_file

# done

# for i in {1..10}; do

#     docker compose -f docker-compose.yaml \
#         run -d --entrypoint="" \
#         --name ${tmci_container_name} \
#         --no-deps \
#         -v "$(readlink -f ./$prior/labels.csv):/tmp/labels.csv:ro" \
#         -v "$(readlink -f ./$prior/cluster_tree.json):/tmp/cluster_tree.json:ro" \
#         -v "$(readlink -f ./blank-config.json):/tmp/config.json:ro" \
#         -v "$(readlink -f /tmp):/tmp/results" \
#         node node dist/exportTree.js \
#         --labelPath /tmp/labels.csv \
#         --treePath /tmp/cluster_tree.json \
#         --outPath /tmp/results/export667.svg \
#         --configPath /tmp/config.json

#     max_mem=$(get_max_mem ${tmci_container_name})

#     docker stop ${tmci_container_name}
#     start_finish=$(docker container inspect ${tmci_container_name} | jq .[0].State.StartedAt,.[0].State.FinishedAt | tr '\n' ',')

#     docker rm $tmci_container_name

#     echo -e "tmci,prune,${start_finish}${max_mem}," >> $result_file

# done

for i in {1..3}; do

    # start postgres in the background
    docker compose -f docker-compose.yaml run --name ${tmci_db_container_name} -d postgres

    #load data
    docker compose -f docker-compose.yaml run -d --name ${tmci_container_name} \
        -v "$(readlink -f ./matrix-10x):/usr/data/matrices:ro" \
        node --init

    echo "starting import"
    date

    max_mem_load=$(get_max_mem_dual ${tmci_db_container_name} ${tmci_container_name})

    echo "finished import"
    date

    docker rm ${tmci_container_name}

    # now run export
    docker compose -f docker-compose.yaml \
        run -d --entrypoint="" \
        --name ${tmci_container_name} \
        -v "$(readlink -f ./$prior/labels.csv):/tmp/labels.csv:ro" \
        -v "$(readlink -f ./$prior/cluster_tree.json):/tmp/cluster_tree.json:ro" \
        -v "$(readlink -f ./blank-config.json):/tmp/config.json:ro" \
        -v "$(readlink -f /tmp):/tmp/results" \
        node node dist/exportTree.js \
        --labelPath /tmp/labels.csv \
        --treePath /tmp/cluster_tree.json \
        --outPath /tmp/results/export667.svg \
        --configPath /tmp/config.json

    max_mem_export=$(get_max_mem_dual ${tmci_db_container_name} ${tmci_container_name})

    docker stop ${tmci_db_container_name} ${tmci_container_name}

    start_finish=$(docker container inspect ${tmci_db_container_name} | jq .[0].State.StartedAt,.[0].State.FinishedAt | tr '\n' ',')

    docker rm ${tmci_db_container_name} ${tmci_container_name} #? seems like down isn't working well w/ aliases?

    docker compose down -v

    if [[ $max_mem_load -gt $max_mem_export ]]; then
        max_mem=$max_mem_load
    else
        max_mem=$max_mem_export
    fi

    echo -e "tmci,with-matrix,${start_finish}${max_mem}," >> $result_file

done

# for i in {1..5}; do

#     # start postgres in the background
#     docker compose -f docker-compose.yaml run --name ${tmci_db_container_name} -d postgres

#     #load data
#     docker compose -f docker-compose.yaml run -d --name ${tmci_container_name} \
#         -v "$(readlink -f ./matrix-10x/):/usr/data/matrices:ro" \
#         node --init

#     echo "starting import"

#     max_mem_load=$(get_max_mem_dual ${tmci_db_container_name} ${tmci_container_name})

#     echo "import done"

#     docker rm ${tmci_container_name}

#     # now run export
#     docker compose -f docker-compose.yaml \
#         run -d --entrypoint="" \
#         --name ${tmci_container_name} \
#         -v "$(readlink -f ./$prior/labels.csv):/tmp/labels.csv:ro" \
#         -v "$(readlink -f ./$prior/cluster_tree.json):/tmp/cluster_tree.json:ro" \
#         -v "$(readlink -f ./feature-config.json):/tmp/config.json:ro" \
#         -v "$(readlink -f ./tmp):/tmp/results" \
#         node node dist/exportTree.js \
#         --labelPath /tmp/labels.csv \
#         --treePath /tmp/cluster_tree.json \
#         --outPath /tmp/results/feature-tree-tmci.svg \
#         --configPath /tmp/config.json

#     max_mem_export=$(get_max_mem_dual ${tmci_db_container_name} ${tmci_container_name})

#     docker stop ${tmci_db_container_name}

#     start_finish=$(docker container inspect ${tmci_db_container_name} | jq .[0].State.StartedAt,.[0].State.FinishedAt | tr '\n' ',')

#     docker rm ${tmci_db_container_name}  ${tmci_container_name}

#     docker compose -f docker-compose.yaml down -v

#     if [[ $max_mem_load -gt $max_mem_export ]]; then
#         max_mem=$max_mem_load
#     else
#         max_mem=$max_mem_export
#     fi

#     echo -e "tmci,with-feature,${start_finish}${max_mem}," >> $result_file

# done


# will run out of memory
# for i in {1..10}; do

#     docker run -d -v `pwd`:`pwd` -w `pwd` \
#         --user 1000 --name ${tmc_container_name} \
#         gregoryschwartz/too-many-cells:3.0.1.0 make-tree \
#         --prior ./$prior \
#         --matrix-path ./matrix-10x \
#         --feature-column 2 \
#         --draw-leaf "DrawItem (DrawContinuous [\"Cd4\"])" \
#         --labels-file ./$prior/labels.csv \
#         --dendrogram-output feature-tree-tmc.svg \
#         --output /tmp

#     max_mem=$(get_max_mem ${tmc_container_name})

#     start_finish=$(docker container inspect ${tmc_container_name} | jq .[0].State.StartedAt,.[0].State.FinishedAt | tr '\n' ',')

#     docker rm ${tmc_container_name}

#     echo -e "tmc,with-feature,${start_finish}${max_mem}," >> $result_file

# done



for i in {1..5}; do

    max_mem="0"

     # start postgres in the background
    docker compose -f docker-compose.yaml run --name ${tmci_db_container_name} -d postgres

    #load data
    docker compose -f docker-compose.yaml run -d --name ${tmci_container_name} \
        -v "$(readlink -f ./matrix-10x/):/usr/data/matrices:ro" \
        node --init

    max_mem_load=$(get_max_mem_dual ${tmci_db_container_name} ${tmci_container_name})

    docker rm ${tmci_container_name}

    sleep .1

    for feature in "${five_features[@]}"; do

        jq --arg feature "$feature" 'setpath(["features"];[$feature])' ./feature-config.json > ./feature-config-t.json

        docker compose -f docker-compose.yaml \
            run -d --entrypoint="" \
            --name ${tmci_container_name} \
            -v "$(readlink -f ./$prior/labels.csv):/tmp/labels.csv:ro" \
            -v "$(readlink -f ./$prior/cluster_tree.json):/tmp/cluster_tree.json:ro" \
            -v "$(readlink -f ./feature-config-t.json):/tmp/config.json:ro" \
            -v "$(readlink -f /tmp):/tmp/results" \
            node node dist/exportTree.js \
            --labelPath /tmp/labels.csv \
            --treePath /tmp/cluster_tree.json \
            --outPath /tmp/results/feature-tree-tmci-multi.svg \
            --configPath /tmp/config.json

        max_mem_export=$(get_max_mem_dual ${tmci_db_container_name} ${tmci_container_name})

        if [[ $max_mem_load -gt $max_mem_export ]]; then
            max_mem=$max_mem_load
        else
            max_mem=$max_mem_export
        fi

        docker rm ${tmci_container_name}

    done

    docker stop ${tmci_db_container_name}

    start_finish=$(docker container inspect ${tmci_db_container_name} | jq .[0].State.StartedAt,.[0].State.FinishedAt | tr '\n' ',')

    echo -e "tmci,with-5-features,${start_finish}${max_mem}," >> $result_file

    docker rm ${tmci_db_container_name}

    docker compose -f docker-compose.yaml down -v


done


# will run out of memory
# for i in {1..5}; do

#     max_mem="0"
#     start_time=$(date +"%Y-%m-%dT%T.%9NZ")

#     for feature in "${five_features[@]}"; do

#        docker run -d -v `pwd`:`pwd` -w `pwd` \
#         --user 1000 --name ${tmc_container_name} \
#         gregoryschwartz/too-many-cells:3.0.1.0 make-tree \
#         --prior ./$prior \
#         --matrix-path ./matrix-10x \
#         --feature-column 2 \
#         --draw-leaf "DrawItem (DrawContinuous [\"Cd4\"])" \
#         --labels-file ./$prior/labels.csv \
#         --dendrogram-output feature-tree-tmc.svg \
#         --output /tmp


#         process_max_mem=$(get_max_mem ${tmc_container_name})

#         if [[ $process_max_mem -gt $max_mem ]]; then
#             max_mem=$process_max_mem
#         fi

#         docker rm ${tmc_container_name}
#     done

#     end_time=$(date +"%Y-%m-%dT%T.%9NZ")

#     echo -e "tmc,with-5-features,${start_time},${end_time},${max_mem}," >> $result_file

# done
