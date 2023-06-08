#! /usr/bin/env bash

set -euo pipefail

# Benchmark TMC and TMI against tm-subset dataset

result_file=_results-tmc-subset.csv
tmc_container_name=tmc
tmci_container_name=tmci_node
tmci_db_container_name=postgres

echo "type,start,finish,RAM(b)," > $result_file

ten_features=("Cd4" "Apoe" "Ahr" "Myc" "Serp1" "C6" "Neb" "Brca" "Brca2" "Cad")


cgroup_mem_root=$(grep cgroup /proc/mounts | grep memory | awk '{print $2}')/docker

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

for i in {1..10}; do
    # we run the process in the background and watch for exit on mem usage loop and remove manually after benchmarking time
    docker run -d -v `pwd`:`pwd` -w `pwd` \
        --user 1000 --name ${tmc_container_name} \
        gregoryschwartz/too-many-cells:3.0.1.0 make-tree \
        --prior ./data/prior/tm-subset \
        --labels-file ./data/prior/tm-subset/labels.csv \
        --dendrogram-output newtree.svg \
        --output ./data

    max_mem=$(get_max_mem ${tmc_container_name})

    start_finish=$(docker container inspect ${tmc_container_name} | jq .[0].State.StartedAt,.[0].State.FinishedAt | tr '\n' ',')

    docker rm ${tmc_container_name}

    echo -e "tmc,basic,${start_finish}${max_mem}," >> $result_file

done

for i in {1..10}; do
    docker run -d -v `pwd`:`pwd` -w `pwd` \
        --user 1000 --name ${tmc_container_name} \
        gregoryschwartz/too-many-cells:3.0.1.0 make-tree \
        --prior ./data/prior/tm-subset/ \
        --min-size 1000 \
        --labels-file ./data/prior/tm-subset/labels.csv \
        --dendrogram-output newtree.svg \
        --output ./data

    max_mem=$(get_max_mem ${tmc_container_name})

    start_finish=$(docker container inspect ${tmc_container_name} | jq .[0].State.StartedAt,.[0].State.FinishedAt | tr '\n' ',')

    docker rm ${tmc_container_name}

    echo -e "tmc,prune,${start_finish}${max_mem}," >> $result_file

done

for i in {1..10}; do
    docker run -d -v `pwd`:`pwd` -w `pwd` \
        --user 1000 --name ${tmc_container_name} \
        gregoryschwartz/too-many-cells:3.0.1.0 make-tree \
        --matrix-path ./data/tm-subset/Bladder-10X_P4_3 \
        --matrix-path ./data/tm-subset/Heart_and_Aorta-10X_P7_4 \
        --matrix-path ./data/tm-subset/Kidney-10X_P7_5 \
        --matrix-path ./data/tm-subset/Limb_Muscle-10X_P7_15 \
        --matrix-path ./data/tm-subset/Liver-10X_P7_1 \
        --matrix-path ./data/tm-subset/Mammary_Gland-10X_P7_13 \
        --matrix-path ./data/tm-subset/Marrow-10X_P7_3 \
        --matrix-path ./data/tm-subset/Spleen-10X_P7_6 \
        --matrix-path ./data/tm-subset/Thymus-10X_P7_11 \
        --matrix-path ./data/tm-subset/Tongue-10X_P7_10 \
        --labels-file ./data/prior/tm-subset/labels.csv \
        --prior ./data/prior/tm-subset \
        --dendrogram-output newtree.svg \
        --output ./data

    max_mem=$(get_max_mem ${tmc_container_name})

    start_finish=$(docker container inspect ${tmc_container_name} | jq .[0].State.StartedAt,.[0].State.FinishedAt | tr '\n' ',')

    docker rm ${tmc_container_name}

    echo -e "tmc,with-matrix,${start_finish}${max_mem}," >> $result_file

done


for i in {1..10}; do

    docker-compose -f docker-compose.yaml \
        run -d --entrypoint="" \
        --name ${tmci_container_name} \
        --no-deps \
        -v "$(readlink -f ./data/prior/tm-subset/labels.csv):/tmp/labels.csv:ro" \
        -v "$(readlink -f ./data/prior/tm-subset/cluster_tree.json):/tmp/cluster_tree.json:ro" \
        -v "$(readlink -f ./data/blank-config.json):/tmp/config.json:ro" \
        -v "$(readlink -f ./data):/tmp/results" \
        node node dist/exportTree.js \
        --labelPath /tmp/labels.csv \
        --treePath /tmp/cluster_tree.json \
        --outPath /tmp/results/export667.svg \
        --configPath /tmp/config.json

    max_mem=$(get_max_mem ${tmci_container_name})

    docker stop ${tmci_container_name}
    start_finish=$(docker container inspect ${tmci_container_name} | jq .[0].State.StartedAt,.[0].State.FinishedAt | tr '\n' ',')
    docker rm $tmci_container_name

    echo -e "tmci,basic,${start_finish}${max_mem}," >> $result_file

done

for i in {1..10}; do

    docker-compose -f docker-compose.yaml \
        run -d --entrypoint="" \
        --name ${tmci_container_name} \
        --no-deps \
        -v "$(readlink -f ./data/prior/tm-subset/labels.csv):/tmp/labels.csv:ro" \
        -v "$(readlink -f ./data/prior/tm-subset/cluster_tree.json):/tmp/cluster_tree.json:ro" \
        -v "$(readlink -f ./data/blank-config.json):/tmp/config.json:ro" \
        -v "$(readlink -f ./data):/tmp/results" \
        node node dist/exportTree.js \
        --labelPath /tmp/labels.csv \
        --treePath /tmp/cluster_tree.json \
        --outPath /tmp/results/export667.svg \
        --configPath /tmp/config.json

    max_mem=$(get_max_mem ${tmci_container_name})

    docker stop ${tmci_container_name}
    start_finish=$(docker container inspect ${tmci_container_name} | jq .[0].State.StartedAt,.[0].State.FinishedAt | tr '\n' ',')

    docker rm $tmci_container_name

    echo -e "tmci,prune,${start_finish}${max_mem}," >> $result_file

done

for i in {1..10}; do

    # start postgres in the background
    docker-compose -f docker-compose.yaml run --name ${tmci_db_container_name} -d postgres

    #load data
    docker-compose -f docker-compose.yaml run -d --name ${tmci_container_name} \
        -v "$(readlink -f ./data/tm-subset):/usr/data/matrices:ro" \
        node --init

    max_mem_load=$(get_max_mem_dual ${tmci_db_container_name} ${tmci_container_name})

    docker rm ${tmci_container_name}

    # now run export
    docker-compose -f docker-compose.yaml \
        run -d --entrypoint="" \
        --name ${tmci_container_name} \
        -v "$(readlink -f ./data/prior/tm-subset/labels.csv):/tmp/labels.csv:ro" \
        -v "$(readlink -f ./data/prior/tm-subset/cluster_tree.json):/tmp/cluster_tree.json:ro" \
        -v "$(readlink -f ./data/blank-config.json):/tmp/config.json:ro" \
        -v "$(readlink -f ./data):/tmp/results" \
        node node dist/exportTree.js \
        --labelPath /tmp/labels.csv \
        --treePath /tmp/cluster_tree.json \
        --outPath /tmp/results/export667.svg \
        --configPath /tmp/config.json

    max_mem_export=$(get_max_mem_dual ${tmci_db_container_name} ${tmci_container_name})

    docker stop ${tmci_db_container_name}

    start_finish=$(docker container inspect ${tmci_db_container_name} | jq .[0].State.StartedAt,.[0].State.FinishedAt | tr '\n' ',')

    docker-compose -f docker-compose.yaml down -v

    if [[ $max_mem_load -gt $max_mem_export ]]; then
        max_mem=$max_mem_load
    else
        max_mem=$max_mem_export
    fi

    echo -e "tmci,with-matrix,${start_finish}${max_mem}," >> $result_file

done

for i in {1..10}; do

    # start postgres in the background
    docker-compose -f docker-compose.yaml run --name ${tmci_db_container_name} -d postgres

    #load data
    docker-compose -f docker-compose.yaml run -d --name ${tmci_container_name} \
        -v "$(readlink -f ./data/tm-subset/):/usr/data/matrices:ro" \
        node --init

    max_mem_load=$(get_max_mem_dual ${tmci_db_container_name} ${tmci_container_name})

    docker rm ${tmci_container_name}

    # now run export
    docker-compose -f docker-compose.yaml \
        run -d --entrypoint="" \
        --name ${tmci_container_name} \
        -v "$(readlink -f ./data/prior/tm-subset/labels.csv):/tmp/labels.csv:ro" \
        -v "$(readlink -f ./data/prior/tm-subset/cluster_tree.json):/tmp/cluster_tree.json:ro" \
        -v "$(readlink -f ./data/feature-config.json):/tmp/config.json:ro" \
        -v "$(readlink -f ./data):/tmp/results" \
        node node dist/exportTree.js \
        --labelPath /tmp/labels.csv \
        --treePath /tmp/cluster_tree.json \
        --outPath /tmp/results/feature-tree-tmci.svg \
        --configPath /tmp/config.json

    max_mem_export=$(get_max_mem_dual ${tmci_db_container_name} ${tmci_container_name})

    docker stop ${tmci_db_container_name}

    start_finish=$(docker container inspect ${tmci_db_container_name} | jq .[0].State.StartedAt,.[0].State.FinishedAt | tr '\n' ',')

    docker-compose -f docker-compose.yaml down -v

    if [[ $max_mem_load -gt $max_mem_export ]]; then
        max_mem=$max_mem_load
    else
        max_mem=$max_mem_export
    fi

    echo -e "tmci,with-feature,${start_finish}${max_mem}," >> $result_file

done

for i in {1..10}; do

    docker run -d -v `pwd`:`pwd` -w `pwd` \
        --user 1000 --name ${tmc_container_name} \
        gregoryschwartz/too-many-cells:3.0.1.0 make-tree \
        --prior ./data/prior/tm-subset \
        --matrix-path ./data/tm-subset/Bladder-10X_P4_3 \
        --matrix-path ./data/tm-subset/Heart_and_Aorta-10X_P7_4 \
        --matrix-path ./data/tm-subset/Kidney-10X_P7_5 \
        --matrix-path ./data/tm-subset/Limb_Muscle-10X_P7_15 \
        --matrix-path ./data/tm-subset/Liver-10X_P7_1 \
        --matrix-path ./data/tm-subset/Mammary_Gland-10X_P7_13 \
        --matrix-path ./data/tm-subset/Marrow-10X_P7_3 \
        --matrix-path ./data/tm-subset/Spleen-10X_P7_6 \
        --matrix-path ./data/tm-subset/Thymus-10X_P7_11 \
        --matrix-path ./data/tm-subset/Tongue-10X_P7_10 \
        --labels-file ./data/prior/tm-subset/labels.csv \
        --feature-column 2 \
        --draw-leaf "DrawItem (DrawContinuous [\"Cd4\"])" \
        --dendrogram-output feature-tree-tmc.svg \
        --output ./data/

    max_mem=$(get_max_mem ${tmc_container_name})

    start_finish=$(docker container inspect ${tmc_container_name} | jq .[0].State.StartedAt,.[0].State.FinishedAt | tr '\n' ',')

    docker rm ${tmc_container_name}

    echo -e "tmc,with-feature,${start_finish}${max_mem}," >> $result_file

done



for i in {1..10}; do

    max_mem="0"

     # start postgres in the background
    docker-compose -f docker-compose.yaml run --name ${tmci_db_container_name} -d postgres

    #load data
    docker-compose -f docker-compose.yaml run -d --name ${tmci_container_name} \
        -v "$(readlink -f ./data/tm-subset/):/usr/data/matrices:ro" \
        node --init

    max_mem_load=$(get_max_mem_dual ${tmci_db_container_name} ${tmci_container_name})

    docker rm ${tmci_container_name}

    sleep .1

    for feature in "${ten_features[@]}"; do

        jq --arg feature "$feature" 'setpath(["features"];[$feature])' ./data/feature-config.json > ./data/feature-config-t.json

        docker-compose -f docker-compose.yaml \
            run -d --entrypoint="" \
            --name ${tmci_container_name} \
            -v "$(readlink -f ./data/prior/tm-subset/labels.csv):/tmp/labels.csv:ro" \
            -v "$(readlink -f ./data/prior/tm-subset/cluster_tree.json):/tmp/cluster_tree.json:ro" \
            -v "$(readlink -f ./data/feature-config-t.json):/tmp/config.json:ro" \
            -v "$(readlink -f ./data):/tmp/results" \
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

    echo -e "tmci,with-10-features,${start_finish}${max_mem}," >> $result_file

    docker-compose -f docker-compose.yaml down -v


done

for i in {1..5}; do

    max_mem="0"
    start_time=$(date +"%Y-%m-%dT%T.%9NZ")

    for feature in "${ten_features[@]}"; do

       docker run -d -v `pwd`:`pwd` -w `pwd` \
        --user 1000 --name ${tmc_container_name} \
        gregoryschwartz/too-many-cells:3.0.1.0 make-tree \
        --prior ./data/prior/tm-subset \
--matrix-path ./data/tm-subset/Bladder-10X_P4_3 \
        --matrix-path ./data/tm-subset/Heart_and_Aorta-10X_P7_4 \
        --matrix-path ./data/tm-subset/Kidney-10X_P7_5 \
        --matrix-path ./data/tm-subset/Limb_Muscle-10X_P7_15 \
        --matrix-path ./data/tm-subset/Liver-10X_P7_1 \
        --matrix-path ./data/tm-subset/Mammary_Gland-10X_P7_13 \
        --matrix-path ./data/tm-subset/Marrow-10X_P7_3 \
        --matrix-path ./data/tm-subset/Spleen-10X_P7_6 \
        --matrix-path ./data/tm-subset/Thymus-10X_P7_11 \
        --matrix-path ./data/tm-subset/Tongue-10X_P7_10 \
        --labels-file ./data/prior/tm-subset/labels.csv \
        --feature-column 2 \
        --draw-leaf "DrawItem (DrawContinuous [\"Cd4\"])" \
        --labels-file ./data/prior/tm-subset/labels.csv \
        --dendrogram-output feature-tree-tmc.svg \
        --output ./data/


        process_max_mem=$(get_max_mem ${tmc_container_name})

        if [[ $process_max_mem -gt $max_mem ]]; then
            max_mem=$process_max_mem
        fi

        docker rm ${tmc_container_name}
    done

    end_time=$(date +"%Y-%m-%dT%T.%9NZ")

    echo -e "tmc,with-10-features,${start_time},${end_time},${max_mem}," >> $result_file

done

