All benchmarks were performed on were performed on an AWS EC2 r6idn.16xlarge instance running Ubuntu 20.04 and
Docker 20.10.17 with 64x Intel Xeon Platinum 8375C CPU @ 2.90 GHz and 534 GB RAM.

Running these scripts requires the latest too-many-cells-interactive Docker image, which can be built by cloning the (public GitHub repository)[https://github.com/schwartzlab-methods/too-many-cells-interactive]. Replace the tag in the [docker-compose.yaml](./docker-compose.yaml) with your local image tag. You will also need to pull in the latest toomanycells image, for instance with the command `docker pull gregoryschwartz/too-many-cells:latest`. 

When comparing performance between toomanycells and toomanycellsinteractive, we benchmarks TMCI as a node application without a frontend. Both applications are run in Docker containers and metrics are gathered from [system control groups](https://docs.docker.com/config/containers/runmetrics/#control-groups). 

When comparing with cellxgene and cirrcoumulus, benchmarks were made with the TMCI webapp. We used the Docker control groups to measure memory usage and journald logs to measure processing times, beginning with the start of execution and concluding with the launch of the webserver. 
