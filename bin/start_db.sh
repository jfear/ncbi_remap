#!/bin/bash

# Change ulimits
ulimit -u 6400 -n 7000

# Start mongo
numactl --interleave=all mongod --dbpath ../output/db --logpath ../output/logs/mongoDB.log --port 27022 --fork
hostname > ../output/.mongodb_host
