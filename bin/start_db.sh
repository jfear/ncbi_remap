#!/bin/bash

# Change ulimits
ulimit -u 4000

# Start mongo
numactl --interleave=all mongod --dbpath ../output/db --logpath ../output/logs/mongoDB.log --port 27022 --fork
hostname > ../output/.mongodb_host
