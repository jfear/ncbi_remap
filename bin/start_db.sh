#!/bin/bash
# Start mongo
mongod --dbpath ../output/db --logpath ../output/logs/mongoDB.log --port 27022 --fork
hostname > .mongodb_host
