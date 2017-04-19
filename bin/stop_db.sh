#!/bin/bash
mongod --dbpath ../output/db --shutdown
rm ../output/.mongodb_host
