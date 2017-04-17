#!/bin/bash
mongod --dbpath ../output/db --shutdown
rm .mongodb_host
