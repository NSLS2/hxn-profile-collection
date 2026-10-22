#!/bin/bash

export AZURE_TESTING=1

sudo mkdir -v -p /nsls2/xf05id1/shared/config/
sudo mkdir -v -p /home/xf03id/
sudo chown -Rv $USER: /home/xf03id/
touch /home/xf03id/benchmark.out
