#!/bin/bash

# cp -v <...> ~/.ipython/profile_${TEST_PROFILE}/...

conda install -y -c ${CONDA_CHANNEL_NAME} 03-id-hxn-collection

sudo mkdir -v -p /home/xf03id/
sudo chown -Rv $USER: /home/xf03id/
touch /home/xf03id/benchmark.out
