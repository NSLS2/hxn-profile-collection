#!/bin/bash

export AZURE_TESTING=1

# TODO: watch for https://github.com/NSLS-II-HXN/hxnfly/pull/8 to be
# merged/released, then this install can be removed from here.
python3 -m pip install --no-deps -vv git+https://github.com/NSLS-II-HXN/hxnfly@master
python3 -m pip install mpl_qtthread

sudo mkdir -v -p /home/xf03id/
sudo chown -Rv $USER: /home/xf03id/
touch /home/xf03id/benchmark.out
