#!/bin/bash

export AZURE_TESTING=1

sudo mkdir -p /etc/bluesky
if [ ! -e /etc/bluesky/kafka.yml ] && [ ! -L /etc/bluesky/kafka.yml ]; then
    sudo sh -c 'umask 022; set -C; printf "%s\n" \
        "bootstrap_servers:" \
        "  - localhost:9092" \
        "runengine_producer_config:" \
        "  sasl.password: ci-placeholder" > /etc/bluesky/kafka.yml' || exit 1
fi

sudo mkdir -v -p /nsls2/xf05id1/shared/config/
sudo mkdir -v -p /home/xf03id/
sudo chown -Rv $USER: /home/xf03id/
touch /home/xf03id/benchmark.out
