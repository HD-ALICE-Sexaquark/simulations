#!/bin/bash

while true; do
    squeue -h -u "$USER" -t PD \
        | grep "(launch failed requeued held)" \
        | awk '{print $1}' \
        | xargs -r scontrol release
    sleep 300
done
