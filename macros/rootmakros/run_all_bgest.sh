#!/bin/bash

for channel in mu ele
do
    for region in SR VR
    do
        root -l -b -q backgroundEstimation.C\(\"$channel\",\"$region\"\)
        root -l -b -q backgroundEstimation.C\(\"$channel\",\"$region\",\"btagging_totalUp\"\)
        root -l -b -q backgroundEstimation.C\(\"$channel\",\"$region\",\"btagging_totalDown\"\)
    done
done