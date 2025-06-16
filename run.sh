#!/bin/bash

#export DPL_DEFAULT_PIPELINE_LENGTH=8

O="-b --configuration json://configuration.json --shm-segment-size 2684354560 --aod-memory-rate-limit 214748360"

o2-analysis-track-propagation $O | o2-analysis-timestamp $O | o2-analysis-cf-filter-correlations $O | o2-analysis-trackselection $O | o2-analysis-event-selection $O | o2-analysis-cf-correlations $O | o2-analysis-tracks-extra-v002-converter $O --aod-file @input.txt



