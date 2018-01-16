#!/bin/bash

rsync -avzhP --exclude='data/burnins/*' bionc04.eva.mpg.de:/mnt/scratch/mp/nea-over-time/data .
