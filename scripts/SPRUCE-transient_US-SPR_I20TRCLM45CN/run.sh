#!/bin/bash
./SPRUCE-transient_US-SPR_I20TRCLM45CN.build
cp env_run.xmlnorestart env_run.xml
./SPRUCE-transient_US-SPR_I20TRCLM45CN.run
cp /home/u9x/output/SPRUCE-transient_US-SPR_I20TRCLM45CN/run/* /home/u9x/output/SPRUCE-transient_US-SPR_I20TRCLM45CN/run-norestart/
cp env_run.xmlrestart env_run.xml
./SPRUCE-transient_US-SPR_I20TRCLM45CN.submit

