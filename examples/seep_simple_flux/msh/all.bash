#!/bin/bash

go run gen2d.go && cp /tmp/gemlab/d2-coarse.msh .
go run gen3d.go && cp /tmp/gemlab/d3-coarse.msh .
