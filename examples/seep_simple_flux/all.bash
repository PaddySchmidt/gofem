#!/bin/bash

mpirun -np 4 gofem d2-simple-flux && GenVtu d2-simple-flux 1
mpirun -np 4 gofem d3-simple-flux && GenVtu d3-simple-flux 1
go run doplot.go
