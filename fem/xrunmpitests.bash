#!/bin/bash

go build -o /tmp/gofem/t_bh16_main t_bh16_main.go && mpirun -np 3 /tmp/gofem/t_bh16_main
go build -o /tmp/gofem/t_spo751_main t_spo751_main.go && mpirun -np 3 /tmp/gofem/t_spo751_main