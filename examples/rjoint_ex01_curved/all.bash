#!/bin/bash

set -e

SIMFILES="rjoint01"

for f in $SIMFILES; do
    echo
    echo "[1;33m>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> $f <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<[0m"
    gofem $f
    go run doplot.go
    GenVtu $f
done
