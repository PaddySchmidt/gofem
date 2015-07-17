#!/bin/bash

set -e

SIMFILES="sg111 sg114 sg1121"

for f in $SIMFILES; do
    echo
    echo "[1;33m>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> $f <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<[0m"
    gofem $f
    go run doplot-"$f".go
done
