#!/bin/bash

FILES="shp.go nurbs.go t_nurbs_test.go"

while true; do
    inotifywait -q -e modify $FILES
    echo
    echo
    go test -test.run="nurbs04"
done
