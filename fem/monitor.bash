#!/bin/bash

FILES="e_u.go e_p.go e_up.go t_nurbs_test.go"

while true; do
    inotifywait -q -e modify $FILES
    echo
    echo
    go test -test.run="nurbs03"
done
