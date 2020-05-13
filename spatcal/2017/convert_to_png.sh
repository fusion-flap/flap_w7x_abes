#!/bin/bash
for file in *.bmp; do
convert "$file" "$(basename "${file/.bmp}")".png
sleep 1
done
