#!/bin/bash

input_files=("Signal_Selection1.txt" "Signal_Selection2.txt" "Signal_Selection3.txt" "Signal_Selection4.txt"
"Signal_Selection5.txt" "Signal_Selection6.txt" "Signal_Selection7.txt" "Signal_Selection8.txt" 
"Signal_Selection9.txt" "Signal_Selection10.txt" "Signal_Selection11.txt" "Signal_Selection12.txt"
 "Signal_Selection13.txt" "Signal_Selection14.txt" "Signal_Selection15.txt" "Signal_Selection16.txt" 
"Signal_Selection17.txt" "Signal_Selection18.txt" "Signal_Selection19.txt"
"Signal_Selection20.txt" "Signal_Selection21.txt" "Signal_Selection22.txt" "Signal_Selection23.txt" "Signal_Selection24.txt"
 "Signal_Selection25.txt" "Signal_Selection26.txt" "Signal_Selection27.txt" "Signal_Selection28.txt"
  "Signal_Selection29.txt" "Signal_Selection30.txt")

  
for file in "${input_files[@]}"
do
  echo "Running ROOT with $file"
  root -b -q "ScalarLept_wNC.C(\"$file\")"
  echo ""
done
