#!/bin/bash

# Folder containing your input files
input_folder="/exp/dune/data/users/rrichi/MLProject/Selection_files"

# Collect all .txt files in that folder into an array
input_files=("$input_folder"/*.txt)

# Loop over each input file
for file in "${input_files[@]}"
do
  echo "Running ROOT with $file"
  root -b -q "VectorLept_wNC.C(\"$file\")"
  echo ""
done
