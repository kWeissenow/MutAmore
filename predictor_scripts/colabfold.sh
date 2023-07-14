#!/bin/bash

# This example uses ColabFold's single-sequence mode.
# In order to generate multiple-sequence alignments, please consider installing MMSeqs2 and using colabfold_search as described in the ColabFold documentation


# run ColabFold
colabfold_batch --msa-mode single_sequence MUTAMORE_INPUT MUTAMORE_OUTPUT

# rename output files
cd MUTAMORE_OUTPUT
for file in *.pdb; do
  if [[ $file == *rank_001* ]]; then
      new_name="${file%%_unrelaxed*}.pdb"
      cp "$file" "$new_name"
  fi
done
