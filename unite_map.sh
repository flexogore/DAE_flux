#!/bin/sh

# Clear the output file before starting
echo -n "" > all_mapping.txt

# Loop through each reaction folder
for rxn in $(ls -1 reaction_intermediates); do
    rxn_folder=reaction_intermediates/$rxn

    # Check if mapping.txt exists in the folder
    if [ ! -f "$rxn_folder/mapping.txt" ]; then
        echo "Skipping $rxn_folder: mapping.txt not found"
        continue
    fi

    # Process the mapping.txt file
    sed 's/,/\n/g' "$rxn_folder/mapping.txt" | sort > "$rxn_folder/mapping.sorted.txt"
    cat "$rxn_folder/mapping.sorted.txt" | sed "s/^/$rxn /" >> all_mapping.txt
done

# Generate sorted and analyzed files
sort all_mapping.txt > all_mapping.sorted.txt

grep ':N#' all_mapping.sorted.txt > all_mapping.N.sorted.txt
grep ':C#' all_mapping.sorted.txt > all_mapping.C.sorted.txt

sed 's/ .*//' all_mapping.N.sorted.txt | sort | uniq -c | sort -n > all_rxn_N_count.txt
sed 's/ .*//' all_mapping.C.sorted.txt | sort | uniq -c | sort -n > all_rxn_C_count.txt

sed 's/ *//; s/ .*//' all_rxn_N_count.txt | sort -n | uniq -c > all_rxn_N_count.histo
sed 's/ *//; s/ .*//' all_rxn_C_count.txt | sort -n | uniq -c > all_rxn_C_count.histo

sed 's/.* //; s/=/\n/' all_mapping.N.sorted.txt | sort -u > all_atoms.N.sorted.txt
#sed 's/.* //; s/=/\n/' all_mapping.C.sorted.txt | sort -u > all_atoms.C.sorted.txt
sed 's/.* //; s/=/\n/' all_mapping.C.sorted.txt | sort -u > all_atoms.C.sorted.txt
#sed 's/.*=//; s/[^[:print:]\n]*//g' all_mapping.C.sorted.txt | sort > all_atoms.C.sorted.txt

sed 's/.* //; s/=/\n/' all_mapping.N.sorted.txt | sort | uniq -c | sort -n > all_atoms.N.count.txt
sed 's/.* //; s/=/\n/' all_mapping.C.sorted.txt | sort | uniq -c | sort -n > all_atoms.C.count.txt


sed 's/ *//; s/ .*//' all_atoms.N.count.txt | sort -n | uniq -c > all_atoms.N.count.histo
sed 's/ *//; s/ .*//' all_atoms.C.count.txt | sort -n | uniq -c > all_atoms.C.count.histo
