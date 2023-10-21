#!/bin/bash

input_file=template.cfg

output_file="$1".cfg

# Replace all of the @uniform tags with a random integer between -100 and 100:
generate_random_number() {
    echo $((RANDOM % 201 - 100))
}

cp "$input_file" "$output_file"

while grep -q '@uniform' "$output_file"; do
    sed -i "0,/@uniform/s//$(generate_random_number)/" "$output_file"
done

# Run a fit
mpirun --mca btl_openib_allow_ib 1 ../bin/fitMPI -c $output_file
