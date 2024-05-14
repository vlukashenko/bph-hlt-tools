#!/bin/bash

modify_csv() {
    local file="$1"
    local name="$2"
    local number="$3"

		echo $file
		echo $name
		echo $number

    # Check if the file exists
    if [ ! -e "$file" ]; then
        echo "Error: File '$file' not found."
        exit 1
    fi

    # Look for the line with the specified Name and modify trailing numbers
    awk -v name="$name" -v number="$number" -F',' '
        BEGIN { OFS=FS }
        $2 == name {
            for (i = 4; i <= NF; i++) {
                $i = number
            }
        }
        { print }
    ' "$file" > "$file.tmp"

    # Rename the temporary file to the original file
    mv "$file.tmp" "$file"

		dos2unix "$file"
}

# Example usage:
# modify_csv your_file.csv L1_SingleMuCosmics 42
modify_csv "$1" "$2" "$3"
