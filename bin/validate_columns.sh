#!/bin/bash
# validate_columns.sh

required_columns="$1"
gwas_input="$2"

header=$(zcat "$gwas_input" | head -n 1)
missing_columns=""

for col in $required_columns; do
    if ! echo "$header" | grep -wq "$col"; then
        missing_columns+="$col "
    fi
done

if [ -n "$missing_columns" ]; then
    echo "Error: The following required columns from $gwas_input are missing: $missing_columns"
    exit 1
fi

