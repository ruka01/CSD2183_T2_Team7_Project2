#!/bin/bash
# Usage: ./run.sh <input.csv> <target_vertices>
# Example: ./run.sh input.csv 17

set -e

INPUT_CSV="$1"
TARGET="$2"

if [ -z "$INPUT_CSV" ] || [ -z "$TARGET" ]; then
    echo "Usage: $0 <input.csv> <target_vertices>"
    exit 1
fi

if [ ! -f "$INPUT_CSV" ]; then
    echo "Error: input file '$INPUT_CSV' not found"
    exit 1
fi

# Folder name = input filename without extension
BASENAME=$(basename "$INPUT_CSV" .csv)
FOLDER="./$BASENAME"
mkdir -p "$FOLDER"

# Copy input CSV and View.html into the folder
cp "$INPUT_CSV" "$FOLDER/input.csv"

# Find View.html next to this script
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cp "$SCRIPT_DIR/View.html" "$FOLDER/View.html"

# Run simplify from inside the folder (it writes after.csv to cwd)
cd "$FOLDER"
"$SCRIPT_DIR/simplify" "input.csv" "$TARGET"

# Rename after.csv -> output.csv
mv after.csv output.csv

echo ""
echo "Done! Folder: $FOLDER"
echo "Open $FOLDER/View.html in your browser."
