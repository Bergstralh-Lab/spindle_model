#!/usr/bin/env bash
# Rename files in a folder by replacing a prefix.
# Example: Mask_1.jpg -> Spindle_1.jpg
#
# Usage:
#   ./rename_files.sh <folder> <old_prefix> <new_prefix> [--dry-run]

set -euo pipefail

if [[ $# -lt 3 ]]; then
    echo "Usage: $0 <folder> <old_prefix> <new_prefix> [--dry-run]"
    echo "Example: $0 ./images Mask_ Spindle_"
    exit 1
fi

FOLDER="$1"
OLD_PREFIX="$2"
NEW_PREFIX="$3"
DRY_RUN=false

if [[ "${4:-}" == "--dry-run" ]]; then
    DRY_RUN=true
fi

if [[ ! -d "$FOLDER" ]]; then
    echo "Error: '$FOLDER' is not a directory"
    exit 1
fi

COUNT=0

for filepath in "$FOLDER"/"$OLD_PREFIX"*; do
    [[ -e "$filepath" ]] || { echo "No files found with prefix '$OLD_PREFIX' in '$FOLDER'"; exit 0; }
    filename=$(basename "$filepath")
    new_name="${NEW_PREFIX}${filename#"$OLD_PREFIX"}"
    echo "  $filename  ->  $new_name"
    if [[ "$DRY_RUN" == false ]]; then
        mv "$filepath" "$FOLDER/$new_name"
    fi
    COUNT=$((COUNT + 1))
done

echo ""
if [[ "$DRY_RUN" == true ]]; then
    echo "[Dry run] $COUNT file(s) would be renamed. Remove --dry-run to apply."
else
    echo "$COUNT file(s) renamed."
fi