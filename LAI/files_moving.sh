#!/bin/bash

# Check if an argument is passed
if [ "$#" -ne 1 ]; then
  echo "Usage: $0 <result_folder>"
  exit 1
fi

result_folder="$1"

# Directory che contiene i modelli, una per ogni chr
SOURCE_DIR="${result_folder}/gnomix_models"

# Directory dove salvare tutti i risultati MSP
DEST_DIR="${result_folder}/msp_folder"

mkdir -p "$DEST_DIR"

for model_dir in "$SOURCE_DIR"/chr*_model; do
  chr=$(basename "$model_dir" | sed 's/_model//')  # estrae "chrN"
  msp_file="$model_dir/query_results.msp"
  new_msp_file="$DEST_DIR/${chr}.msp"
  
  if [[ -f "$msp_file" ]]; then
    cp "$msp_file" "$new_msp_file"
    echo "✅ Copiato e rinominato $msp_file → $new_msp_file"
  else
    echo "❌ $chr: file MSP mancante"
  fi
done

echo "$DEST_DIR"
