#!/bin/bash

# Directory che contiene i modelli, una per ogni chr
SOURCE_DIR="/private/groups/ioannidislab/smeriglio/tests_files/results/gnomix_models"

# Directory dove salvare tutti i risultati MSP
DEST_DIR="/private/groups/ioannidislab/smeriglio/tests_files/results/msp_folder"

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
