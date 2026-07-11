#!/bin/bash
# Genera TUTTI i mountain plot nella cartella unica ld_mountain_plots_all/:
#   - 6 varianti LD per hit: {allsamples, hgdp} x {OR, BYp, P}      (veloce)
#   - 1 plot entropia per hit: LD globale + entropia aplotipica     (pesante)
# Da lanciare in uno screen. Esegue tutto in sequenza.
#
#   screen -S ldplots
#   bash run_all_ld_plots.sh

cd /private/home/rsmerigl/codes/cleaned_codes/Admixture_mapping/post_analysis

echo "=== [1/2] 6 plot LD (allsamples + HGDP x OR/BYp/P) ==="
python ld_mountain_plots_all.py

echo ""
echo "=== [2/2] plot entropia (estrazione aplotipi: pesante, ~5-10 min/hit) ==="
python ld_entropy_plots.py

echo ""
echo "=== consolidamento: rimuovo le vecchie cartelle (ora superate) ==="
for d in ld_mountain_plots_6wind ld_mountain_plots_6wind_signif ld_mountain_plots_6wind_rawp; do
  if [ -d "$d" ]; then echo "  rimuovo $d/"; rm -rf "$d"; fi
done

echo ""
echo "DONE. Tutti i plot sono in: ld_mountain_plots_all/"
ls ld_mountain_plots_all/*.png 2>/dev/null | wc -l | xargs echo "PNG totali:"
