import argparse
import glob
import os
from collections import OrderedDict, defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from snputils import MSPReader


class DatasetAncestrySummary:
    """
    Utility that streams MSP files for a dataset and extracts the reviewer-facing
    metrics: global ancestry proportions, tract length distributions, and tract counts.
    """

    def __init__(self, dataset_name: str):
        self.dataset_name = dataset_name
        # id -> label map (maintain insertion order for presentation)
        self.ancestry_map: OrderedDict[int, str] = OrderedDict()

        self.n_samples: Optional[int] = None
        self.n_haplotypes: Optional[int] = None
        self.sample_ids: Optional[Tuple[str, ...]] = None

        self.total_window_bp: float = 0.0  # Haploid genome length covered by windows

        # Arrays keyed by ancestry label (not id) for easier presentation.
        self.bp_per_sample: Dict[str, np.ndarray] = {}
        self.tract_counts: Dict[str, np.ndarray] = {}
        self.tract_lengths_bp: Dict[str, List[float]] = defaultdict(list)

        # Running state for tract construction
        self.current_ancestry: Optional[np.ndarray] = None
        self.current_length_bp: Optional[np.ndarray] = None
        self.last_chromosome: Optional[str] = None

    def consume_file(self, path: str) -> None:
        reader = MSPReader(path)
        lai_obj = reader.read()

        self._bootstrap_ancestry_map(lai_obj)
        self._initialize_arrays(lai_obj)

        lai_matrix = np.asarray(lai_obj.lai, dtype=np.int32)
        window_lengths = _extract_window_lengths(lai_obj)
        chromosomes = (
            np.asarray(lai_obj.chromosomes)
            if getattr(lai_obj, "chromosomes", None) is not None
            else None
        )

        if window_lengths.shape[0] != lai_matrix.shape[0]:
            raise ValueError(
                f"Window length size mismatch in {path}: "
                f"{window_lengths.shape[0]} rows vs {lai_matrix.shape[0]} LAI windows."
            )

        self.total_window_bp += float(window_lengths.sum())
        self._update_global_bp(lai_matrix, window_lengths)
        self._update_tracts(lai_matrix, window_lengths, chromosomes)

        # Reset run state so that files can be processed independently.
        self._flush_all_runs()

    def summarize_rows(self) -> List[Dict[str, object]]:
        if not self.ancestry_map:
            return []

        diploid_bp = self.total_window_bp * 2.0
        rows: List[Dict[str, object]] = []

        for ancestry_id, label in self.ancestry_map.items():
            bp_array = self.bp_per_sample.get(label)
            if bp_array is None and self.n_samples is not None:
                bp_array = np.zeros(self.n_samples, dtype=np.float64)

            if bp_array is not None and diploid_bp > 0:
                percentages = (bp_array / diploid_bp) * 100.0
            else:
                percentages = np.array([], dtype=np.float64)

            percent_mean, percent_iqr = _mean_and_iqr(percentages)
            percent_iqr_fmt = _format_interval(percent_iqr, suffix="%")

            lengths_bp = np.asarray(self.tract_lengths_bp.get(label, []), dtype=np.float64)
            lengths_mb = lengths_bp / 1e6 if lengths_bp.size else lengths_bp
            tract_mean, tract_iqr = _mean_and_iqr(lengths_mb)
            tract_iqr_fmt = _format_interval(tract_iqr, suffix=" Mb")

            counts = self.tract_counts.get(label)
            if counts is None and self.n_samples is not None:
                counts = np.zeros(self.n_samples, dtype=np.int64)
            tract_count_mean = float(np.mean(counts)) if counts is not None else np.nan

            rows.append(
                {
                    "Dataset": self.dataset_name,
                    "Ancestry": label,
                    "Mean_percent_genome": None
                    if np.isnan(percent_mean)
                    else round(percent_mean, 2),
                    "Percent_IQR": percent_iqr_fmt,
                    "Mean_tract_length_Mb": None
                    if np.isnan(tract_mean)
                    else round(tract_mean, 3),
                    "Tract_length_IQR_Mb": tract_iqr_fmt,
                    "Mean_tracts_per_individual": None
                    if np.isnan(tract_count_mean)
                    else round(tract_count_mean, 1),
                }
            )
        return rows

    def _bootstrap_ancestry_map(self, lai_obj) -> None:
        if self.ancestry_map:
            return

        file_map = getattr(lai_obj, "ancestry_map", None)
        if isinstance(file_map, dict):
            for key, label in file_map.items():
                try:
                    idx = int(key)
                except (TypeError, ValueError):
                    continue
                self.ancestry_map[idx] = str(label)
        elif isinstance(file_map, (list, tuple)):
            for idx, label in enumerate(file_map):
                self.ancestry_map[idx] = str(label)

    def _initialize_arrays(self, lai_obj) -> None:
        lai_matrix = np.asarray(lai_obj.lai)
        n_haplotypes = lai_matrix.shape[1]
        if n_haplotypes % 2 != 0:
            raise ValueError("The LAI matrix must contain an even number of haplotypes (diploid data).")

        n_samples = n_haplotypes // 2

        if self.n_samples is None:
            self.n_samples = n_samples
            self.n_haplotypes = n_haplotypes
            self.sample_ids = tuple(lai_obj.samples)
            self.current_ancestry = np.full(n_haplotypes, -1, dtype=np.int32)
            self.current_length_bp = np.zeros(n_haplotypes, dtype=np.float64)

            for label in self.ancestry_map.values():
                self._ensure_label_allocations(label)
        else:
            if self.n_samples != n_samples:
                raise ValueError(
                    f"Inconsistent number of samples for dataset {self.dataset_name}: "
                    f"expected {self.n_samples}, found {n_samples}."
                )
            if self.n_haplotypes != n_haplotypes:
                raise ValueError(
                    f"Inconsistent number of haplotypes for dataset {self.dataset_name}."
                )
            if self.sample_ids != tuple(lai_obj.samples):
                raise ValueError("Samples differ across MSP files for the same dataset.")

    def _ensure_label_allocations(self, label: str) -> None:
        if self.n_samples is None:
            return

        if label not in self.bp_per_sample:
            self.bp_per_sample[label] = np.zeros(self.n_samples, dtype=np.float64)
        if label not in self.tract_counts:
            self.tract_counts[label] = np.zeros(self.n_samples, dtype=np.int64)
        if label not in self.tract_lengths_bp:
            self.tract_lengths_bp[label] = []

    def _label_for_ancestry(self, ancestry_id: int) -> Optional[str]:
        if ancestry_id < 0:
            return None
        if ancestry_id not in self.ancestry_map:
            self.ancestry_map[ancestry_id] = str(ancestry_id)
            self._ensure_label_allocations(self.ancestry_map[ancestry_id])
        return self.ancestry_map[ancestry_id]

    def _update_global_bp(self, lai_matrix: np.ndarray, window_lengths: np.ndarray) -> None:
        if self.n_samples is None or self.current_ancestry is None:
            raise RuntimeError("Aggregator not initialized before updating statistics.")

        hap_length = window_lengths[:, None].astype(np.float64)

        for ancestry_id, label in self.ancestry_map.items():
            self._ensure_label_allocations(label)
            match = (lai_matrix == ancestry_id).astype(np.float64)
            bp_per_hap = (match * hap_length).sum(axis=0)
            bp_per_sample = bp_per_hap.reshape(self.n_samples, 2).sum(axis=1)
            self.bp_per_sample[label] += bp_per_sample

    def _update_tracts(
        self, lai_matrix: np.ndarray, window_lengths: np.ndarray, chromosomes: Optional[np.ndarray]
    ) -> None:
        if self.current_ancestry is None or self.current_length_bp is None:
            raise RuntimeError("Aggregator not initialized before updating statistics.")

        n_windows = lai_matrix.shape[0]

        for idx in range(n_windows):
            chrom = chromosomes[idx] if chromosomes is not None else None
            chrom = str(chrom) if chrom is not None else chrom

            if chrom is not None:
                if self.last_chromosome is None:
                    self.last_chromosome = chrom
                elif chrom != self.last_chromosome:
                    self._flush_all_runs()
                    self.last_chromosome = chrom

            window_bp = float(window_lengths[idx])
            row = lai_matrix[idx]

            same_mask = row == self.current_ancestry
            if np.any(same_mask):
                self.current_length_bp[same_mask] += window_bp

            changed_idx = np.flatnonzero(~same_mask)
            if changed_idx.size:
                self._flush_runs(changed_idx)
                self.current_ancestry[changed_idx] = row[changed_idx]
                self.current_length_bp[changed_idx] = window_bp

    def _flush_runs(self, hap_indices: np.ndarray) -> None:
        if hap_indices.size == 0 or self.current_length_bp is None or self.current_ancestry is None:
            return

        lengths = self.current_length_bp[hap_indices]
        ancestries = self.current_ancestry[hap_indices]

        valid_mask = (lengths > 0) & (ancestries >= 0)
        if not np.any(valid_mask):
            self.current_length_bp[hap_indices] = 0.0
            return

        hap_valid = hap_indices[valid_mask]
        lengths = lengths[valid_mask]
        ancestries = ancestries[valid_mask].astype(int)

        sample_ids = (hap_valid // 2).astype(int)
        order = np.argsort(ancestries)
        ancestries = ancestries[order]
        lengths = lengths[order]
        sample_ids = sample_ids[order]

        unique_ids, start_idx = np.unique(ancestries, return_index=True)
        for i, ancestry_id in enumerate(unique_ids):
            label = self._label_for_ancestry(int(ancestry_id))
            if label is None:
                continue
            idx_start = start_idx[i]
            idx_end = start_idx[i + 1] if i + 1 < len(start_idx) else lengths.size
            tract_lengths = lengths[idx_start:idx_end]
            related_samples = sample_ids[idx_start:idx_end]

            self._ensure_label_allocations(label)
            self.tract_lengths_bp[label].extend(float(val) for val in tract_lengths)

            counts = np.bincount(related_samples, minlength=self.n_samples)
            self.tract_counts[label] += counts

        self.current_length_bp[hap_indices] = 0.0

    def _flush_all_runs(self) -> None:
        if self.current_length_bp is None or self.current_ancestry is None:
            return
        active_idx = np.flatnonzero(self.current_length_bp > 0)
        if active_idx.size:
            self._flush_runs(active_idx)
        self.current_ancestry.fill(-1)
        self.current_length_bp.fill(0.0)
        self.last_chromosome = None


def _extract_window_lengths(lai_obj) -> np.ndarray:
    physical_pos = getattr(lai_obj, "physical_pos", None)
    if physical_pos is not None and len(physical_pos):
        pos = np.asarray(physical_pos, dtype=np.float64)
        if pos.ndim != 2 or pos.shape[1] != 2:
            raise ValueError("physical_pos must be a (n_windows, 2) array of (start, end).")
        lengths = pos[:, 1] - pos[:, 0]
    else:
        window_sizes = getattr(lai_obj, "window_sizes", None)
        if window_sizes is None or not len(window_sizes):
            raise ValueError("Unable to determine window sizes from the MSP object.")
        lengths = np.asarray(window_sizes, dtype=np.float64)

    if np.any(lengths <= 0):
        raise ValueError("Window lengths must be positive.")
    return lengths


def _mean_and_iqr(values: np.ndarray) -> Tuple[float, Tuple[float, float]]:
    if values is None or values.size == 0:
        return np.nan, (np.nan, np.nan)
    mean = float(np.mean(values))
    quartiles = np.percentile(values, [25, 75])
    return mean, (float(quartiles[0]), float(quartiles[1]))


def _format_interval(interval: Tuple[float, float], suffix: str = "") -> str:
    lower, upper = interval
    if np.isnan(lower) or np.isnan(upper):
        return "NA"
    return f"{lower:.2f}-{upper:.2f}{suffix}"


def _expand_msp_inputs(target: str, default_pattern: str) -> List[str]:
    files: List[str]
    if any(char in target for char in "*?[]"):
        files = sorted(glob.glob(target))
    else:
        candidate = Path(target)
        if candidate.is_dir():
            files = sorted(str(p) for p in candidate.glob(default_pattern))
        elif candidate.is_file():
            files = [str(candidate)]
        else:
            files = sorted(glob.glob(target))

    if not files:
        raise FileNotFoundError(f"No MSP files matched '{target}'.")
    return files


def build_summary_table(files: List[str]) -> pd.DataFrame:
    rows: List[Dict[str, object]] = []
    aggregator = DatasetAncestrySummary(dataset_name="UKBB")
    for path in files:
        print(f"[UKBB] Processing {path}")
        aggregator.consume_file(path)
    rows.extend(aggregator.summarize_rows())

    if not rows:
        return pd.DataFrame(
            columns=[
                "Dataset",
                "Ancestry",
                "Mean_percent_genome",
                "Percent_IQR",
                "Mean_tract_length_Mb",
                "Tract_length_IQR_Mb",
                "Mean_tracts_per_individual",
            ]
        )

    return pd.DataFrame(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create a reviewer-ready LAI summary table (mean % genome, tract lengths, tract counts)."
    )
    parser.add_argument(
        "--msp-input",
        required=True,
        help="Path to the MSP file, directory, or glob pattern containing UKBB MSP files.",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Path to the Excel file that will contain the summary table.",
    )
    parser.add_argument(
        "--sheet-name",
        default="LAI_summary",
        help="Worksheet name for the Excel output.",
    )
    parser.add_argument(
        "--msp-pattern",
        default="*.msp.tsv",
        help="Glob pattern used when a dataset points to a directory (default: *.msp.tsv).",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    files = _expand_msp_inputs(args.msp_input, args.msp_pattern)
    table = build_summary_table(files)
    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    sheet_name = f"{args.sheet_name}_ukbb"
    table.to_excel(args.output, index=False, sheet_name=sheet_name)
    print(f"Summary written to {args.output}")


if __name__ == "__main__":
    main()
