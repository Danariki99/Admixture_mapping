import os
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


def load_pheno_table(path, sheet_name):
    if not os.path.exists(path):
        raise FileNotFoundError(f"Phenotype table not found: {path}")
    return pd.read_excel(path, sheet_name=sheet_name)


def get_pheno_name(df_pheno: pd.DataFrame, pheno_id: str) -> str:
    pheno_row = df_pheno[df_pheno["ID"] == pheno_id]["ID2"]
    if not pheno_row.empty:
        return str(pheno_row.iloc[0])
    return "Unknown"

def plot_locuszoom(
    df_add: pd.DataFrame,
    snp: str,
    chromosome: str,
    title: str,
    out_png: str,
    window_bp: int = 400_000,   # matches your example PDF (~ +/- 0.4 Mb)
    figsize=(10.5, 4.2),        # close to the PDF aspect
    point_size: int = 14,
):
    """
    Minimal LocusZoom-like plot from PLINK2 .glm.logistic.hybrid output.

    Inputs:
      - df_add: already filtered to TEST == 'ADD'
      - snp: lead SNP (ID)
      - chromosome: e.g. 'chr6' (used only for labels)
      - out_png: path to save png
    """

    # --- basic sanity + required columns ---
    required = {"POS", "ID", "P"}
    missing = required - set(df_add.columns)
    if missing:
        raise ValueError(f"Missing required columns in df_add: {sorted(missing)}")

    d = df_add.copy()

    # Ensure numeric + valid p-values
    d["POS"] = pd.to_numeric(d["POS"], errors="coerce")
    d["P"] = pd.to_numeric(d["P"], errors="coerce")
    d = d.dropna(subset=["POS", "P"])
    d = d[(d["P"] > 0) & np.isfinite(d["P"])]

    # Find lead SNP position
    lead_rows = d[d["ID"].astype(str) == str(snp)]
    if len(lead_rows) == 0:
        # fallback: use the best (smallest) p-value in region
        lead_idx = d["P"].idxmin()
        lead_pos = float(d.loc[lead_idx, "POS"])
        lead_snp = str(d.loc[lead_idx, "ID"])
    else:
        lead_pos = float(lead_rows.iloc[0]["POS"])
        lead_snp = str(snp)

    # Window around lead
    lo = lead_pos - window_bp
    hi = lead_pos + window_bp
    d = d[(d["POS"] >= lo) & (d["POS"] <= hi)].copy()

    if d.empty:
        raise ValueError(f"No variants left after window filter around {lead_snp} ({chromosome}:{lead_pos}).")

    # Prepare axes values
    d["pos_mb"] = d["POS"] / 1e6
    d["mlog10p"] = -np.log10(d["P"])

    # --- plotting ---
    fig, ax = plt.subplots(figsize=figsize)

    # Other SNPs
    is_lead = d["ID"].astype(str) == str(lead_snp)
    ax.scatter(
        d.loc[~is_lead, "pos_mb"],
        d.loc[~is_lead, "mlog10p"],
        s=point_size,
        alpha=0.9,
        linewidths=0.0,
        label="Other SNPs",
    )

    # Lead SNP (highlight)
    ax.scatter(
        d.loc[is_lead, "pos_mb"],
        d.loc[is_lead, "mlog10p"],
        s=point_size * 2.2,
        alpha=1.0,
        linewidths=0.8,
        edgecolors="black",
        label=f"Lead: {lead_snp}",
        zorder=5,
    )

    # Label lead SNP
    if is_lead.any():
        x0 = float(d.loc[is_lead, "pos_mb"].iloc[0])
        y0 = float(d.loc[is_lead, "mlog10p"].iloc[0])
        ax.text(
            x0,
            y0,
            lead_snp,
            fontsize=9,
            ha="center",
            va="bottom",
        )

    ax.set_title(title)
    ax.set_xlabel(f"Position on {chromosome} (Mb)")
    ax.set_ylabel(r"$-\log_{10}(p\mathrm{-value})$")

    # Legend styling similar to locuszoom-ish
    ax.legend(frameon=True, fontsize=9, loc="upper left")

    # Nice x-limits tightly around data (still within window)
    ax.set_xlim(d["pos_mb"].min(), d["pos_mb"].max())

    # A bit of headroom
    ymax = d["mlog10p"].max()
    ax.set_ylim(0, max(1.0, ymax * 1.05))

    fig.tight_layout()

    out_dir = os.path.dirname(out_png)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    fig.savefig(out_png, dpi=300)
    plt.close(fig)


if __name__ == "__main__":
 #define all the data
    HIT_FOLDER = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/fine_mapping_ancestries_PCA_verbose'
    OUTPUT_FOLDER ='/private/groups/ioannidislab/smeriglio/out_cleaned_codes/locuszoom_plots'
    PHENO_TABLE = 'ukbb_v1.xlsx'
    PHENO_SHEET = 'first_batch'
    df_first_batch = load_pheno_table(PHENO_TABLE, PHENO_SHEET)
    counter = 0

    snp_to_info = {
  'rs72826020': [{'chromosome': 'chr10',
                  'ancestry_tested': 'WAS',
                  'ancestry_of_population': 'EUR',
                  'phenotype': 'HC221'}],
  'rs41307444': [{'chromosome': 'chr9',
                  'ancestry_tested': 'EUR',
                  'ancestry_of_population': 'AFR',
                  'phenotype': 'HC1007'}],
  'rs6931921': [{'chromosome': 'chr6',
                 'ancestry_tested': 'WAS',
                 'ancestry_of_population': 'EUR',
                 'phenotype': 'HC643'}],
  'rs887468': [{'chromosome': 'chr6',
                'ancestry_tested': 'WAS',
                'ancestry_of_population': 'WAS',
                'phenotype': 'HC643'}],
  'rs28383172': [{'chromosome': 'chr6',
                  'ancestry_tested': 'SAS',
                  'ancestry_of_population': 'EAS',
                  'phenotype': 'HC1581'},
                 {'chromosome': 'chr6',
                  'ancestry_tested': 'SAS',
                  'ancestry_of_population': 'EAS',
                  'phenotype': 'HC1036'}],
  'rs4248166': [{'chromosome': 'chr6',
                 'ancestry_tested': 'SAS',
                 'ancestry_of_population': 'EUR',
                 'phenotype': 'HC1581'},
                {'chromosome': 'chr6',
                 'ancestry_tested': 'SAS',
                 'ancestry_of_population': 'EUR',
                 'phenotype': 'HC1036'}],
  'rs2844477': [{'chromosome': 'chr6',
                 'ancestry_tested': 'SAS',
                 'ancestry_of_population': 'WAS',
                 'phenotype': 'HC1581'}],
  'rs7758128': [{'chromosome': 'chr6',
                 'ancestry_tested': 'SAS',
                 'ancestry_of_population': 'NAT',
                 'phenotype': 'HC1581'}],
  'rs4128469': [{'chromosome': 'chr8',
                 'ancestry_tested': 'EUR',
                 'ancestry_of_population': 'SAS',
                 'phenotype': 'HC221'}],
  'rs4870843': [{'chromosome': 'chr8',
                 'ancestry_tested': 'EUR',
                 'ancestry_of_population': 'NAT',
                 'phenotype': 'HC221'}],
  'rs2248372': [{'chromosome': 'chr6',
                 'ancestry_tested': 'EUR',
                 'ancestry_of_population': 'EUR',
                 'phenotype': 'HC382'},
                {'chromosome': 'chr6',
                 'ancestry_tested': 'EUR',
                 'ancestry_of_population': 'EUR',
                 'phenotype': 'HC1581'}],
  'rs3131633': [{'chromosome': 'chr6',
                 'ancestry_tested': 'EUR',
                 'ancestry_of_population': 'SAS',
                 'phenotype': 'HC382'}],
  'rs2596548': [{'chromosome': 'chr6',
                 'ancestry_tested': 'EUR',
                 'ancestry_of_population': 'NAT',
                 'phenotype': 'HC382'},
                {'chromosome': 'chr6',
                 'ancestry_tested': 'EUR',
                 'ancestry_of_population': 'NAT',
                 'phenotype': 'HC1581'}],
  'Affx-28441669': [{'chromosome': 'chr6',
                     'ancestry_tested': 'WAS',
                     'ancestry_of_population': 'EUR',
                     'phenotype': 'HC219'}],
  'rs10484554': [{'chromosome': 'chr6',
                  'ancestry_tested': 'SAS',
                  'ancestry_of_population': 'EAS',
                  'phenotype': 'HC382'}],
  'rs9266490': [{'chromosome': 'chr6',
                 'ancestry_tested': 'SAS',
                 'ancestry_of_population': 'EUR',
                 'phenotype': 'HC382'}],
  'rs4394275': [{'chromosome': 'chr6',
                 'ancestry_tested': 'SAS',
                 'ancestry_of_population': 'WAS',
                 'phenotype': 'HC382'}],
  'rs13198903': [{'chromosome': 'chr6',
                  'ancestry_tested': 'SAS',
                  'ancestry_of_population': 'NAT',
                  'phenotype': 'HC382'}],
  'rs903160': [{'chromosome': 'chr17',
                'ancestry_tested': 'SAS',
                'ancestry_of_population': 'EUR',
                'phenotype': 'FH1220'}],
  'rs62067003': [{'chromosome': 'chr17',
                  'ancestry_tested': 'SAS',
                  'ancestry_of_population': 'SAS',
                  'phenotype': 'FH1220'}],
  'rs9469093': [{'chromosome': 'chr6',
                 'ancestry_tested': 'SAS',
                 'ancestry_of_population': 'WAS',
                 'phenotype': 'HC1036'}],
  'rs379464': [{'chromosome': 'chr6',
                'ancestry_tested': 'SAS',
                'ancestry_of_population': 'NAT',
                'phenotype': 'HC1036'}],
  'rs9972241': [{'chromosome': 'chr14',
                 'ancestry_tested': 'EUR',
                 'ancestry_of_population': 'EUR',
                 'phenotype': 'HC221'}]
}

    counter = 30
    for snp, info_list in snp_to_info.items():
        
        for info in info_list:
            folder = f'{HIT_FOLDER}/{info["ancestry_tested"]}_{info["phenotype"]}_{info["chromosome"]}'
            p_file = os.path.join(folder, f'{info["ancestry_tested"]}_{info["phenotype"]}_{info["chromosome"]}_output.{info["ancestry_of_population"]}.{info["phenotype"]}.glm.logistic.hybrid')

            df = pd.read_csv(p_file, sep="\t", dtype=str, low_memory=False)

            # 2) Keep only ADD rows
            df_add = df[df["TEST"] == "ADD"].copy()

            # 3) Output name (Supplementary Figure n + cleaned phenotype name)
            pheno_name = get_pheno_name(df_first_batch, info["phenotype"])
            if "/" in pheno_name:
                pheno_name = pheno_name.replace("/", "_")
            out_png = os.path.join(
                OUTPUT_FOLDER,
                f'Supplementary Figure {counter}.png'
            )
            if pheno_name == 'HC1007_TTE_acute_upper_respiratory_infections_of_multiple_and_unspecified_sites':
                print('hereeeee')
                pheno_name = 'HC1007_TTE_acute_upper_respiratory_infections'

            # 4) Title similar to your plot naming
            title = f'{info["ancestry_tested"]} | {"_".join(pheno_name.split("_")[1:])} | {info["ancestry_of_population"]} population | {info["chromosome"]} | {snp}'
            title = title.replace('NAT', 'AMR')

            # 5) Plot
            plot_locuszoom(
                df_add=df_add,
                snp=snp,
                chromosome=info["chromosome"],
                title=title,
                out_png=out_png,
            )
            counter += 1
    print(f'Total plots generated: {counter}')
