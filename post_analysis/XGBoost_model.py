import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.metrics import (balanced_accuracy_score, roc_auc_score,
                              f1_score, precision_score, recall_score,
                              confusion_matrix, ConfusionMatrixDisplay)
from xgboost import XGBClassifier
import shap

# Paths
base_xgboost  = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/XGBoost_test_snps'
base_output   = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/XGBoost_output'
phe_folder    = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/phe_files/ukbb'
wind_folder   = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/FUMA/ukbb/wind'
covar_folder  = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/wind_covar_files_new'

NON_SNP = {'FID', 'IID', 'PAT', 'MAT', 'SEX', 'PHENOTYPE'}

# Columns always excluded (LAI is window-specific, EUR is reference category)
ALWAYS_EXCLUDE = {'FID', 'IID', 'LAI', 'EUR'}

# Ancestry proportion columns
PROP_COLS = {'AFR', 'AHG', 'EAS', 'NAT', 'OCE', 'SAS', 'WAS'}

XGB_PARAMS = dict(
    n_estimators=300,
    max_depth=4,
    learning_rate=0.05,
    subsample=0.8,
    colsample_bytree=0.8,
    random_state=42,
    eval_metric='logloss',
    verbosity=0,
)


def run_xgboost(X, y, feature_cols, label, output_folder, hit_label):
    """Train 80/20, evaluate, save confusion matrix, retrain full + SHAP.
    Returns a dict of metrics for CSV summary."""

    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, random_state=42, stratify=y
    )
    n_train_cases    = int(y_train.sum())
    n_train_controls = int((y_train == 0).sum())
    n_test_cases     = int(y_test.sum())
    n_test_controls  = int((y_test == 0).sum())

    print(f"\n  [{label}] Train: {len(y_train):,} ({n_train_cases:,} cases) | "
          f"Test: {len(y_test):,} ({n_test_cases:,} cases)")

    model = XGBClassifier(**XGB_PARAMS)
    model.fit(X_train, y_train, verbose=False)

    y_pred = model.predict(X_test)
    y_prob = model.predict_proba(X_test)[:, 1]

    auc  = roc_auc_score(y_test, y_prob)
    bac  = balanced_accuracy_score(y_test, y_pred)
    f1   = f1_score(y_test, y_pred, zero_division=0)
    prec = precision_score(y_test, y_pred, zero_division=0)
    rec  = recall_score(y_test, y_pred, zero_division=0)

    print(f"  [{label}] Test-set metrics:")
    print(f"    {'balanced_accuracy':<22} {bac:.4f}")
    print(f"    {'AUC':<22} {auc:.4f}")
    print(f"    {'F1':<22} {f1:.4f}")
    print(f"    {'precision':<22} {prec:.4f}")
    print(f"    {'recall':<22} {rec:.4f}")

    # Confusion matrix
    cm = confusion_matrix(y_test, y_pred)
    fig, ax = plt.subplots(figsize=(5, 4))
    ConfusionMatrixDisplay(confusion_matrix=cm,
                           display_labels=['Control', 'Case']).plot(
        ax=ax, colorbar=False, cmap='Blues'
    )
    ax.set_title(f'{hit_label}\n{label}  —  Test set (20%)  AUC={auc:.3f}')
    plt.tight_layout()
    cm_path = os.path.join(output_folder, f'{hit_label}_{label}_confusion_matrix.png')
    plt.savefig(cm_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  [{label}] Confusion matrix → {cm_path}")

    # Retrain on full data → SHAP
    print(f"  [{label}] Training full model for SHAP...")
    model_full = XGBClassifier(**XGB_PARAMS)
    model_full.fit(X, y, verbose=False)

    explainer = shap.TreeExplainer(model_full)
    shap_vals = explainer.shap_values(X)

    max_display = min(30, len(feature_cols))
    plt.figure(figsize=(10, max(4, max_display * 0.35)))
    shap.summary_plot(shap_vals, X, feature_names=feature_cols,
                      show=False, max_display=max_display)
    shap_path = os.path.join(output_folder, f'{hit_label}_{label}_shap.png')
    plt.savefig(shap_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  [{label}] SHAP summary     → {shap_path}")

    return {
        'n_train':          len(y_train),
        'n_train_cases':    n_train_cases,
        'n_train_controls': n_train_controls,
        'n_test':           len(y_test),
        'n_test_cases':     n_test_cases,
        'n_test_controls':  n_test_controls,
        'balanced_accuracy': round(bac,  4),
        'AUC':               round(auc,  4),
        'F1':                round(f1,   4),
        'precision':         round(prec, 4),
        'recall':            round(rec,  4),
    }


# ── Main loop: iterate over all hits ──────────────────────────────────────────
all_results = []
for wind_filename in sorted(os.listdir(wind_folder)):
    if not wind_filename.endswith('_wind.txt'):
        continue

    ancestry = wind_filename.split('_')[0]
    pheno    = wind_filename.split('_')[1]

    # Load phenotype
    phe_file = os.path.join(phe_folder, f'{pheno}.phe')
    if not os.path.exists(phe_file):
        print(f"Skipping {wind_filename}: phenotype file not found")
        continue

    phe_df = pd.read_csv(phe_file, sep=r'\s+')
    phe_df.columns = [c.lstrip('#') for c in phe_df.columns]
    phe_df['IID']  = phe_df['IID'].astype(str)
    pheno_col = phe_df.columns[-1]
    phe_df    = phe_df[phe_df[pheno_col].isin([1, 2])].copy()
    phe_df['y'] = phe_df[pheno_col].astype(int) - 1

    wind_file = os.path.join(wind_folder, wind_filename)
    wind_df   = pd.read_csv(wind_file, sep='\t')

    for chr_val in wind_df['chr'].unique():
        wind_chr = wind_df[wind_df['chr'] == chr_val]
        hit_label = f'{ancestry}_{pheno}_chr{chr_val}'

        input_folder  = os.path.join(base_xgboost, hit_label)
        output_folder = os.path.join(base_output,  hit_label)
        os.makedirs(output_folder, exist_ok=True)

        print(f"\n{'='*65}")
        print(f"{hit_label}")

        # ── Load single combined .raw file ─────────────────────────────────────
        raw_file = os.path.join(input_folder, f'{hit_label}_snps.raw')
        if not os.path.exists(raw_file):
            print(f"  Skipping: raw file not found ({raw_file})")
            continue

        raw_df = pd.read_csv(raw_file, sep='\t')
        raw_df.columns = [c.lstrip('#') for c in raw_df.columns]
        raw_df['IID']  = raw_df['IID'].astype(str)
        snp_cols = [c for c in raw_df.columns if c.upper() not in NON_SNP]
        print(f"  SNPs: {len(snp_cols)}")

        # ── Load covariates from first available window covar file ─────────────
        covar_df = None
        for _, row in wind_chr.iterrows():
            cf = os.path.join(
                covar_folder,
                f'{ancestry}_{pheno}_chr{chr_val}_{int(row["start"])}_{int(row["end"])}_covar.tsv'
            )
            if os.path.exists(cf):
                covar_df = pd.read_csv(cf, sep='\t')
                covar_df.columns = [c.lstrip('#') for c in covar_df.columns]
                covar_df['IID']   = covar_df['IID'].astype(str)
                break

        if covar_df is None:
            print(f"  Skipping: no covariate file found")
            continue

        # All covar columns except always-excluded ones
        all_covar_cols  = [c for c in covar_df.columns if c.upper() not in ALWAYS_EXCLUDE]
        # Covar columns without ancestry proportions
        basic_covar_cols = [c for c in all_covar_cols if c.upper() not in PROP_COLS]

        # ── Merge SNPs + all covariates + ancestry col for ranking + phenotype ─
        # Always load the focal ancestry column for balancing
        balance_col  = [ancestry] if ancestry not in all_covar_cols and ancestry in covar_df.columns else []
        load_cols    = list(set(all_covar_cols + balance_col))

        merged = (raw_df[['IID'] + snp_cols]
                  .merge(covar_df[['IID'] + load_cols], on='IID', how='inner')
                  .merge(phe_df[['IID', 'y']], on='IID', how='inner'))

        n_cases_orig    = int((merged['y'] == 1).sum())
        n_controls_orig = int((merged['y'] == 0).sum())
        print(f"  Before balancing — Cases: {n_cases_orig:,}  |  Controls: {n_controls_orig:,}")

        # ── Balance: keep all cases, top N controls by focal ancestry prop ─────
        cases    = merged[merged['y'] == 1]
        controls = merged[merged['y'] == 0]
        controls_ranked = controls.sort_values(ancestry, ascending=False).head(len(cases))
        merged_bal = pd.concat([cases, controls_ranked]).sample(frac=1, random_state=42)

        n_bal = int((merged_bal['y'] == 1).sum())
        print(f"  After balancing  — Cases: {n_bal:,}  |  Controls (top {ancestry}): {n_bal:,}")

        # ── Prepare X for both model variants ─────────────────────────────────
        def build_X(df, feat_cols):
            X = df[feat_cols].replace('NA', np.nan).values.astype(float)
            col_means = np.nanmean(X, axis=0)
            nan_mask  = np.isnan(X)
            X[nan_mask] = np.take(col_means, np.where(nan_mask)[1])
            return X

        y = merged_bal['y'].values.astype(int)

        base_row = {
            'hit':                hit_label,
            'ancestry':           ancestry,
            'pheno':              pheno,
            'chr':                chr_val,
            'n_snps':             len(snp_cols),
            'n_total_balanced':   len(merged_bal),
            'n_cases_balanced':   n_bal,
            'n_controls_balanced': n_bal,
        }

        # Model A: SNPs + all ancestry proportions (except EUR) + age/sex/BMI
        feat_A = snp_cols + all_covar_cols
        X_A    = build_X(merged_bal, feat_A)
        print(f"\n  Model with_prop  — {len(feat_A)} features ({len(snp_cols)} SNPs + {len(all_covar_cols)} covariates)")
        metrics_A = run_xgboost(X_A, y, feat_A, 'with_prop', output_folder, hit_label)

        # Model B: SNPs + age/sex/BMI only (no ancestry proportions)
        feat_B = snp_cols + basic_covar_cols
        X_B    = build_X(merged_bal, feat_B)
        print(f"\n  Model no_prop    — {len(feat_B)} features ({len(snp_cols)} SNPs + {len(basic_covar_cols)} covariates)")
        metrics_B = run_xgboost(X_B, y, feat_B, 'no_prop', output_folder, hit_label)

        # One row per experiment (long format)
        for label, metrics in [('with_prop', metrics_A), ('no_prop', metrics_B)]:
            row = {**base_row, 'experiment': label, **metrics}
            all_results.append(row)

        # Save per-hit CSV (2 rows: one per experiment)
        hit_csv = os.path.join(output_folder, f'{hit_label}_results.csv')
        hit_rows = [r for r in all_results if r['hit'] == hit_label]
        pd.DataFrame(hit_rows).to_csv(hit_csv, index=False)
        print(f"\n  Results saved → {hit_csv}")

print(f"\n{'='*65}")
if all_results:
    summary_path = os.path.join(base_output, 'XGBoost_summary.csv')
    pd.DataFrame(all_results).to_csv(summary_path, index=False)
    print(f"Summary saved → {summary_path}")
print("Done.")
