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
base_xgboost = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/XGBoost_test_snps'
phe_folder   = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/phe_files/ukbb'
wind_folder  = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/FUMA/ukbb/wind'
covar_folder = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/wind_covar_files_new'

# Hit to analyze
wind_filename = 'EAS_HC1158_wind.txt'
ancestry = wind_filename.split('_')[0]   # EAS
pheno    = wind_filename.split('_')[1]   # HC1158

# ── Load phenotype ─────────────────────────────────────────────────────────────
phe_file = os.path.join(phe_folder, f'{pheno}.phe')
phe_df   = pd.read_csv(phe_file, sep='\t')
phe_df.columns = [c.lstrip('#') for c in phe_df.columns]
phe_df['IID']  = phe_df['IID'].astype(str)

pheno_col = phe_df.columns[-1]
phe_df    = phe_df[phe_df[pheno_col].isin([1, 2])].copy()
phe_df['y'] = phe_df[pheno_col].astype(int) - 1   # 1/2 → 0/1

# ── Windows ────────────────────────────────────────────────────────────────────
wind_file = os.path.join(wind_folder, wind_filename)
wind_df   = pd.read_csv(wind_file, sep='\t')

NON_SNP   = {'FID', 'IID', 'PAT', 'MAT', 'SEX', 'PHENOTYPE'}
NON_COVAR = {'FID', 'IID'}

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

for chr_val in wind_df['chr'].unique():
    wind_chr     = wind_df[wind_df['chr'] == chr_val]
    input_folder = os.path.join(base_xgboost, f'{ancestry}_{pheno}_chr{chr_val}')

    for _, row in wind_chr.iterrows():
        window_start = int(row['start'])
        window_end   = int(row['end'])

        raw_file   = os.path.join(
            input_folder,
            f'{ancestry}_{pheno}_chr{chr_val}_{window_start}_{window_end}_snps.raw'
        )
        covar_file = os.path.join(
            covar_folder,
            f'{ancestry}_{pheno}_chr{chr_val}_{window_start}_{window_end}_covar.tsv'
        )

        if not os.path.exists(raw_file):
            print(f"Skipping (raw not found): {raw_file}")
            continue
        if not os.path.exists(covar_file):
            print(f"Skipping (covar not found): {covar_file}")
            continue

        # ── Load SNP dosages ───────────────────────────────────────────────────
        raw_df = pd.read_csv(raw_file, sep='\t')
        raw_df.columns = [c.lstrip('#') for c in raw_df.columns]
        raw_df['IID']  = raw_df['IID'].astype(str)
        snp_cols = [c for c in raw_df.columns if c.upper() not in NON_SNP]

        # ── Load covariates (age, sex, BMI, 8 ancestry props, LAI) ────────────
        covar_df = pd.read_csv(covar_file, sep='\t')
        covar_df.columns = [c.lstrip('#') for c in covar_df.columns]
        covar_df['IID']  = covar_df['IID'].astype(str)
        covar_cols = [c for c in covar_df.columns if c.upper() not in NON_COVAR]

        # ── Merge SNPs + covariates + phenotype ────────────────────────────────
        merged = (raw_df[['IID'] + snp_cols]
                  .merge(covar_df[['IID'] + covar_cols], on='IID', how='inner')
                  .merge(phe_df[['IID', 'y']], on='IID', how='inner'))

        feature_cols = snp_cols + covar_cols
        X = merged[feature_cols].replace('NA', np.nan).values.astype(float)
        y = merged['y'].values.astype(int)

        # Impute missing with column mean
        col_means = np.nanmean(X, axis=0)
        nan_mask  = np.isnan(X)
        X[nan_mask] = np.take(col_means, np.where(nan_mask)[1])

        n_cases    = int(y.sum())
        n_controls = int((y == 0).sum())

        print(f"\n{'='*65}")
        print(f"{ancestry}_{pheno}  |  chr{chr_val}:{window_start}-{window_end}")
        print(f"  Samples: {len(y):,}  |  Cases: {n_cases:,}  |  Controls: {n_controls:,}")
        print(f"  SNPs: {len(snp_cols)}  |  Covariates: {len(covar_cols)}")

        if n_cases < 10 or n_controls < 10:
            print("  Too few samples in one class — skipping")
            continue

        spw = n_controls / n_cases   # class imbalance weight

        # ── Stratified 80/20 split ─────────────────────────────────────────────
        X_train, X_test, y_train, y_test = train_test_split(
            X, y, test_size=0.2, random_state=42, stratify=y
        )
        print(f"  Train: {len(y_train):,} ({y_train.sum():,} cases) | "
              f"Test: {len(y_test):,} ({y_test.sum():,} cases)")

        # ── Train on 80% ───────────────────────────────────────────────────────
        model = XGBClassifier(**XGB_PARAMS, scale_pos_weight=spw)
        model.fit(X_train, y_train, verbose=False)

        # ── Evaluate on 20% ────────────────────────────────────────────────────
        y_pred = model.predict(X_test)
        y_prob = model.predict_proba(X_test)[:, 1]

        print(f"\n  Test-set metrics:")
        print(f"    {'balanced_accuracy':<22} {balanced_accuracy_score(y_test, y_pred):.4f}")
        print(f"    {'AUC':<22} {roc_auc_score(y_test, y_prob):.4f}")
        print(f"    {'F1':<22} {f1_score(y_test, y_pred, zero_division=0):.4f}")
        print(f"    {'precision':<22} {precision_score(y_test, y_pred, zero_division=0):.4f}")
        print(f"    {'recall':<22} {recall_score(y_test, y_pred, zero_division=0):.4f}")

        # ── Confusion matrix (test set) ────────────────────────────────────────
        cm = confusion_matrix(y_test, y_pred)
        fig, ax = plt.subplots(figsize=(5, 4))
        ConfusionMatrixDisplay(confusion_matrix=cm,
                               display_labels=['Control', 'Case']).plot(
            ax=ax, colorbar=False, cmap='Blues'
        )
        ax.set_title(
            f'{ancestry}_{pheno}  chr{chr_val}:{window_start}-{window_end}\n'
            f'Test set (20%,  stratified)'
        )
        plt.tight_layout()
        cm_path = os.path.join(
            input_folder,
            f'{ancestry}_{pheno}_chr{chr_val}_{window_start}_{window_end}_confusion_matrix.png'
        )
        plt.savefig(cm_path, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"\n  Confusion matrix  → {cm_path}")

        # ── Retrain on full data for SHAP ──────────────────────────────────────
        print(f"  Training full model for SHAP...")
        model_full = XGBClassifier(**XGB_PARAMS, scale_pos_weight=spw)
        model_full.fit(X, y, verbose=False)

        explainer = shap.TreeExplainer(model_full)
        shap_vals = explainer.shap_values(X)

        max_display = min(30, len(feature_cols))
        plt.figure(figsize=(10, max(4, max_display * 0.35)))
        shap.summary_plot(shap_vals, X, feature_names=feature_cols,
                          show=False, max_display=max_display)
        shap_path = os.path.join(
            input_folder,
            f'{ancestry}_{pheno}_chr{chr_val}_{window_start}_{window_end}_shap.png'
        )
        plt.savefig(shap_path, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"  SHAP summary      → {shap_path}")

print(f"\n{'='*65}")
print("Done.")
