#!/usr/bin/env python3
# ==============================================================================
# audit_data.py — Auditoria reprodutível de integridade de dados (pré-commit)
#
# Reproduz as checagens descritas em AUDIT_REPORT.md:
#   1. Tipos mistos / colunas-índice residuais
#   2. Missing (% por variável, alerta > 5%)
#   3. Duplicatas nas tabelas primárias
#   4. Consistência amostra × metadado (RPKM vs pheno)
#   5. Benchmarks: FC do NQO1, FC da doxorrubicina, outliers (Z-score), normalização
#
# Uso:
#   python scripts/audit_data.py            # imprime relatório no stdout
#   python scripts/audit_data.py --json      # grava output/audit/audit_metrics.json
# ==============================================================================

from __future__ import annotations

import argparse
import gzip
import json
import os
import re
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]

# ----------------------------------------------------------------------------
# Configuração (espelha config.yaml; idealmente ler com yaml.safe_load)
# ----------------------------------------------------------------------------
CFG = {
    "missing_threshold": 0.05,
    "z_score_threshold": 3.0,
    "fold_change_significant": 1.5,
    "log2fc_significant": 0.585,  # log2(1.5)
}


def pct_na(df: pd.DataFrame) -> pd.Series:
    return df.isna().mean().mul(100).round(2)


def z_outliers(s: pd.Series, thr: float = 3.0) -> int:
    s = pd.to_numeric(s, errors="coerce")
    if s.nunique() <= 1:
        return 0
    return int((((s - s.mean()) / s.std(ddof=0)).abs() > thr).sum())


def classify_sample(name: str) -> str:
    m = re.match(r"([A-Za-z]+)", str(name))
    return m.group(1) if m else str(name)


def main() -> dict:
    report: dict = {}

    # ------------------------------------------------------------------ 1. NQO1
    dgb_path = ROOT / "output" / "tables" / "nqo1" / "DGB_results_NQO1_from_xlsx_clean.csv"
    dgb = pd.read_csv(dgb_path, low_memory=False)

    report["nqo1_clean"] = {
        "shape": list(dgb.shape),
        "duplicated_rows": int(dgb.duplicated().sum()),
        "missing_pct_gt5": pct_na(dgb)[pct_na(dgb) > 5].round(2).to_dict(),
        "index_col_X1_present": "X1" in dgb.columns,
        "index_col_X1_unique": int(dgb["X1"].nunique()) if "X1" in dgb.columns else 0,
        "index_col_X1_duplicated": int(dgb["X1"].duplicated().sum()) if "X1" in dgb.columns else 0,
        "q_value_zero": int((dgb["q_value"] == 0).sum()),
        "p_value_zero": int((dgb["p_value"] == 0).sum()),
        "has_effect_score": "effect_score" in dgb.columns,
        "has_log2FC_derived": "log2FC_derived" in dgb.columns,
        "has_q_value_raw": "q_value_raw" in dgb.columns,
    }

    doxo = dgb[dgb["Drug_Name_Norm"].str.contains("doxo", na=False)]
    if not doxo.empty:
        report["doxorubicin_nqo1"] = {
            "n_tests": int(len(doxo)),
            "mean_fc": float(doxo["Fold_Change"].mean()),
            "median_fc": float(doxo["Fold_Change"].median()),
            "max_abs_fc": float(doxo["Fold_Change"].abs().max()),
            "min_q": float(doxo["q_value"].min()),
        }

    # -------------------------------------------------------------- 2. DEG NQO1
    deg = pd.read_csv(ROOT / "output" / "tables" / "ampk" / "DEG_full_table.csv")
    nqo1 = deg[deg["gene"] == "NQO1"]
    report["deg_full"] = {
        "shape": list(deg.shape),
        "duplicated_rows": int(deg.duplicated().sum()),
        "duplicated_genes": int(deg["gene"].duplicated().sum()),
        "significant": int(deg["significant"].sum()),
    }
    if not nqo1.empty:
        r = nqo1.iloc[0]
        report["nqo1_deg"] = {
            "log2FoldChange": float(r["log2FoldChange"]),
            "padj": float(r["padj"]),
            "direction": str(r["direction"]),
            "passes_fc_threshold": abs(float(r["log2FoldChange"])) > CFG["log2fc_significant"],
        }

    # ------------------------------------------------------------ 3. Z-outliers
    z = {}
    for col in ["p_value", "q_value", "Fold_Change", "Specificity", "score", "neg_log10_q_capped"]:
        if col in dgb.columns:
            z[col] = z_outliers(dgb[col], CFG["z_score_threshold"])
    report["z_outliers"] = z

    # ---------------------------------------------------- 4. RPKM missing + NQO1
    rpkm_path = ROOT / "data" / "raw" / "GSE116250_rpkm.txt.gz"
    with gzip.open(rpkm_path, "rt") as f:
        header = f.readline().strip().split("\t")
        samples = header[2:]
        n_missing = 0
        total_cells = 0
        n_genes = 0
        nqo1_vals = None
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            n_genes += 1
            vals = parts[2:]
            n_missing += sum(1 for v in vals if v.strip() == "")
            total_cells += len(vals)
            if parts[1] == "NQO1":
                nqo1_vals = [float(v) for v in vals if v.strip() != ""]

    report["rpkm"] = {
        "n_samples": len(samples),
        "class_counts": dict(Counter(classify_sample(s) for s in samples)),
        "n_genes": n_genes,
        "missing_cells": n_missing,
        "total_cells": total_cells,
        "missing_pct": round(100 * n_missing / total_cells, 4) if total_cells else None,
        "nqo1_present": nqo1_vals is not None,
    }
    if nqo1_vals:
        report["rpkm"]["nqo1_mean"] = round(float(np.mean(nqo1_vals)), 3)
        report["rpkm"]["nqo1_median"] = round(float(np.median(nqo1_vals)), 3)

    # --------------------------------------------- 5. Consistência amostra×metadado
    pheno = pd.read_csv(ROOT / "output" / "tables" / "ampk" / "pheno_data_raw.csv")
    rpk_set = set(samples)
    pheno_set = set(pheno["title"].astype(str))
    report["sample_label_consistency"] = {
        "rpk_n": len(rpk_set),
        "pheno_n": len(pheno_set),
        "in_rpk_not_pheno": sorted(rpk_set - pheno_set),
        "in_pheno_not_rpk": sorted(pheno_set - rpk_set),
        "disease_counts": pheno["disease:ch1"].value_counts().to_dict(),
    }

    # ---------------------------------------------------------- 6. Duplicatas PPI
    edges = pd.read_csv(ROOT / "output" / "tables" / "ppi" / "PPI_AMPK_STRING_edges.csv")
    nodes = pd.read_csv(ROOT / "output" / "tables" / "ppi" / "PPI_AMPK_STRING_nodes.csv")
    report["ppi"] = {
        "edges_shape": list(edges.shape),
        "edges_dup": int(edges.duplicated().sum()),
        "nodes_shape": list(nodes.shape),
        "nodes_dup": int(nodes.duplicated().sum()),
    }

    # ----------------------------------------------------------- 7. GSEA (AMPK)
    gsea = pd.read_csv(ROOT / "output" / "tables" / "ampk" / "GSEA_AMPK_results.csv")
    if not gsea.empty:
        report["gsea_ampk"] = {
            "NES": float(gsea["NES"].iloc[0]),
            "padj": float(gsea["padj"].iloc[0]),
            "significant": bool(float(gsea["padj"].iloc[0]) < CFG["missing_threshold"]),
        }

    # ------------------------------------------- 8. Integridade dos outputs finais
    ampk_dir = ROOT / "output" / "tables" / "ampk"
    nqo1_dir = ROOT / "output" / "tables" / "nqo1"

    def _csv_stats(path):
        if not path.exists():
            return {"present": False}
        try:
            df = pd.read_csv(path)
            return {"present": True, "rows": int(len(df)), "cols": int(df.shape[1])}
        except Exception as e:  # noqa: BLE001
            return {"present": True, "error": str(e)}

    gsva_path = ampk_dir / "GSVA_AMPK_results.csv"

    report["outputs_integrity"] = {
        "DEG_full_table.csv": _csv_stats(ampk_dir / "DEG_full_table.csv"),
        "NQO1_results.csv": _csv_stats(ampk_dir / "NQO1_results.csv"),
        "DEG_AMPK_pathway.csv": _csv_stats(ampk_dir / "DEG_AMPK_pathway.csv"),
        "GSEA_AMPK_results.csv": _csv_stats(ampk_dir / "GSEA_AMPK_results.csv"),
        "GSVA_AMPK_results.csv": _csv_stats(gsva_path),
        "DGB_results_NQO1.xlsx": {"present": (nqo1_dir / "DGB_results_NQO1.xlsx").exists()},
        "doxorubicin_NQO1_summary.csv": _csv_stats(nqo1_dir / "doxorubicin_NQO1_summary.csv"),
    }

    return report


def print_report(rep: dict) -> None:
    def sec(title: str) -> None:
        print(f"\n=== {title} ===")

    sec("1. NQO1 clean CSV")
    print(json.dumps(rep["nqo1_clean"], indent=2, ensure_ascii=False))
    sec("2. Doxorrubicina -> NQO1")
    print(json.dumps(rep.get("doxorubicin_nqo1", {}), indent=2, ensure_ascii=False))
    sec("3. DEG NQO1")
    print(json.dumps(rep.get("nqo1_deg", {}), indent=2, ensure_ascii=False))
    print(json.dumps(rep.get("deg_full", {}), indent=2, ensure_ascii=False))
    sec("4. Outliers (|Z|>3)")
    print(json.dumps(rep.get("z_outliers", {}), indent=2))
    sec("5. RPKM")
    print(json.dumps(rep.get("rpkm", {}), indent=2, ensure_ascii=False))
    sec("6. Consistência amostra x metadado")
    print(json.dumps(rep.get("sample_label_consistency", {}), indent=2, ensure_ascii=False))
    sec("7. PPI")
    print(json.dumps(rep.get("ppi", {}), indent=2))
    sec("8. GSEA AMPK")
    print(json.dumps(rep.get("gsea_ampk", {}), indent=2))
    sec("9. Integridade dos outputs finais")
    print(json.dumps(rep.get("outputs_integrity", {}), indent=2, ensure_ascii=False))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", action="store_true", help="grava output/audit/audit_metrics.json")
    args = ap.parse_args()

    rep = main()

    if args.json:
        out_dir = ROOT / "output" / "audit"
        out_dir.mkdir(parents=True, exist_ok=True)
        out = out_dir / "audit_metrics.json"
        out.write_text(json.dumps(rep, indent=2, ensure_ascii=False), encoding="utf-8")
        print(f"Gravado: {out}")

    print_report(rep)
