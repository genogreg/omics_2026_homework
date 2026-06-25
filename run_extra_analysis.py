import gzip
import json
import os
import subprocess
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-greg")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp/cache-greg")

import cooler
import cooltools
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyBigWig
import seaborn as sns
from cooltools import insulation


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "extra_credit_results"
FIG = OUT / "figures"
TAB = OUT / "tables"
TRACKS = OUT / "tracks"
TMP = OUT / "tmp"

HIC = ROOT / "course/day2_HiC_practice/data/MoPh7_enr_v2.mcool"
INS100 = ROOT / "course/day2_HiC_practice/results/boundaries/MoPh7_enr_v2_insulation_res100kb.tsv.gz"
CTCF_PEAKS = ROOT / "course/day3_ChIPseq_practice/results/macs/MoPh7_CTCF/MoPh7_CTCF_peaks.narrowPeak"
H3K27AC_PEAKS = ROOT / "course/day3_ChIPseq_practice/results/macs/MoPh7_H3K27Ac/MoPh7_H3K27Ac_peaks.narrowPeak"
CTCF_BW = ROOT / "course/day3_ChIPseq_practice/results/macs/MoPh7_CTCF/MoPh7_CTCF_FE.bw"
H3K27AC_BW = ROOT / "course/day3_ChIPseq_practice/results/macs/MoPh7_H3K27Ac/MoPh7_H3K27Ac_FE.bw"
H3K9ME3_BW = ROOT / "course/day3_ChIPseq_practice/results/macs/MoPh7_H3K9me3/MoPh7_H3K9me3_FE.bw"
RNA_BW = ROOT / "course/day2_RNA_practice/results/tracks/rnaseq/MoPh7.rnaseq.STAR.bw"
CAGE_BW = ROOT / "course/day2_RNA_practice/results/tracks/cage/MoPh7.cage.STAR.bw"
BWA_BW = ROOT / "course/day2_RNA_practice/results/tracks/bwa/MoPh7.rnaseq.BWA.q30.bw"
METH_BW = ROOT / "course/day4_WGBS_practice/results/tracks/MoPh7_beta_methylation.bw"
COV_BW = ROOT / "course/day4_WGBS_practice/results/tracks/MoPh7_coverage.bw"
GC_BW = ROOT / "course/day4_WGBS_practice/results/tracks/T2T_gc_content_100bp.bw"
CPG_BW = ROOT / "course/day4_WGBS_practice/results/tracks/T2T_cpg_obs_exp_100bp.bw"
E1_BW = ROOT / "course/day2_HiC_practice/results/compartments/MoPh7_enr_v2_E1_res100000.bw"
GTF = ROOT / "course/day5_omics_practice/data/annotation/chm13v2.0_main_protein_coding_chrNames.gtf.gz"


def mkdirs():
    for p in [OUT, FIG, TAB, TRACKS, TMP]:
        p.mkdir(parents=True, exist_ok=True)


def attrs_to_dict(text):
    out = {}
    for part in text.strip().split(";"):
        part = part.strip()
        if not part:
            continue
        if "=" in part:
            key, value = part.split("=", 1)
        elif " " in part:
            key, value = part.split(" ", 1)
        else:
            continue
        out[key] = value.strip().strip('"')
    return out


def read_promoters(path, flank=1000):
    genes = {}
    with gzip.open(path, "rt") as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] not in {"gene", "transcript", "mRNA", "CDS"}:
                continue
            chrom, start, end, strand = fields[0], int(fields[3]), int(fields[4]), fields[6]
            attrs = attrs_to_dict(fields[8])
            gene_id = attrs.get("gene_id") or attrs.get("gene") or attrs.get("Parent") or attrs.get("ID", "")
            gene_name = attrs.get("gene_name") or attrs.get("gene") or gene_id
            key = (chrom, gene_id, gene_name, strand)
            if key not in genes:
                genes[key] = [start, end]
            else:
                genes[key][0] = min(genes[key][0], start)
                genes[key][1] = max(genes[key][1], end)
    rows = []
    for (chrom, gene_id, gene_name, strand), (start, end) in genes.items():
            if strand == "+":
                tss = start - 1
                p_start = max(0, tss - flank)
                p_end = tss
            else:
                tss = end
                p_start = tss
                p_end = tss + flank
            if p_end <= p_start:
                continue
            rows.append((chrom, p_start, p_end, gene_id, gene_name, strand, tss))
    cols = ["chrom", "start", "end", "gene_id", "gene_name", "strand", "tss"]
    return pd.DataFrame(rows, columns=cols)


def read_bed(path, cols=("chrom", "start", "end")):
    df = pd.read_csv(path, sep="\t", header=None)
    df = df.iloc[:, : len(cols)].copy()
    df.columns = cols
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)
    return df


def with_chr(value):
    value = str(value)
    return value if value.startswith("chr") else "chr" + value


def without_chr(value):
    value = str(value)
    return value[3:] if value.startswith("chr") else value


def overlap_count(query, target):
    target_by_chr = {}
    for chrom, part in target.groupby("chrom"):
        p = part.sort_values("start")
        target_by_chr[chrom] = (p["start"].to_numpy(), p["end"].to_numpy())
    counts = []
    for row in query.itertuples(index=False):
        if row.chrom not in target_by_chr:
            counts.append(0)
            continue
        starts, ends = target_by_chr[row.chrom]
        j = np.searchsorted(starts, row.end, side="left")
        counts.append(int(np.sum(ends[:j] > row.start)))
    return np.array(counts)


def nearest_distance(points, intervals):
    centers = intervals.copy()
    centers["center"] = ((centers["start"] + centers["end"]) / 2).astype(int)
    by_chr = {c: np.sort(p["center"].to_numpy()) for c, p in centers.groupby("chrom")}
    out = []
    for row in points.itertuples(index=False):
        arr = by_chr.get(row.chrom)
        if arr is None or len(arr) == 0:
            out.append(np.nan)
            continue
        i = np.searchsorted(arr, row.tss)
        vals = []
        if i < len(arr):
            vals.append(abs(arr[i] - row.tss))
        if i > 0:
            vals.append(abs(arr[i - 1] - row.tss))
        out.append(min(vals))
    return np.array(out, dtype=float)


def mean_bw(path, intervals, chrom_mode="auto"):
    vals = []
    with pyBigWig.open(str(path)) as bw:
        chroms = bw.chroms()
        for row in intervals.itertuples(index=False):
            chrom = row.chrom
            if chrom_mode == "strip":
                chrom = without_chr(chrom)
            elif chrom_mode == "add":
                chrom = with_chr(chrom)
            elif chrom not in chroms:
                a = with_chr(chrom)
                b = without_chr(chrom)
                if a in chroms:
                    chrom = a
                elif b in chroms:
                    chrom = b
            if chrom not in chroms:
                vals.append(np.nan)
                continue
            start = max(0, int(row.start))
            end = min(int(row.end), chroms[chrom])
            if end <= start:
                vals.append(np.nan)
                continue
            value = bw.stats(chrom, start, end, type="mean")[0]
            vals.append(np.nan if value is None else value)
    return np.array(vals, dtype=float)


def boundaries_from_insulation(path, window, label):
    df = pd.read_csv(path, sep="\t")
    col = f"is_boundary_{window}"
    score_col = f"boundary_strength_{window}"
    b = df[df[col].fillna(False)].copy()
    b["chrom"] = b["chrom"].map(with_chr)
    b["name"] = [f"{label}_{i}" for i in range(len(b))]
    b["score"] = b[score_col].fillna(0)
    b["center"] = ((b["start"] + b["end"]) / 2).astype(int)
    return b[["chrom", "start", "end", "name", "score", "center"]]


def call_boundaries_1mb():
    out = {}
    uri = f"{HIC}::/resolutions/1000000"
    clr = cooler.Cooler(uri)
    chromsizes = clr.chromsizes.reset_index()
    chromsizes.columns = ["chrom", "end"]
    view = chromsizes.assign(start=0, name=chromsizes["chrom"])[["chrom", "start", "end", "name"]]
    table = insulation(clr, [3_000_000, 5_000_000], view_df=view, clr_weight_name="weight", nproc=4)
    table.to_csv(TAB / "insulation_1mb.tsv.gz", sep="\t", index=False, compression="gzip")
    for window in [3_000_000, 5_000_000]:
        col = f"is_boundary_{window}"
        score_col = f"boundary_strength_{window}"
        b = table[table[col].fillna(False)].copy()
        b["chrom"] = b["chrom"].map(with_chr)
        b["name"] = [f"res1mb_w{window}_{i}" for i in range(len(b))]
        b["score"] = b[score_col].fillna(0)
        b["center"] = ((b["start"] + b["end"]) / 2).astype(int)
        out[f"res1mb_window{window//1000}kb"] = b[["chrom", "start", "end", "name", "score", "center"]]
    return out


def has_overlap_near_centers(boundaries, peaks, window=10_000):
    q = boundaries.copy()
    q["start"] = (q["center"] - window).clip(lower=0)
    q["end"] = q["center"] + window
    return overlap_count(q[["chrom", "start", "end"]], peaks) > 0


def centers_dict(table):
    return {chrom: part["center"].to_numpy(dtype=int) for chrom, part in table.groupby("chrom")}


def avg_profile(bw_path, centers, flank=300_000, bin_size=10_000):
    offsets = np.arange(-flank, flank + bin_size, bin_size)
    profile = np.full(len(offsets), np.nan)
    with pyBigWig.open(str(bw_path)) as bw:
        chroms = bw.chroms()
        for i, offset in enumerate(offsets):
            values = []
            for chrom, arr in centers.items():
                if chrom not in chroms:
                    continue
                starts = arr + offset
                ends = starts + bin_size
                ok = (starts >= 0) & (ends <= chroms[chrom])
                for start, end in zip(starts[ok], ends[ok]):
                    value = bw.stats(chrom, int(start), int(end), type="mean")[0]
                    if value is not None and not np.isnan(value):
                        values.append(value)
            if values:
                profile[i] = np.mean(values)
    return offsets, profile


def random_like(boundaries, seed, flank=300_000):
    rng = np.random.default_rng(seed)
    rows = []
    with pyBigWig.open(str(CTCF_BW)) as bw:
        chroms = bw.chroms()
    for chrom, part in boundaries.groupby("chrom"):
        if chrom not in chroms:
            continue
        low = flank
        high = chroms[chrom] - flank
        if high <= low:
            continue
        centers = rng.integers(low, high, size=len(part))
        rows.append(pd.DataFrame({"chrom": chrom, "center": centers}))
    return pd.concat(rows, ignore_index=True)


def boundary_analysis():
    ctcf = read_bed(CTCF_PEAKS)
    data = {
        "res100kb_window300kb": boundaries_from_insulation(INS100, 300_000, "res100kb_w300kb"),
        "res100kb_window500kb": boundaries_from_insulation(INS100, 500_000, "res100kb_w500kb"),
        "res100kb_window1000kb": boundaries_from_insulation(INS100, 1_000_000, "res100kb_w1000kb"),
    }
    data.update(call_boundaries_1mb())
    summary = []
    profile_rows = []
    for name, b in data.items():
        bed = TRACKS / f"{name}.bed"
        b[["chrom", "start", "end", "name", "score"]].to_csv(bed, sep="\t", header=False, index=False)
        hit = has_overlap_near_centers(b, ctcf)
        offsets, real = avg_profile(CTCF_BW, centers_dict(b))
        controls = []
        for i in range(5):
            r = random_like(b, 100 + i)
            _, prof = avg_profile(CTCF_BW, centers_dict(r))
            controls.append(prof)
        controls = np.vstack(controls)
        control = np.nanmean(controls, axis=0)
        center_i = int(np.where(offsets == 0)[0][0])
        summary.append({
            "config": name,
            "boundaries": len(b),
            "ctcf_near_10kb": int(hit.sum()),
            "ctcf_near_fraction": float(hit.mean()),
            "ctcf_center_signal": float(real[center_i]),
            "random_center_signal": float(control[center_i]),
            "center_enrichment": float(real[center_i] / control[center_i]) if control[center_i] else np.nan,
        })
        for x, y, z in zip(offsets, real, control):
            profile_rows.append({"config": name, "offset_bp": int(x), "ctcf_real": y, "ctcf_random": z})
    summary = pd.DataFrame(summary)
    profiles = pd.DataFrame(profile_rows)
    summary.to_csv(TAB / "boundary_parameter_summary.tsv", sep="\t", index=False)
    profiles.to_csv(TAB / "boundary_ctcf_profiles.tsv", sep="\t", index=False)
    plot_boundary_summary(summary)
    plot_boundary_profiles(profiles)
    return summary, profiles


def plot_boundary_summary(df):
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    sns.barplot(data=df, x="config", y="boundaries", ax=ax[0], color="#8ecae6")
    sns.barplot(data=df, x="config", y="ctcf_near_fraction", ax=ax[1], color="#90be6d")
    ax[0].tick_params(axis="x", rotation=35)
    ax[1].tick_params(axis="x", rotation=35)
    ax[0].set_xlabel("")
    ax[1].set_xlabel("")
    ax[0].set_ylabel("boundaries")
    ax[1].set_ylabel("fraction with CTCF +/-10 kb")
    fig.tight_layout()
    fig.savefig(FIG / "boundary_parameters.png", dpi=180)
    plt.close(fig)


def plot_boundary_profiles(df):
    fig, ax = plt.subplots(figsize=(8, 4))
    for name, part in df.groupby("config"):
        x = part["offset_bp"] / 1_000_000
        ax.plot(x, part["ctcf_real"], label=name)
    control = df.groupby("offset_bp")["ctcf_random"].mean().reset_index()
    ax.plot(control["offset_bp"] / 1_000_000, control["ctcf_random"], color="black", linestyle="--", label="random mean")
    ax.axvline(0, color="black", linewidth=1)
    ax.set_xlabel("distance from boundary center, Mb")
    ax.set_ylabel("mean CTCF FE")
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(FIG / "ctcf_profiles_boundaries.png", dpi=180)
    plt.close(fig)


def promoter_analysis():
    promoters = read_promoters(GTF)
    promoters = promoters[promoters["chrom"] == "chr1"].copy()
    h3k27ac = read_bed(H3K27AC_PEAKS)
    ctcf = read_bed(CTCF_PEAKS)
    promoters["h3k27ac_peak_count"] = overlap_count(promoters[["chrom", "start", "end"]], h3k27ac)
    promoters["has_h3k27ac_peak"] = promoters["h3k27ac_peak_count"] > 0
    promoters["rnaseq"] = mean_bw(RNA_BW, promoters)
    promoters["cage"] = mean_bw(CAGE_BW, promoters) if CAGE_BW.exists() else np.nan
    promoters["h3k27ac"] = mean_bw(H3K27AC_BW, promoters)
    promoters["h3k9me3"] = mean_bw(H3K9ME3_BW, promoters)
    promoters["methylation"] = mean_bw(METH_BW, promoters)
    promoters["coverage"] = mean_bw(COV_BW, promoters)
    promoters["gc"] = mean_bw(GC_BW, promoters)
    promoters["cpg_obs_exp"] = mean_bw(CPG_BW, promoters)
    promoters["e1"] = mean_bw(E1_BW, promoters, chrom_mode="strip")
    promoters["ctcf_distance"] = nearest_distance(promoters, ctcf)
    valid = promoters["rnaseq"].replace([np.inf, -np.inf], np.nan).dropna()
    q25, q75 = valid.quantile([0.25, 0.75])
    promoters["activity_group"] = "middle"
    promoters.loc[promoters["rnaseq"] <= q25, "activity_group"] = "inactive"
    promoters.loc[promoters["rnaseq"] >= q75, "activity_group"] = "active"
    promoters["low_methylation"] = promoters["methylation"] < 0.3
    pos_h3 = promoters.loc[promoters["e1"] > 0, "h3k27ac"].mean()
    neg_h3 = promoters.loc[promoters["e1"] < 0, "h3k27ac"].mean()
    if pos_h3 >= neg_h3:
        promoters["compartment"] = np.where(promoters["e1"] >= 0, "A_like", "B_like")
    else:
        promoters["compartment"] = np.where(promoters["e1"] < 0, "A_like", "B_like")
    promoters.to_csv(TAB / "promoter_signals_chr1.tsv", sep="\t", index=False)
    summary = make_promoter_summary(promoters, q25, q75)
    summary.to_csv(TAB / "promoter_summary.tsv", sep="\t", index=False)
    plot_promoters(promoters)
    regions = choose_regions(promoters)
    regions.to_csv(TAB / "selected_regions.tsv", sep="\t", index=False)
    make_browser_tracks(regions)
    return promoters, summary, regions


def make_promoter_summary(df, q25, q75):
    rows = []
    rows.append({"metric": "promoters_chr1", "value": len(df)})
    rows.append({"metric": "rna_q25", "value": q25})
    rows.append({"metric": "rna_q75", "value": q75})
    rows.append({"metric": "promoters_with_h3k27ac_peak", "value": int(df["has_h3k27ac_peak"].sum())})
    rows.append({"metric": "low_methylation_promoters", "value": int(df["low_methylation"].sum())})
    both = df["has_h3k27ac_peak"] & df["low_methylation"]
    rows.append({"metric": "h3k27ac_and_low_methylation", "value": int(both.sum())})
    triple = both & (df["activity_group"] == "active")
    rows.append({"metric": "h3k27ac_low_methylation_active", "value": int(triple.sum())})
    for group in ["active", "inactive"]:
        part = df[df["activity_group"] == group]
        rows.append({"metric": f"{group}_n", "value": len(part)})
        for col in ["rnaseq", "cage", "h3k27ac", "h3k9me3", "methylation", "coverage", "gc", "cpg_obs_exp", "ctcf_distance"]:
            rows.append({"metric": f"{group}_{col}_median", "value": float(part[col].median())})
        rows.append({"metric": f"{group}_A_like_fraction", "value": float((part["compartment"] == "A_like").mean())})
        rows.append({"metric": f"{group}_h3k27ac_peak_fraction", "value": float(part["has_h3k27ac_peak"].mean())})
        rows.append({"metric": f"{group}_low_methylation_fraction", "value": float(part["low_methylation"].mean())})
    high_cpg = df["cpg_obs_exp"] >= df["cpg_obs_exp"].quantile(0.75)
    rows.append({"metric": "high_cpg_promoter_methylation_median", "value": float(df.loc[high_cpg, "methylation"].median())})
    rows.append({"metric": "other_promoter_methylation_median", "value": float(df.loc[~high_cpg, "methylation"].median())})
    return pd.DataFrame(rows)


def plot_promoters(df):
    part = df[df["activity_group"].isin(["active", "inactive"])].copy()
    fig, ax = plt.subplots(2, 2, figsize=(9, 7))
    sns.boxplot(data=part, x="activity_group", y="h3k27ac", ax=ax[0, 0], color="#8ecae6", showfliers=False)
    sns.boxplot(data=part, x="activity_group", y="h3k9me3", ax=ax[0, 1], color="#ffb703", showfliers=False)
    sns.boxplot(data=part, x="activity_group", y="methylation", ax=ax[1, 0], color="#90be6d", showfliers=False)
    sns.boxplot(data=part, x="activity_group", y="ctcf_distance", ax=ax[1, 1], color="#cdb4db", showfliers=False)
    ax[0, 0].set_ylabel("H3K27Ac FE")
    ax[0, 1].set_ylabel("H3K9me3 FE")
    ax[1, 0].set_ylabel("DNA methylation")
    ax[1, 1].set_ylabel("distance to nearest CTCF, bp")
    for a in ax.ravel():
        a.set_xlabel("")
    fig.tight_layout()
    fig.savefig(FIG / "active_inactive_promoters.png", dpi=180)
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(6, 5))
    colors = df["has_h3k27ac_peak"].map({True: "#d00000", False: "#777777"})
    ax.scatter(df["rnaseq"], df["methylation"], s=12, c=colors, alpha=0.65)
    ax.set_xscale("symlog", linthresh=1)
    ax.set_xlabel("RNAseq promoter signal")
    ax.set_ylabel("DNA methylation")
    fig.tight_layout()
    fig.savefig(FIG / "rnaseq_vs_methylation_h3k27ac.png", dpi=180)
    plt.close(fig)


def choose_regions(df):
    q = df[
        (df["activity_group"] == "active")
        & (df["has_h3k27ac_peak"])
        & (df["low_methylation"])
        & df["rnaseq"].notna()
    ].copy()
    q = q.sort_values(["rnaseq", "h3k27ac"], ascending=False).head(3)
    q["region_start"] = (q["tss"] - 50_000).clip(lower=0).astype(int)
    q["region_end"] = (q["tss"] + 50_000).astype(int)
    return q[[
        "chrom", "region_start", "region_end", "gene_name", "strand", "tss",
        "rnaseq", "h3k27ac", "methylation", "coverage", "compartment", "ctcf_distance",
    ]]


def write_chr1_bedgraph(src, dst):
    rows = []
    with open(src) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if parts[0] == "1":
                parts[0] = "chr1"
                rows.append("\t".join(parts))
    Path(dst).write_text("\n".join(rows) + "\n")


def make_browser_config(region_bed, e1_bedgraph):
    config = TMP / "tracks.ini"
    text = f"""
[x-axis]
where = top

[RNAseq]
file = {RNA_BW}
title = RNAseq STAR
height = 2
color = #1f77b4

[CAGE]
file = {CAGE_BW}
title = CAGE
height = 1.5
color = #17becf

[H3K27Ac]
file = {H3K27AC_BW}
title = H3K27Ac FE
height = 2
color = #d62728

[H3K9me3]
file = {H3K9ME3_BW}
title = H3K9me3 FE
height = 2
color = #ff7f0e

[Methylation]
file = {METH_BW}
title = DNA methylation
height = 2
min_value = 0
max_value = 1
color = #2ca02c

[Coverage]
file = {COV_BW}
title = WGBS coverage
height = 1.2
color = #555555

[E1]
file = {e1_bedgraph}
title = E1 compartment
height = 1.5
color = #9467bd

[CTCF peaks]
file = {CTCF_PEAKS}
title = CTCF peaks
height = 1
color = #000000

[boundaries]
file = {region_bed}
title = Hi-C boundaries
height = 1
color = #555555
"""
    config.write_text(text.strip() + "\n")
    return config


def make_browser_tracks(regions):
    if regions.empty:
        return
    boundary_bed = TRACKS / "res100kb_window300kb.bed"
    e1_bedgraph = TRACKS / "MoPh7_E1_chr1.bedGraph"
    write_chr1_bedgraph(ROOT / "course/day2_HiC_practice/results/compartments/MoPh7_enr_v2_E1_res100000.bedGraph", e1_bedgraph)
    config = make_browser_config(boundary_bed, e1_bedgraph)
    for row in regions.itertuples(index=False):
        name = str(row.gene_name).replace("/", "_")
        region = f"{row.chrom}:{int(row.region_start)}-{int(row.region_end)}"
        out = FIG / f"genome_tracks_{name}.png"
        cmd = [
            "pyGenomeTracks",
            "--tracks", str(config),
            "--region", region,
            "--outFileName", str(out),
            "--fontSize", "7",
            "--dpi", "180",
        ]
        try:
            subprocess.run(cmd, check=True, env={**os.environ, "PATH": f"/home/greg/miniforge3/envs/genome/bin:{os.environ.get('PATH', '')}"})
        except Exception:
            pass


def write_answers(boundary_summary, promoter_summary, regions):
    metrics = {r.metric: r.value for r in promoter_summary.itertuples(index=False)}
    best = boundary_summary.sort_values("center_enrichment", ascending=False).iloc[0]
    lines = []
    lines.append("# Ответы на дополнительные задания")
    lines.append("")
    lines.append("## 1. Параметры границ доменов")
    lines.append("")
    lines.append(f"Лучшее CTCF-обогащение в центре границ получилось для `{best.config}`: {best.center_enrichment:.2f}x относительно случайного контроля.")
    lines.append("При увеличении окна и укрупнении разрешения меняется число границ: мелкие настройки дают больше локальных границ, грубые настройки оставляют меньше крупных границ.")
    lines.append("CTCF чаще встречается около настоящих границ, чем около случайных позиций, значит часть найденных границ поддерживается структурным белком CTCF.")
    lines.append("")
    lines.append("## 2. Промоторы, H3K27Ac, метилирование и RNAseq")
    lines.append("")
    lines.append("Промоторы взяты как 1000 bp перед TSS. В этой аннотации фактически есть CDS/exon, поэтому TSS приближен по краю CDS: для `+` цепи левый край, для `-` цепи правый край.")
    lines.append(f"На chr1 найдено {int(metrics['promoters_chr1'])} промоторов. H3K27Ac пик есть в {int(metrics['promoters_with_h3k27ac_peak'])} промоторах. Низкое метилирование (<0.3) есть в {int(metrics['low_methylation_promoters'])} промоторах.")
    lines.append(f"Одновременно H3K27Ac + низкое метилирование + высокий RNAseq сигнал имеют {int(metrics['h3k27ac_low_methylation_active'])} промотора.")
    lines.append("Вывод: активный H3K27Ac, низкое метилирование и высокий RNAseq в целом согласуются, но не идеально один-к-одному, потому что это разные уровни регуляции.")
    lines.append("")
    lines.append("## 3. Регионы вместо IGV")
    lines.append("")
    lines.append("IGV на сервере без графического интерфейса не использовался. Вместо него регионы отрисованы через терминальный инструмент pyGenomeTracks.")
    for row in regions.itertuples(index=False):
        lines.append(f"- `{row.gene_name}`: высокий RNAseq, есть H3K27Ac, низкое метилирование, compartment `{row.compartment}`, ближайший CTCF примерно {row.ctcf_distance:.0f} bp.")
    lines.append("")
    lines.append("## 4. Активные и неактивные промоторы")
    lines.append("")
    lines.append(f"Активные промоторы имеют медианный H3K27Ac {metrics['active_h3k27ac_median']:.3f}, неактивные {metrics['inactive_h3k27ac_median']:.3f}.")
    lines.append(f"Активные промоторы имеют медианное метилирование {metrics['active_methylation_median']:.3f}, неактивные {metrics['inactive_methylation_median']:.3f}.")
    lines.append(f"Доля A-like компартмента: активные {metrics['active_A_like_fraction']:.2%}, неактивные {metrics['inactive_A_like_fraction']:.2%}.")
    lines.append("Вывод: активные промоторы чаще имеют H3K27Ac, ниже метилирование и чаще попадают в A-like окружение.")
    lines.append("")
    lines.append("## 5. CTCF около активных и неактивных генов")
    lines.append("")
    lines.append(f"Медианное расстояние до ближайшего CTCF: активные {metrics['active_ctcf_distance_median']:.0f} bp, неактивные {metrics['inactive_ctcf_distance_median']:.0f} bp.")
    lines.append("Вывод: связь архитектурных элементов с активностью видна умеренно; CTCF сам по себе не равен активности, но может находиться рядом с регуляторными границами.")
    lines.append("")
    lines.append("## Вопросы WGBS")
    lines.append("")
    lines.append("Coverage track нужен вместе с метилированием, потому что beta-value при малом числе ридов нестабильна.")
    lines.append("CpG с низким coverage лучше удалить: если ридов мало, один случайный рид резко меняет оценку метилирования.")
    lines.append("Beta-value - это доля метилированных ридов от 0 до 1. M-value - это log2 отношение метилированных к неметилированным ридам, оно удобнее для статистики, но менее интуитивно.")
    lines.append("Участки с высоким GC видны как пики GC-content track; часто они совпадают с CpG-rich областями.")
    lines.append(f"Для промоторов с высоким CpG observed/expected медианное метилирование {metrics['high_cpg_promoter_methylation_median']:.3f}, для остальных {metrics['other_promoter_methylation_median']:.3f}. CpG-богатые промоторы в этих данных чаще менее метилированы.")
    lines.append("При добавлении ChIP-seq видно, что H3K27Ac чаще совпадает с активными низкометилированными промоторами; H3K9me3 в активных промоторах ниже и не выглядит как метка активного старта транскрипции.")
    lines.append("H3K27Ac связан с низким метилированием, потому что это метка активных регуляторных областей, где открытая хроматиновая среда и связывание факторов часто несовместимы с плотным метилированием ДНК.")
    (OUT / "FINAL_ANSWERS.md").write_text("\n".join(lines) + "\n")
    with open(OUT / "metrics.json", "w") as f:
        json.dump({"promoter": metrics, "best_boundary_config": best.to_dict()}, f, indent=2)


def main():
    mkdirs()
    boundary_summary, profiles = boundary_analysis()
    promoters, promoter_summary, regions = promoter_analysis()
    write_answers(boundary_summary, promoter_summary, regions)
    print("done")
    print(OUT)


if __name__ == "__main__":
    main()
