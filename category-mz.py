# -*- coding: utf-8 -*-
# 分类规则：按 H/C 与 O/C 将分子分为 6 类（未命中归为 Other）
# 统计：按 m/z 分箱（默认 25 Da），计算每箱各类相对丰度（%）
# 输出：两个 CSV（counts 与 percent）+ 一张堆叠面积图 PNG

import os, math
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# ======== 配置 ========
INPUT_PATH = r"D:\data\category-mz\NOM-14.txt"  # ← 改成你的文件
SAMPLE_TAG = "NOM-14"
BIN_WIDTH_DA = 25                               # 分箱宽度（Da）
# =====================

in_path = Path(INPUT_PATH)
out_dir = in_path.parent / f"{SAMPLE_TAG}_mw_class_outputs"
out_dir.mkdir(parents=True, exist_ok=True)
base = out_dir / f"{SAMPLE_TAG}_mw_class"

# 1) 读取数据（假设是制表符分隔）
df = pd.read_csv(INPUT_PATH, sep="\t")

# 2) 列名标准化
def norm_name(c):
    c2 = c.strip().lower().replace(" ", "").replace("/", "_")
    return c2
col_map = {c: norm_name(c) for c in df.columns}
df = df.rename(columns=col_map)

# 尝试识别常见写法
name_map = {}
for c in df.columns:
    lc = c
    if lc in ["m_z", "mz", "mass", "masses"]:
        name_map[c] = "mz"
    elif lc in ["o_c", "oc"]:
        name_map[c] = "o_c"
    elif lc in ["h_c", "hc"]:
        name_map[c] = "h_c"

df = df.rename(columns=name_map)
required = {"mz", "o_c", "h_c"}
missing = required - set(df.columns)
if missing:
    raise ValueError(f"缺少必要列: {missing}；现有列: {list(df.columns)}")

# 3) 数值清洗
df["mz"]  = pd.to_numeric(df["mz"],  errors="coerce")
df["o_c"] = pd.to_numeric(df["o_c"], errors="coerce")
df["h_c"] = pd.to_numeric(df["h_c"], errors="coerce")
df = df.dropna(subset=["mz", "o_c", "h_c"]).copy()

# 4) 分类函数（与你给的阈值完全一致）
def classify_row(hc, oc):
    if (hc >= 1.5 and hc <= 2.0) and (oc >= 0 and oc < 0.3):
        return "Lipids"
    elif (hc >= 1.5 and hc <= 2.2) and (oc >= 0.3 and oc < 0.67):
        return "Proteins"
    elif (hc >= 1.5 and hc <= 2.2) and (oc >= 0.67 and oc <= 1.2):
        return "Carbohydrates"
    elif (hc >= 0.7 and hc < 1.5) and (oc >= 0.1 and oc <= 0.67):
        return "Lignins"
    elif (hc >= 0.2 and hc < 0.7) and (oc >= 0 and oc <= 0.67):
        return "Condensed Aromatics"
    elif (hc >= 0.5 and hc < 1.5) and (oc >= 0.67 and oc <= 1.2):
        return "Tannins"
    else:
        return "Other"

df["Class"] = [classify_row(h, o) for h, o in zip(df["h_c"], df["o_c"])]

# 5) m/z 分箱
min_mz = math.floor(df["mz"].min() / BIN_WIDTH_DA) * BIN_WIDTH_DA
max_mz = math.ceil(df["mz"].max() / BIN_WIDTH_DA) * BIN_WIDTH_DA
bins = np.arange(min_mz, max_mz + BIN_WIDTH_DA, BIN_WIDTH_DA)
labels = [f"{int(b)}-{int(b+BIN_WIDTH_DA)}" for b in bins[:-1]]  # 用连字符避免字体问题
df["MW_bin"] = pd.cut(df["mz"], bins=bins, labels=labels, include_lowest=True, right=False)

# 6) 聚合计数、补全空格子
keep_classes = ["Lipids","Proteins","Carbohydrates","Lignins","Condensed Aromatics","Tannins","Other"]
counts = df.groupby(["MW_bin","Class"]).size().reset_index(name="Count")
grid = pd.MultiIndex.from_product([pd.Index(labels, name="MW_bin"), keep_classes],
                                  names=["MW_bin","Class"])
counts_full = counts.set_index(["MW_bin","Class"]).reindex(grid, fill_value=0).reset_index()

# 7) 透视为宽表
wide_counts = counts_full.pivot(index="MW_bin", columns="Class", values="Count").fillna(0)
wide_counts = wide_counts[keep_classes]  # 保持列顺序
wide_counts.insert(0, "MW_bin", wide_counts.index.astype(str))

# 8) 百分比表
wide_percent = wide_counts.copy()
tot = wide_percent[keep_classes].sum(axis=1).replace(0, np.nan)
for c in keep_classes:
    wide_percent[c] = (wide_percent[c] / tot) * 100.0
wide_percent = wide_percent.fillna(0)

# 9) 导出
counts_csv = f"{base}_counts.csv"
percent_csv = f"{base}_percent.csv"
wide_counts.to_csv(counts_csv, index=False, encoding="utf-8-sig")
wide_percent.to_csv(percent_csv, index=False, encoding="utf-8-sig")

# 10) 作图（堆叠面积）
plt.figure(figsize=(10, 6))
x = np.arange(len(wide_percent))
y_stack = [wide_percent[c].values for c in keep_classes]
plt.stackplot(x, y_stack, labels=keep_classes)
plt.xticks(x, wide_percent["MW_bin"], rotation=60, ha="right")
plt.ylabel("Relative Abundance (%)")
plt.xlabel("Molecular weight (Da)")
plt.title(f"{SAMPLE_TAG}: Class composition across molecular-weight bins")
plt.legend(loc="upper right", ncol=2)
plt.tight_layout()
png_path = f"{base}_stacked_area.png"
plt.savefig(png_path, dpi=200)
plt.show()

print("已导出：")
print(" -", counts_csv)
print(" -", percent_csv)
print(" -", png_path)
