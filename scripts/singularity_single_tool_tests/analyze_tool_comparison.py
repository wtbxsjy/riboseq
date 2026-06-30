#!/usr/bin/env python3
"""
工具方法学比较分析
分析不同工具对 sORF 检测的一致性和特异性
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict

# ========== 配置 ==========
INPUT_FILE = "results_all_tools_integrated/Mouse_AllTools_Comparison.orfs.out"
OUTPUT_PREFIX = "tool_comparison"

# 样本列表
SAMPLES = [f"Sample{i:02d}" for i in range(1, 37)]  # Sample01-Sample36

# 工具列表
TOOLS = ["ORFquant", "RiboTISH", "Ribotricer"]

print("=" * 60)
print("工具方法学比较分析")
print("=" * 60)

# ========== 1. 加载数据 ==========
print("\n[1] 加载数据...")
df = pd.read_csv(INPUT_FILE, sep='\t')
print(f"   总 ORF 数: {len(df)}")

# ========== 2. 提取工具检测信息 ==========
print("\n[2] 提取工具检测矩阵...")

# 为每个工具创建检测矩阵
tool_detection = {}

for tool in TOOLS:
    # 查找包含工具名的列
    tool_cols = [col for col in df.columns if tool in col]

    if len(tool_cols) == 0:
        print(f"   警告：未找到 {tool} 的列")
        continue

    # 每个 ORF 在该工具的检测情况（任一样本检测到即为1）
    tool_detection[tool] = df[tool_cols].max(axis=1)

    detected_count = tool_detection[tool].sum()
    print(f"   {tool}: {detected_count} ORFs ({detected_count/len(df)*100:.1f}%)")

# 转换为 DataFrame
detection_df = pd.DataFrame(tool_detection)
detection_df.index = df['orf_id']

# ========== 3. 工具交集分析（维恩图统计）==========
print("\n[3] 工具交集分析...")

def get_tool_combinations():
    """计算所有工具组合的 ORF 数量"""
    results = {}

    # 单个工具独有
    for tool in TOOLS:
        mask = (detection_df[tool] == 1)
        for other_tool in TOOLS:
            if other_tool != tool:
                mask &= (detection_df[other_tool] == 0)
        results[f"{tool}_only"] = mask.sum()

    # 两两交集
    from itertools import combinations
    for tool1, tool2 in combinations(TOOLS, 2):
        mask = (detection_df[tool1] == 1) & (detection_df[tool2] == 1)
        # 排除第三个工具
        other_tool = [t for t in TOOLS if t not in [tool1, tool2]][0]
        mask &= (detection_df[other_tool] == 0)
        results[f"{tool1}_{tool2}_only"] = mask.sum()

    # 三工具共同
    mask = (detection_df[TOOLS[0]] == 1) & \
           (detection_df[TOOLS[1]] == 1) & \
           (detection_df[TOOLS[2]] == 1)
    results["all_three"] = mask.sum()

    return results

combinations = get_tool_combinations()

print("\n   工具组合统计：")
for combo, count in combinations.items():
    print(f"   {combo:30s}: {count:6d} ORFs")

# ========== 4. ORF 类型分布（按工具）==========
print("\n[4] ORF 类型分布分析...")

orf_type_stats = []

for tool in TOOLS:
    if tool not in tool_detection:
        continue

    # 该工具检测到的 ORF
    tool_orfs = df[tool_detection[tool] == 1]

    # 按类型统计
    type_counts = tool_orfs['orf_biotype'].value_counts()

    for orf_type, count in type_counts.items():
        orf_type_stats.append({
            'Tool': tool,
            'ORF_Type': orf_type,
            'Count': count
        })

orf_type_df = pd.DataFrame(orf_type_stats)
print(orf_type_df.pivot(index='ORF_Type', columns='Tool', values='Count').fillna(0))

# ========== 5. 样本重复性分析 ==========
print("\n[5] 样本重复性分析...")

for tool in TOOLS:
    tool_cols = [col for col in df.columns if tool in col]
    if len(tool_cols) == 0:
        continue

    # 计算每个 ORF 在多少个样本中被检测到
    sample_support = df[tool_cols].sum(axis=1)

    print(f"\n   {tool} 样本支持度分布：")
    print(f"   1-5 样本:   {(sample_support <= 5).sum():6d} ORFs ({(sample_support <= 5).sum()/len(df)*100:.1f}%)")
    print(f"   6-10 样本:  {((sample_support > 5) & (sample_support <= 10)).sum():6d} ORFs")
    print(f"   11-20 样本: {((sample_support > 10) & (sample_support <= 20)).sum():6d} ORFs")
    print(f"   21-30 样本: {((sample_support > 20) & (sample_support <= 30)).sum():6d} ORFs")
    print(f"   >30 样本:   {(sample_support > 30).sum():6d} ORFs ({(sample_support > 30).sum()/len(df)*100:.1f}%)")

# ========== 6. 高置信度 ORF 筛选 ==========
print("\n[6] 高置信度 ORF 筛选...")

# 标准 1：被所有三个工具检测到
all_tools = (detection_df[TOOLS[0]] == 1) & \
            (detection_df[TOOLS[1]] == 1) & \
            (detection_df[TOOLS[2]] == 1)

high_conf_all_tools = df[all_tools]
print(f"\n   标准 1（所有工具）: {len(high_conf_all_tools)} ORFs")

# 标准 2：被至少2个工具检测到
at_least_two = detection_df[TOOLS].sum(axis=1) >= 2
high_conf_two_tools = df[at_least_two]
print(f"   标准 2（≥2工具）:   {len(high_conf_two_tools)} ORFs")

# 标准 3：被至少2个工具 + 至少10个样本检测到
tool_cols_all = []
for tool in TOOLS:
    tool_cols_all.extend([col for col in df.columns if tool in col])

sample_support_all = df[tool_cols_all].sum(axis=1)
high_conf_stringent = at_least_two & (sample_support_all >= 10)
high_conf_stringent_orfs = df[high_conf_stringent]
print(f"   标准 3（≥2工具+≥10样本）: {len(high_conf_stringent_orfs)} ORFs")

# 保存高置信度 ORF
high_conf_stringent_orfs.to_csv(f"{OUTPUT_PREFIX}_high_confidence.tsv", sep='\t', index=False)
print(f"\n   ✅ 高置信度 ORF 已保存: {OUTPUT_PREFIX}_high_confidence.tsv")

# ========== 7. 工具特异性 ORF ==========
print("\n[7] 工具特异性 ORF 分析...")

for tool in TOOLS:
    if tool not in tool_detection:
        continue

    # 仅被该工具检测到的 ORF
    only_this = tool_detection[tool] == 1
    for other_tool in TOOLS:
        if other_tool != tool and other_tool in tool_detection:
            only_this &= (tool_detection[other_tool] == 0)

    specific_orfs = df[only_this]

    print(f"\n   {tool} 特异性 ORF: {len(specific_orfs)} 个")

    # 分析类型分布
    if len(specific_orfs) > 0:
        type_dist = specific_orfs['orf_biotype'].value_counts().head(5)
        print(f"   主要类型:")
        for orf_type, count in type_dist.items():
            print(f"     - {orf_type}: {count}")

    # 保存特异性 ORF
    specific_orfs.to_csv(f"{OUTPUT_PREFIX}_{tool}_specific.tsv", sep='\t', index=False)

# ========== 8. 生成总结报告 ==========
print("\n[8] 生成总结报告...")

with open(f"{OUTPUT_PREFIX}_summary_report.txt", 'w') as f:
    f.write("=" * 80 + "\n")
    f.write("工具方法学比较分析报告\n")
    f.write("=" * 80 + "\n\n")

    f.write(f"总样本数: {len(SAMPLES)}\n")
    f.write(f"工具数: {len(TOOLS)}\n")
    f.write(f"工具: {', '.join(TOOLS)}\n\n")

    f.write("=" * 80 + "\n")
    f.write("1. 总体统计\n")
    f.write("=" * 80 + "\n")
    f.write(f"唯一 ORF 数（合并后）: {len(df)}\n\n")

    for tool in TOOLS:
        if tool in tool_detection:
            count = tool_detection[tool].sum()
            f.write(f"{tool:20s}: {count:6d} ORFs ({count/len(df)*100:.1f}%)\n")

    f.write("\n" + "=" * 80 + "\n")
    f.write("2. 工具交集分析\n")
    f.write("=" * 80 + "\n")
    for combo, count in combinations.items():
        f.write(f"{combo:30s}: {count:6d} ORFs\n")

    f.write("\n" + "=" * 80 + "\n")
    f.write("3. 高置信度 ORF（推荐用于下游分析）\n")
    f.write("=" * 80 + "\n")
    f.write(f"标准 1（所有工具）:          {len(high_conf_all_tools):6d} ORFs\n")
    f.write(f"标准 2（≥2工具）:            {len(high_conf_two_tools):6d} ORFs\n")
    f.write(f"标准 3（≥2工具+≥10样本）:   {len(high_conf_stringent_orfs):6d} ORFs\n")

print(f"\n✅ 总结报告已保存: {OUTPUT_PREFIX}_summary_report.txt")

print("\n" + "=" * 60)
print("分析完成！")
print("=" * 60)
