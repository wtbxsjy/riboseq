import os
import glob
import csv
import json
import argparse
import re
import sys

def parse_args():
    parser = argparse.ArgumentParser(description="根据 fastq 文件生成 Sample Sheet CSV")
    
    parser.add_argument('-i', '--input_dir', required=True, help="Fastq/fq 文件所在的目录路径")
    parser.add_argument('-o', '--output_csv', required=True, help="输出 CSV 文件的路径 (例如: samples.csv)")
    parser.add_argument('--strandedness', default='auto', help="CSV 中的 strandedness 列的值 (默认: auto)")
    parser.add_argument('--type', default='riboseq', help="CSV 中的 type 列的值 (默认: riboseq)")
    parser.add_argument(
        '--group-map',
        default=None,
        help=(
            "可选：指定样本到分组的映射，用于生成 samplesheet 中的 group 列。\n"
            "接受两种格式：\n"
            "  1. JSON 文件：{\"sample_id\": \"group_name\", ...}\n"
            "  2. CSV 文件（两列，无表头）：sample_id,group_name\n"
            "未出现在映射中的样本不写 group 列（留空）。\n"
            "示例 JSON: {\"rep1\": \"treatment_A\", \"rep2\": \"treatment_A\"}"
        )
    )
    
    return parser.parse_args()


def load_group_map(group_map_path):
    """Load sample-to-group mapping from a JSON or two-column CSV file."""
    if not os.path.exists(group_map_path):
        print(f"错误: group-map 文件不存在: {group_map_path}", file=sys.stderr)
        sys.exit(1)

    ext = os.path.splitext(group_map_path)[1].lower()
    mapping = {}

    if ext == '.json':
        with open(group_map_path) as f:
            mapping = json.load(f)
        if not isinstance(mapping, dict):
            print("错误: JSON group-map 必须是 {sample_id: group_name} 格式", file=sys.stderr)
            sys.exit(1)
    else:
        # Treat as two-column CSV (no header): sample_id,group_name
        with open(group_map_path, newline='') as f:
            reader = csv.reader(f)
            for row in reader:
                if len(row) >= 2:
                    mapping[row[0].strip()] = row[1].strip()

    return mapping

def get_sample_name(filename, mode):
    """
    根据文件名和模式(PE/SE)提取样品名
    去除 _R1, _1, .fastq.gz, .fq.gz 等后缀
    """
    # 先去除扩展名
    base = re.sub(r'\.(fastq|fq)\.gz$', '', filename)
    
    if mode == 'PE':
        # 如果是双端，去除结尾的 _R1, _1 等标记
        # 匹配 _R1, _1, _R1_001 等，且位于字符串末尾
        base = re.sub(r'(_R1|_1)(_001)?$', '', base)
    
    return base

def main():
    args = parse_args()
    
    input_dir = os.path.abspath(args.input_dir)
    if not os.path.exists(input_dir):
        print(f"错误: 目录 {input_dir} 不存在")
        sys.exit(1)

    # Load optional group mapping
    group_map = {}
    if args.group_map:
        group_map = load_group_map(args.group_map)
        print(f"已加载 group 映射，共 {len(group_map)} 个样本")

    # 寻找所有的 .fastq.gz 和 .fq.gz 文件
    patterns = ['*.fastq.gz', '*.fq.gz']
    files = []
    for p in patterns:
        files.extend(glob.glob(os.path.join(input_dir, p)))
    
    files.sort() # 排序以确保稳定性
    
    samples = {} # 字典用于临时存储配对信息
    processed_files = set()

    print(f"扫描目录: {input_dir}")
    print(f"找到 {len(files)} 个 fastq 文件")

    # 准备 CSV 数据列表
    csv_rows = []

    for f_path in files:
        if f_path in processed_files:
            continue
            
        filename = os.path.basename(f_path)
        dir_path = os.path.dirname(f_path)
        
        # 逻辑：判断是否为 Read 1
        # 常见模式: _R1.fastq.gz, _R1_001.fastq.gz, _1.fq.gz
        is_r1 = False
        r2_filename = None
        
        # 尝试推测 R2 文件名
        if '_R1' in filename:
            potential_r2 = filename.replace('_R1', '_R2')
            is_r1 = True
        elif '_1.' in filename: # 匹配 _1.fq.gz 或 _1.fastq.gz
            # 小心不要匹配到 sample_10.fq.gz 这种，所以匹配 _1.
            potential_r2 = filename.replace('_1.', '_2.')
            is_r1 = True
        else:
            potential_r2 = None
        
        # 检查 R2 是否存在
        r2_path = os.path.join(dir_path, potential_r2) if potential_r2 else None
        
        if is_r1 and r2_path and os.path.exists(r2_path):
            # === Paired End (PE) ===
            mode = 'PE'
            fastq_1 = f_path
            fastq_2 = r2_path
            
            # 标记 R2 已处理，避免下次循环重复
            processed_files.add(r2_path)
            
        else:
            # === Single End (SE) ===
            # 或者是 Read 1 文件找不到对应的 Read 2
            # 或者是文件名里根本没有 _R1/_1
            mode = 'SE'
            fastq_1 = f_path
            fastq_2 = "" # CSV留空
        
        # 提取样品名
        sample_name = get_sample_name(filename, mode)
        
        # Resolve group for this sample (empty string if not mapped)
        group = group_map.get(sample_name, '')

        # 添加到结果列表
        csv_rows.append({
            'sample': sample_name,
            'fastq_1': fastq_1,
            'fastq_2': fastq_2,
            'strandedness': args.strandedness,
            'type': args.type,
            'group': group,
        })
        
        processed_files.add(f_path)

    # Write CSV — include 'group' column only if any sample has a group assigned
    has_groups = any(row.get('group') for row in csv_rows)
    if has_groups:
        headers = ['sample', 'fastq_1', 'fastq_2', 'strandedness', 'type', 'group']
    else:
        headers = ['sample', 'fastq_1', 'fastq_2', 'strandedness', 'type']
    
    try:
        with open(args.output_csv, 'w', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=headers, extrasaction='ignore')
            writer.writeheader()
            for row in csv_rows:
                writer.writerow(row)
        
        print(f"------------------------------------------------")
        print(f"成功生成 CSV 文件: {args.output_csv}")
        print(f"共包含 {len(csv_rows)} 个样品")
        
    except IOError as e:
        print(f"写入 CSV 失败: {e}")

if __name__ == "__main__":
    main()
