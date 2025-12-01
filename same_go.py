import pandas as pd
import argparse
from pathlib import Path
import shutil
from datetime import datetime

def excel_filter_gcf_override(
    input_file: str,
    strain_col: str = "Strain_ID",  # 存储GCA/GCF标识的列名（默认Strain_ID）
    subset: list = None,
    keep: str = "first",
    ignore_null: bool = True
) -> None:
    """
    Excel筛选工具：保留GCF_前缀行，删除GCA_前缀行（直接覆盖源文件）
    :param input_file: 输入Excel文件路径（支持.xlsx/.xls）
    :param strain_col: 存储GCA/GCF标识的列名（默认：Strain_ID）
    :param subset: 筛选后额外去重的基准列（列表格式，默认：不额外去重）
    :param keep: 重复行保留策略（"first"/"last"/False，默认：first）
    :param ignore_null: 是否忽略空值差异（默认：忽略）
    """
    # 1. 校验输入文件
    input_path = Path(input_file)
    if not input_path.exists():
        raise FileNotFoundError(f"输入文件不存在：{input_file}")
    if input_path.suffix not in [".xlsx", ".xls"]:
        raise ValueError("仅支持.xlsx和.xls格式的Excel文件")
    
    # 2. 生成源文件备份（防止数据丢失）
    backup_suffix = datetime.now().strftime("%Y%m%d%H%M%S")
    backup_file = input_path.parent / f"{input_path.stem}_backup_{backup_suffix}{input_path.suffix}"
    try:
        shutil.copy2(input_file, backup_file)
        print(f"🔧 已生成源文件备份：{backup_file}")
    except Exception as e:
        raise RuntimeError(f"生成备份文件失败：{str(e)}")
    
    # 3. 读取Excel文件
    try:
        engine = "openpyxl" if input_path.suffix == ".xlsx" else "xlrd"
        df = pd.read_excel(input_file, engine=engine)
    except Exception as e:
        raise RuntimeError(f"读取Excel失败：{str(e)}")
    
    # 4. 校验Strain_ID列是否存在
    if strain_col not in df.columns:
        raise KeyError(f"Excel中未找到列名：{strain_col}，请确认存储GCA/GCF标识的列名")
    
    # 5. 数据预处理（清洗Strain_ID列，避免前缀识别失败）
    df[strain_col] = df[strain_col].astype(str).str.strip()  # 转为字符串+去前后空格
    print(f"✅ 成功读取文件：{input_file}")
    print(f"📊 原始总行数：{len(df)}")
    
    # 6. 筛选：保留GCF_前缀行，删除GCA_前缀行
    # 统计各前缀行数
    gcf_count = df[df[strain_col].str.startswith("GCF_", na=False)].shape[0]
    gca_count = df[df[strain_col].str.startswith("GCA_", na=False)].shape[0]
    other_count = len(df) - gcf_count - gca_count  # 既不是GCF也不是GCA的行
    
    print(f"\n📋 前缀分布统计：")
    print(f"   GCF_前缀行数：{gcf_count}（将保留）")
    print(f"   GCA_前缀行数：{gca_count}（将删除）")
    print(f"   其他前缀行数：{other_count}（将保留）")
    
    # 执行筛选
    df_filtered = df[
        (df[strain_col].str.startswith("GCF_", na=False)) |  # 保留GCF_
        (~df[strain_col].str.startswith(("GCA_", "GCF_"), na=False))  # 保留其他前缀
    ].copy()
    
    print(f"\n📊 筛选后总行数：{len(df_filtered)}")
    print(f"🗑️  已删除GCA_前缀行数：{gca_count}")
    
    # 7. 可选：筛选后进行额外去重（如需去重则执行）
    if subset is not None:
        # 数据清洗（解决去重不生效问题）
        for col in df_filtered.columns:
            if df_filtered[col].dtype == "object":
                df_filtered[col] = df_filtered[col].astype(str).str.strip()
        if ignore_null:
            df_filtered = df_filtered.replace(["nan", ""], None)
        else:
            df_filtered = df_filtered.replace("nan", None)
        
        # 执行去重
        original_filtered_rows = len(df_filtered)
        df_filtered = df_filtered.drop_duplicates(
            subset=subset,
            keep=keep,
            ignore_index=True
        )
        deduplicated_rows = len(df_filtered)
        removed_dup_rows = original_filtered_rows - deduplicated_rows
        
        print(f"\n🔍 额外去重（基于列：{subset}）：")
        print(f"   去重前行数：{original_filtered_rows}")
        print(f"   去重后行数：{deduplicated_rows}")
        print(f"   移除重复行数：{removed_dup_rows}")
    
    # 8. 覆盖源文件
    try:
        df_filtered.to_excel(input_file, index=False, engine="openpyxl")
        print(f"\n✅ 已直接覆盖源文件：{input_file}")
    except Exception as e:
        # 覆盖失败时恢复备份
        shutil.copy2(backup_file, input_file)
        raise RuntimeError(f"覆盖源文件失败！已自动恢复原文件：{str(e)}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Excel筛选工具：保留GCF_前缀，删除GCA_前缀（直接覆盖源文件）")
    parser.add_argument("-i", "--input", required=True, help="输入Excel文件路径（.xlsx/.xls）")
    parser.add_argument("-c", "--strain-col", default="Strain_ID", help="存储GCA/GCF标识的列名（默认：Strain_ID）")
    parser.add_argument("-s", "--subset", nargs="+", help="筛选后额外去重的基准列（如：_serotype serotypon，默认：不去重）")
    parser.add_argument("-k", "--keep", choices=["first", "last", "False"], default="first", 
                        help="重复行保留策略（first=保留首次，last=保留最后，False=删除所有，默认：first）")
    parser.add_argument("-n", "--no-ignore-null", action="store_true", help="不忽略空值差异（默认：忽略）")
    
    args = parser.parse_args()
    
    # 处理keep参数
    keep_param = args.keep if args.keep != "False" else False
    
    # 执行筛选
    try:
        excel_filter_gcf_override(
            input_file=args.input,
            strain_col=args.strain_col,
            subset=args.subset,
            keep=keep_param,
            ignore_null=not args.no_ignore_null
        )
        print("\n🎉 处理完成！已保留GCF_前缀行，删除GCA_前缀行（源文件已更新）")
    except Exception as e:
        print(f"\n❌ 出错：{str(e)}")
