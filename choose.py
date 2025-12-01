import os
import shutil
import zipfile
import hashlib
from typing import List

def calculate_md5(file_path: str) -> str:
    """计算文件MD5值，用于精准判断重复文件"""
    md5_hash = hashlib.md5()
    with open(file_path, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            md5_hash.update(chunk)
    return md5_hash.hexdigest()

def unzip_all(zip_dir: str, extract_temp_dir: str) -> None:
    """批量解压zip压缩包到临时目录"""
    os.makedirs(extract_temp_dir, exist_ok=True)
    zip_files = [f for f in os.listdir(zip_dir) if f.endswith(".zip")]
    
    if not zip_files:
        print(f"警告：{zip_dir} 中未找到任何.zip压缩包")
        return

    for zip_file in zip_files:
        zip_file_path = os.path.join(zip_dir, zip_file)
        print(f"正在解压：{zip_file}")
        try:
            with zipfile.ZipFile(zip_file_path, "r") as zf:
                zf.extractall(extract_temp_dir)
            print(f"解压完成：{zip_file}")
        except Exception as e:
            print(f"解压失败 {zip_file}：{str(e)}")

def extract_fna(extract_temp_dir: str, target_fna_dir: str) -> List[str]:
    """从临时解压目录提取所有.fna文件到目标目录（raw_extract）"""
    os.makedirs(target_fna_dir, exist_ok=True)
    fna_file_paths = []

    for root, _, files in os.walk(extract_temp_dir):
        for file in files:
            if file.endswith(".fna"):
                source_path = os.path.join(root, file)
                target_path = os.path.join(target_fna_dir, file)
                shutil.copy2(source_path, target_path)
                fna_file_paths.append(target_path)
                print(f"已提取：{file}")
    
    if not fna_file_paths:
        print(f"警告：{extract_temp_dir} 中未找到任何.fna文件")
    else:
        print(f"\n共提取 {len(fna_file_paths)} 个.fna文件到 {target_fna_dir}")
    return fna_file_paths

def remove_duplicates(target_fna_dir: str) -> tuple[int, int, int]:
    """删除目标目录中的重复.fna文件，保留首次出现的文件"""
    md5_record = {}  # 存储MD5与首次出现文件路径的映射
    duplicate_count = 0
    total_count = 0

    for file in os.listdir(target_fna_dir):
        if file.endswith(".fna"):
            total_count += 1
            file_path = os.path.join(target_fna_dir, file)
            md5_val = calculate_md5(file_path)

            if md5_val in md5_record:
                os.remove(file_path)
                duplicate_count += 1
                print(f"删除重复文件：{file}（与 {os.path.basename(md5_record[md5_val])} 重复）")
            else:
                md5_record[md5_val] = file_path
    
    unique_count = total_count - duplicate_count
    print(f"\n去重完成：处理 {total_count} 个文件，删除 {duplicate_count} 个重复文件，保留 {unique_count} 个唯一文件")
    return total_count, duplicate_count, unique_count

def generate_summary(target_fna_dir: str, total_extracted: int, total_duplicated: int, total_unique: int) -> None:
    """生成FNA文件提取总结报告，保存在raw_extract目录下"""
    summary_path = os.path.join(target_fna_dir, "fna_extract_summary.txt")
    total_size = 0
    unique_files = [f for f in os.listdir(target_fna_dir) if f.endswith(".fna")]
    
    for file in unique_files:
        total_size += os.path.getsize(os.path.join(target_fna_dir, file))
    total_size_mb = total_size / (1024 * 1024)

    with open(summary_path, "w", encoding="utf-8") as f:
        f.write("=" * 50 + "\n")
        f.write("FNA文件提取与去重总结报告\n")
        f.write("=" * 50 + "\n")
        f.write(f"1. 提取统计：共提取 {total_extracted} 个.fna文件\n")
        f.write(f"2. 去重统计：删除 {total_duplicated} 个重复文件，保留 {total_unique} 个唯一文件\n")
        f.write(f"3. 大小统计：唯一文件总大小 {total_size_mb:.2f} MB\n")
        f.write(f"4. 存储路径：{os.path.abspath(target_fna_dir)}\n")
        f.write(f"5. 唯一文件列表：\n")
        for i, file in enumerate(unique_files, 1):
            file_size_mb = os.path.getsize(os.path.join(target_fna_dir, file)) / (1024 * 1024)
            f.write(f"   {i:2d}. {file} （{file_size_mb:.2f} MB）\n")
        f.write("=" * 50 + "\n")
    
    print(f"\n总结报告已生成：{summary_path}")

def main():
    # -------------------------- 核心路径配置（无需修改，适配TOOLS文件夹结构） --------------------------
    # 脚本所在目录（TOOLS文件夹）的上级目录，即包含rawsource、raw_extract、TOOLS的项目根目录
    PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    ZIP_DIR = os.path.join(PROJECT_ROOT, "rawsource")  # 压缩包目录（rawsource）
    EXTRACT_TEMP_DIR = os.path.join(PROJECT_ROOT, "temp_extract")  # 临时解压目录（项目根目录下）
    TARGET_FNA_DIR = os.path.join(PROJECT_ROOT, "raw_extract")  # 目标FNA目录（raw_extract）
    # ---------------------------------------------------------------------------------------------------

    print("=" * 60)
    print("FNA文件提取工具（适配TOOLS文件夹结构）")
    print(f"项目根目录：{PROJECT_ROOT}")
    print(f"压缩包目录：{ZIP_DIR}")
    print(f"临时解压目录：{EXTRACT_TEMP_DIR}")
    print(f"目标FNA目录：{TARGET_FNA_DIR}")
    print("=" * 60 + "\n")

    # 步骤1：批量解压压缩包
    print("【步骤1/4】批量解压压缩包...")
    unzip_all(ZIP_DIR, EXTRACT_TEMP_DIR)

    # 步骤2：提取.fna文件到raw_extract
    print("\n【步骤2/4】提取.fna文件...")
    extracted_files = extract_fna(EXTRACT_TEMP_DIR, TARGET_FNA_DIR)
    total_extracted = len(extracted_files)

    # 步骤3：删除重复文件
    print("\n【步骤3/4】删除重复文件...")
    _, total_duplicated, total_unique = remove_duplicates(TARGET_FNA_DIR)

    # 步骤4：生成总结报告
    print("\n【步骤4/4】生成总结报告...")
    generate_summary(TARGET_FNA_DIR, total_extracted, total_duplicated, total_unique)

    print(f"\n流程结束！临时解压目录 {EXTRACT_TEMP_DIR} 可手动删除以释放空间。")

if __name__ == "__main__":
    main()
