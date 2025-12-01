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

def get_existing_fna_md5(target_dir: str) -> set:
    """获取目标目录中已存在的所有FNA文件的MD5集合"""
    existing_md5 = set()
    if not os.path.exists(target_dir):
        return existing_md5
    
    for file in os.listdir(target_dir):
        if file.endswith(".fna"):
            file_path = os.path.join(target_dir, file)
            existing_md5.add(calculate_md5(file_path))
    print(f"检测到目标目录中已有 {len(existing_md5)} 个唯一FNA文件，将跳过这些文件")
    return existing_md5

def unzip_all(zip_dir: str, extract_temp_dir: str) -> None:
    """批量解压zip压缩包到临时目录（包括子文件夹中的压缩包）"""
    os.makedirs(extract_temp_dir, exist_ok=True)
    zip_files = []
    # 遍历所有子目录查找zip文件
    for root, _, files in os.walk(zip_dir):
        for file in files:
            if file.endswith(".zip"):
                zip_files.append(os.path.join(root, file))
    
    if not zip_files:
        print(f"警告：{zip_dir} 中未找到任何.zip压缩包")
        return

    for zip_file in zip_files:
        print(f"正在解压：{zip_file}")
        try:
            with zipfile.ZipFile(zip_file, "r") as zf:
                zf.extractall(extract_temp_dir)
            print(f"解压完成：{zip_file}")
        except Exception as e:
            print(f"解压失败 {zip_file}：{str(e)}")

def process_existing_folders(source_dir: str, extract_temp_dir: str) -> None:
    """处理已解压的文件夹，复制到临时目录统一处理"""
    os.makedirs(extract_temp_dir, exist_ok=True)
    
    for item in os.listdir(source_dir):
        item_path = os.path.join(source_dir, item)
        if os.path.isdir(item_path) and not item.endswith(".zip"):
            dest_path = os.path.join(extract_temp_dir, item)
            if os.path.exists(dest_path):
                shutil.rmtree(dest_path)
            try:
                shutil.copytree(item_path, dest_path)
                print(f"已复制文件夹：{item} 到临时目录")
            except Exception as e:
                print(f"复制文件夹 {item} 失败：{str(e)}")

def extract_fna(extract_temp_dir: str, target_fna_dir: str, existing_md5: set) -> List[str]:
    """从临时目录提取所有新的.fna文件（排除已存在的）到目标目录"""
    os.makedirs(target_fna_dir, exist_ok=True)
    new_fna_files = []

    for root, _, files in os.walk(extract_temp_dir):
        for file in files:
            if file.endswith(".fna"):
                source_path = os.path.join(root, file)
                # 计算当前文件MD5
                current_md5 = calculate_md5(source_path)
                
                # 检查是否已存在
                if current_md5 in existing_md5:
                    print(f"跳过已存在文件：{file}")
                    continue

                # 处理同名文件
                target_path = os.path.join(target_fna_dir, file)
                counter = 1
                while os.path.exists(target_path):
                    name, ext = os.path.splitext(file)
                    target_path = os.path.join(target_fna_dir, f"{name}_{counter}{ext}")
                    counter += 1
                
                shutil.copy2(source_path, target_path)
                new_fna_files.append(target_path)
                existing_md5.add(current_md5)  # 加入已存在集合，避免后续重复
                print(f"已提取新文件：{file} {'-> ' + os.path.basename(target_path) if counter > 1 else ''}")
    
    if not new_fna_files:
        print(f"警告：{extract_temp_dir} 中未找到新的.fna文件")
    else:
        print(f"\n共提取 {len(new_fna_files)} 个新的.fna文件到 {target_fna_dir}")
    return new_fna_files

def remove_duplicates(target_fna_dir: str) -> tuple[int, int, int]:
    """删除目标目录中的重复.fna文件，保留首次出现的文件"""
    md5_record = {}
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
    """生成FNA文件提取总结报告"""
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
        f.write(f"1. 提取统计：本次共提取 {total_extracted} 个新的.fna文件\n")
        f.write(f"2. 去重统计：删除 {total_duplicated} 个重复文件，当前保留 {total_unique} 个唯一文件\n")
        f.write(f"3. 大小统计：唯一文件总大小 {total_size_mb:.2f} MB\n")
        f.write(f"4. 存储路径：{os.path.abspath(target_fna_dir)}\n")
        f.write(f"5. 唯一文件列表：\n")
        for i, file in enumerate(unique_files, 1):
            file_size_mb = os.path.getsize(os.path.join(target_fna_dir, file)) / (1024 * 1024)
            f.write(f"   {i:2d}. {file} （{file_size_mb:.2f} MB）\n")
        f.write("=" * 50 + "\n")
    
    print(f"\n总结报告已生成：{summary_path}")

def main():
    # 路径配置（核心修改：目标目录改为 clean_genome）
    TOOLS_DIR = os.path.dirname(os.path.abspath(__file__))
    PROJECT_ROOT = os.path.dirname(TOOLS_DIR)
    RAWSOURCE_DIR = os.path.join(PROJECT_ROOT, "rawsource")
    TEMP_DIR = os.path.join(PROJECT_ROOT, "temp_process")
    TARGET_DIR = os.path.join(PROJECT_ROOT, "clean_genome")  # 已改为 clean_genome

    print("=" * 60)
    print("FNA文件增量提取工具（NCBI数据专用）")
    print(f"项目根目录：{PROJECT_ROOT}")
    print(f"原始数据目录：{RAWSOURCE_DIR}")
    print(f"临时处理目录：{TEMP_DIR}")
    print(f"目标基因组目录：{TARGET_DIR}")
    print("=" * 60 + "\n")

    # 关键改进：获取已存在的FNA文件MD5
    print("【预处理】检查目标目录中已存在的FNA文件...")
    existing_md5 = get_existing_fna_md5(TARGET_DIR)

    # 清空并重建临时目录
    if os.path.exists(TEMP_DIR):
        shutil.rmtree(TEMP_DIR)
    os.makedirs(TEMP_DIR, exist_ok=True)

    # 步骤1：处理压缩包
    print("\n【步骤1/5】解压所有压缩包...")
    unzip_all(RAWSOURCE_DIR, TEMP_DIR)

    # 步骤2：处理已解压的文件夹
    print("\n【步骤2/5】处理已解压的文件夹...")
    process_existing_folders(RAWSOURCE_DIR, TEMP_DIR)

    # 步骤3：提取新的.fna文件（跳过已存在的）
    print("\n【步骤3/5】提取新的.fna文件...")
    extracted_files = extract_fna(TEMP_DIR, TARGET_DIR, existing_md5)
    total_extracted = len(extracted_files)

    # 步骤4：删除重复文件（主要处理本次新增的重复）
    print("\n【步骤4/5】删除重复文件...")
    _, total_duplicated, total_unique = remove_duplicates(TARGET_DIR)

    # 步骤5：生成总结报告
    print("\n【步骤5/5】生成总结报告...")
    generate_summary(TARGET_DIR, total_extracted, total_duplicated, total_unique)

    # 清理临时目录
    if os.path.exists(TEMP_DIR):
        shutil.rmtree(TEMP_DIR)
        print(f"\n已自动清理临时目录：{TEMP_DIR}")

    print(f"\n流程结束！新增FNA文件已提取至 {TARGET_DIR}")

if __name__ == "__main__":
    main()