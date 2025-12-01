import os
import shutil
import zipfile
import hashlib
import argparse
from typing import List, Dict, Optional
from pathlib import Path

# ===================== 通用配置（可根据实际需求调整）=====================
SUPPORTED_ARCHIVES = (".zip", ".gz", ".tar.gz", ".tgz")  # 支持的压缩格式（扩展通用）
SUPPORTED_EXTENSIONS = (".fna", ".fasta", ".fa")  # 支持的序列文件格式（扩展通用）
TEMP_DIR_NAME = "temp_process"  # 临时目录名称
SUMMARY_FILE_NAME = "fna_extract_summary.txt"  # 总结报告名称

# ===================== 通用工具函数 =====================
def calculate_md5(file_path: str) -> str:
    """通用MD5计算函数，兼容大文件"""
    md5_hash = hashlib.md5()
    try:
        with open(file_path, "rb") as f:
            for chunk in iter(lambda: f.read(4096), b""):
                md5_hash.update(chunk)
        return md5_hash.hexdigest()
    except Exception as e:
        print(f"⚠️  计算文件MD5失败：{file_path} -> {str(e)}")
        return ""

def validate_dir_structure(root_dir: str) -> Dict[str, str]:
    """
    自动校验目录结构，返回核心目录路径
    核心约定：root_dir 下必须包含 rawsource 文件夹，脚本所在目录为 tools
    """
    # 自动识别核心目录（兼容脚本放在 tools 或其他位置，只要 root_dir 正确）
    dirs = {
        "rawsource": os.path.join(root_dir, "rawsource"),
        "clean_genome": os.path.join(root_dir, "clean_genome"),
        "temp": os.path.join(root_dir, TEMP_DIR_NAME)
    }

    # 校验原始数据目录是否存在
    if not os.path.exists(dirs["rawsource"]):
        raise FileNotFoundError(f"❌ 原始数据目录不存在：{dirs['rawsource']}\n请确保 root_dir 下包含 'rawsource' 文件夹")

    # 自动创建目标目录和临时目录
    for dir_path in dirs.values():
        os.makedirs(dir_path, exist_ok=True)
    
    print("✅ 目录结构校验通过：")
    for name, path in dirs.items():
        print(f"  - {name}: {path}")
    return dirs

def get_existing_file_md5(target_dir: str) -> set:
    """通用函数：获取目标目录中所有支持格式文件的MD5集合"""
    existing_md5 = set()
    for root, _, files in os.walk(target_dir):
        for file in files:
            if file.lower().endswith(SUPPORTED_EXTENSIONS):
                file_path = os.path.join(root, file)
                md5_val = calculate_md5(file_path)
                if md5_val:
                    existing_md5.add(md5_val)
    print(f"✅ 检测到目标目录中已有 {len(existing_md5)} 个唯一序列文件，将跳过重复文件")
    return existing_md5

def extract_archive(archive_path: str, extract_dir: str) -> None:
    """通用解压函数：支持 zip/gz/tar.gz/tgz 格式"""
    archive_ext = os.path.splitext(archive_path.lower())
    try:
        # 处理 zip 格式
        if archive_ext[-1] == ".zip":
            with zipfile.ZipFile(archive_path, "r") as zf:
                zf.extractall(extract_dir)
        # 处理 gz/tar.gz/tgz 格式
        elif archive_ext[-1] == ".gz" or (len(archive_ext) > 1 and archive_ext[-2] == ".tar"):
            import tarfile  # 延迟导入，避免未使用时加载
            mode = "r:gz" if archive_path.lower().endswith((".tar.gz", ".tgz")) else "r:gz"
            with tarfile.open(archive_path, mode) as tf:
                tf.extractall(extract_dir)
        else:
            print(f"⚠️  不支持的压缩格式：{archive_path}")
            return
        print(f"✅ 解压成功：{os.path.basename(archive_path)}")
    except Exception as e:
        print(f"❌ 解压失败 {os.path.basename(archive_path)}：{str(e)}")

def find_all_archives(source_dir: str) -> List[str]:
    """通用函数：递归查找所有支持的压缩包（包括子文件夹）"""
    archives = []
    for root, _, files in os.walk(source_dir):
        for file in files:
            if file.lower().endswith(SUPPORTED_ARCHIVES):
                archives.append(os.path.join(root, file))
    return archives

def copy_existing_folders(source_dir: str, target_dir: str) -> None:
    """通用函数：复制已解压的文件夹（跳过压缩包）"""
    for item in os.listdir(source_dir):
        item_path = os.path.join(source_dir, item)
        # 跳过压缩包，只处理文件夹
        if os.path.isdir(item_path) and not any(item.lower().endswith(ext) for ext in SUPPORTED_ARCHIVES):
            dest_path = os.path.join(target_dir, item)
            # 避免重复复制，先删除已存在的同名文件夹
            if os.path.exists(dest_path):
                shutil.rmtree(dest_path)
            try:
                shutil.copytree(item_path, dest_path, ignore=shutil.ignore_patterns(*[f"*{ext}" for ext in SUPPORTED_ARCHIVES]))
                print(f"✅ 复制文件夹成功：{item}")
            except Exception as e:
                print(f"❌ 复制文件夹失败 {item}：{str(e)}")

def extract_target_files(extract_dir: str, target_dir: str, existing_md5: set) -> List[str]:
    """通用函数：提取所有支持格式的文件（去重、同名处理）"""
    new_files = []
    for root, _, files in os.walk(extract_dir):
        for file in files:
            if file.lower().endswith(SUPPORTED_EXTENSIONS):
                source_path = os.path.join(root, file)
                current_md5 = calculate_md5(source_path)
                
                # 跳过已存在的文件
                if not current_md5 or current_md5 in existing_md5:
                    print(f"⚠️  跳过重复文件：{file}")
                    continue
                
                # 处理同名文件（添加序号）
                target_path = os.path.join(target_dir, file)
                counter = 1
                while os.path.exists(target_path):
                    name, ext = os.path.splitext(file)
                    target_path = os.path.join(target_dir, f"{name}_{counter}{ext}")
                    counter += 1
                
                # 复制文件（保留元数据）
                shutil.copy2(source_path, target_path)
                new_files.append(target_path)
                existing_md5.add(current_md5)
                print(f"✅ 提取新文件：{os.path.basename(source_path)} {'-> ' + os.path.basename(target_path) if counter > 1 else ''}")
    
    if new_files:
        print(f"\n✅ 本次共提取 {len(new_files)} 个新文件到 {target_dir}")
    else:
        print(f"\n⚠️  未找到新的序列文件")
    return new_files

def remove_duplicate_files(target_dir: str) -> tuple[int, int, int]:
    """通用去重函数：删除目标目录中重复的支持格式文件"""
    md5_record = {}
    duplicate_count = 0
    total_count = 0

    for root, _, files in os.walk(target_dir):
        for file in files:
            if file.lower().endswith(SUPPORTED_EXTENSIONS):
                total_count += 1
                file_path = os.path.join(root, file)
                md5_val = calculate_md5(file_path)

                if md5_val in md5_record:
                    try:
                        os.remove(file_path)
                        duplicate_count += 1
                        print(f"❌ 删除重复文件：{file}（与 {os.path.basename(md5_record[md5_val])} 重复）")
                    except Exception as e:
                        print(f"⚠️  删除重复文件失败 {file}：{str(e)}")
                elif md5_val:
                    md5_record[md5_val] = file_path
    
    unique_count = total_count - duplicate_count
    print(f"\n✅ 去重完成：共处理 {total_count} 个文件，删除 {duplicate_count} 个重复文件，保留 {unique_count} 个唯一文件")
    return total_count, duplicate_count, unique_count

def generate_summary_report(target_dir: str, total_extracted: int, total_duplicated: int, total_unique: int) -> None:
    """通用总结报告函数：生成详细的提取报告"""
    summary_path = os.path.join(target_dir, SUMMARY_FILE_NAME)
    total_size = 0
    unique_files = []

    # 统计唯一文件信息
    for root, _, files in os.walk(target_dir):
        for file in files:
            if file.lower().endswith(SUPPORTED_EXTENSIONS):
                file_path = os.path.join(root, file)
                unique_files.append((file, os.path.getsize(file_path), root))
                total_size += os.path.getsize(file_path)
    
    # 单位转换（字节 -> MB）
    total_size_mb = total_size / (1024 * 1024)

    # 写入报告
    with open(summary_path, "w", encoding="utf-8") as f:
        f.write("=" * 60 + "\n")
        f.write("序列文件提取与去重总结报告\n")
        f.write("=" * 60 + "\n")
        f.write(f"报告生成时间：{os.popen('date').read().strip() if os.name != 'nt' else os.popen('date /t && time /t').read().strip()}\n")
        f.write(f"1. 提取统计：本次共提取 {total_extracted} 个新序列文件\n")
        f.write(f"2. 去重统计：删除 {total_duplicated} 个重复文件，当前保留 {total_unique} 个唯一文件\n")
        f.write(f"3. 大小统计：唯一文件总大小 {total_size_mb:.2f} MB\n")
        f.write(f"4. 存储路径：{os.path.abspath(target_dir)}\n")
        f.write(f"5. 支持格式：{', '.join(SUPPORTED_EXTENSIONS)}\n")
        f.write(f"6. 唯一文件列表（共 {len(unique_files)} 个）：\n")
        for i, (file, size, root) in enumerate(unique_files, 1):
            size_mb = size / (1024 * 1024)
            f.write(f"   {i:3d}. 文件名：{file} | 大小：{size_mb:.2f} MB | 路径：{os.path.relpath(root, target_dir)}\n")
        f.write("=" * 60 + "\n")
    
    print(f"\n✅ 总结报告已生成：{summary_path}")

def clean_temp_dir(temp_dir: str) -> None:
    """通用清理函数：安全删除临时目录"""
    if os.path.exists(temp_dir):
        try:
            shutil.rmtree(temp_dir)
            print(f"\n✅ 已清理临时目录：{temp_dir}")
        except Exception as e:
            print(f"⚠️  清理临时目录失败：{str(e)}")

# ===================== 主函数（核心逻辑）=====================
def main(root_dir: Optional[str] = None):
    """
    主函数：支持两种运行方式
    1. 直接运行：自动识别脚本所在目录的上级目录为 root_dir
    2. 命令行指定：通过 --root 参数指定 root_dir
    """
    # 自动识别 root_dir（如果未指定）
    if not root_dir:
        # 脚本所在目录（tools文件夹）
        script_dir = os.path.dirname(os.path.abspath(__file__))
        # root_dir 为 tools 的上级目录（包含 rawsource 和 clean_genome）
        root_dir = os.path.dirname(script_dir)
    
    print("=" * 70)
    print("🎯 通用序列文件增量提取工具（支持 NCBI 数据）")
    print(f"📁 根目录：{root_dir}")
    print("=" * 70 + "\n")

    try:
        # 步骤1：校验目录结构
        print("【步骤1/6】校验目录结构...")
        dirs = validate_dir_structure(root_dir)
        rawsource_dir = dirs["rawsource"]
        clean_dir = dirs["clean_genome"]
        temp_dir = dirs["temp"]

        # 步骤2：获取已存在文件的MD5（增量提取核心）
        print("\n【步骤2/6】检测已存在的序列文件...")
        existing_md5 = get_existing_file_md5(clean_dir)

        # 步骤3：解压所有压缩包
        print("\n【步骤3/6】解压所有支持的压缩包...")
        archives = find_all_archives(rawsource_dir)
        if archives:
            print(f"找到 {len(archives)} 个压缩包，开始解压...")
            for archive in archives:
                extract_archive(archive, temp_dir)
        else:
            print("⚠️  未找到任何支持的压缩包，跳过解压步骤")

        # 步骤4：复制已解压的文件夹
        print("\n【步骤4/6】复制已解压的文件夹...")
        copy_existing_folders(rawsource_dir, temp_dir)

        # 步骤5：提取新文件（去重、同名处理）
        print("\n【步骤5/6】提取新的序列文件...")
        new_files = extract_target_files(temp_dir, clean_dir, existing_md5)
        total_extracted = len(new_files)

        # 步骤6：去重 + 生成报告
        print("\n【步骤6/6】去重并生成总结报告...")
        total_count, total_duplicated, total_unique = remove_duplicate_files(clean_dir)
        generate_summary_report(clean_dir, total_extracted, total_duplicated, total_unique)

        # 清理临时目录
        clean_temp_dir(temp_dir)

        print(f"\n🎉 所有流程完成！最终结果保存在：{clean_dir}")

    except Exception as e:
        print(f"\n❌ 程序运行失败：{str(e)}")
        # 异常时也尝试清理临时目录
        if 'temp_dir' in locals() and os.path.exists(temp_dir):
            clean_temp_dir(temp_dir)
        exit(1)

# ===================== 命令行支持（增强通用性）=====================
if __name__ == "__main__":
    # 解析命令行参数
    parser = argparse.ArgumentParser(description="通用序列文件增量提取工具（支持 FNA/FASTA/FA 格式，自动解压压缩包）")
    parser.add_argument("--root", type=str, help="指定根目录（包含 rawsource 和 clean_genome 的目录，可选）")
    args = parser.parse_args()

    # 运行主函数
    main(root_dir=args.root)