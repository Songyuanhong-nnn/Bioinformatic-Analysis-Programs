import os
import hashlib
import sys

def get_file_hash(file_path, block_size=65536):
    """计算文件的MD5哈希值（判断内容是否相同）"""
    hash_md5 = hashlib.md5()
    try:
        with open(file_path, "rb") as f:
            for block in iter(lambda: f.read(block_size), b""):
                hash_md5.update(block)
        return hash_md5.hexdigest()
    except Exception as e:
        print(f"跳过无法读取的文件：{file_path}（错误：{e}）")
        return None

def delete_duplicates():
    current_dir = os.getcwd()  # 获取程序所在的当前文件夹
    file_hashes = {}  # 存储哈希值与文件路径的映射：{哈希值: [文件路径列表]}
    program_name = os.path.basename(sys.argv[0])  # 自身程序的文件名

    # 扫描当前文件夹中的所有文件（不包括子目录）
    for filename in os.listdir(current_dir):
        file_path = os.path.join(current_dir, filename)
        # 跳过自身程序和目录
        if filename == program_name or os.path.isdir(file_path):
            continue
        # 计算文件哈希值
        file_hash = get_file_hash(file_path)
        if file_hash:
            if file_hash not in file_hashes:
                file_hashes[file_hash] = []
            file_hashes[file_hash].append(file_path)

    # 删除重复文件（保留第一个，其余彻底删除）
    deleted_count = 0
    for hash_val, file_paths in file_hashes.items():
        if len(file_paths) > 1:
            # 保留第一个文件，删除其余
            for path in file_paths[1:]:
                try:
                    os.remove(path)  # 直接彻底删除（不进回收站）
                    print(f"已删除重复文件：{path}")
                    deleted_count += 1
                except Exception as e:
                    print(f"删除失败：{path}（错误：{e}）")

    print(f"\n操作完成，共删除 {deleted_count} 个重复文件")

    # 自我销毁（删除程序自身）
    try:
        os.remove(os.path.abspath(sys.argv[0]))
        print("程序已自动删除自身")
    except Exception as e:
        print(f"程序自我删除失败：{e}")

if __name__ == "__main__":
    delete_duplicates()