import requests
import csv

# NCBI Assembly API查询血清型（输入GCA菌株ID，输出血清型）
def get_serotype_from_ncbi(isolate_id):
    # 提取GCA编号（如GCA_000196095.1_ASM19609v1 → GCA_000196095.1）
    gca_id = isolate_id.split('_')[0]
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=assembly&id={gca_id}&retmode=json"
    
    try:
        response = requests.get(url, timeout=10)
        data = response.json()
        # 从NCBI返回结果中提取血清型（Serotype）
        for key in data['result']:
            if key != 'uids':
                biosample = data['result'][key].get('biosample', {})
                attributes = biosample.get('attributes', [])
                for attr in attributes:
                    if attr.get('name') == 'serotype':
                        return attr.get('value', 'Unknown')
        return 'Unknown'  # 未查询到血清型
    except Exception as e:
        print(f"查询{isolate_id}血清型失败：{e}")
        return 'Unknown'

# 基于毒力基因划分危害等级（文献公认标准）
def classify_hazard_level(row):
    has_tdh = int(row['Has_TDH'])
    has_trh = int(row['Has_TRH'])
    has_t3ss2 = int(row['Has_T3SS2'])
    
    # 高危害：TDH/TRH + T3SS2（核心致病组合）
    if (has_tdh == 1 or has_trh == 1) and has_t3ss2 == 1:
        return '高危害'
    # 中危害：仅TDH/TRH 或 仅T3SS2
    elif (has_tdh == 1 or has_trh == 1) or has_t3ss2 == 1:
        return '中危害'
    # 低危害：无核心毒力基因
    else:
        return '低危害'

# 读取原始表型文件，补充血清型和危害等级
input_file = 'phenotype.txt'
output_file = 'phenotype_completed.txt'

with open(input_file, 'r', newline='', encoding='utf-8') as infile, \
     open(output_file, 'w', newline='', encoding='utf-8') as outfile:
    # 读取原始表型（制表符分隔）
    reader = csv.DictReader(infile, delimiter='\t')
    # 新增Serotype和Hazard_Level列，生成新表头
    fieldnames = ['Isolate', 'Serotype', 'Hazard_Level'] + reader.fieldnames[1:]
    writer = csv.DictWriter(outfile, fieldnames=fieldnames, delimiter='\t')
    writer.writeheader()
    
    # 遍历每一行菌株，自动补充信息
    for row in reader:
        isolate_id = row['Isolate']
        print(f"正在处理菌株：{isolate_id}")
        
        # 1. 自动查询血清型
        serotype = get_serotype_from_ncbi(isolate_id)
        # 2. 自动划分危害等级
        hazard_level = classify_hazard_level(row)
        
        # 3. 写入补充后的行
        writer.writerow({
            'Isolate': isolate_id,
            'Serotype': serotype,
            'Hazard_Level': hazard_level,
            'Has_T3SS1': row['Has_T3SS1'],
            'Has_T3SS2': row['Has_T3SS2'],
            'Has_MAM7': row['Has_MAM7'],
            'Has_TLH': row['Has_TLH'],
            'Has_TDH': row['Has_TDH'],
            'Has_TRH': row['Has_TRH'],
            'Has_Hemolysin_B_C': row['Has_Hemolysin_B_C'],
            'Has_VPA0450': row['Has_VPA0450']
        })

print(f"自动补充完成！输出文件：{output_file}")
