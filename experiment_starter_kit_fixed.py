import pandas as pd

def create_primer_sequences():
    """创建示例引物序列（需要根据实际序列调整）"""
    
    print("=== 🧬 引物序列示例 ==")
    print("注意: 这些是示例序列，需要根据实际基因序列调整")
    print("\n推荐使用以下工具进行引物设计:")
    print("   - Primer3 (https://primer3.org/)")
    print("   - NCBI Primer-BLAST (https://www.ncbi.nlm.nih.gov/tools/primer-blast/)")
    print("   - OligoCalc (用于引物参数计算)")
    
    primer_examples = {
        'group_9360': {
            'forward': 'ATGAGCGTCAACAGCCTGA',  # 示例序列
            'reverse': 'TCAGCTTGTCGATCGCTAG',
            'product_size': '350 bp',
            'notes': '靶向分泌系统蛋白E的保守区域'
        },
        'group_3365': {
            'forward': 'GTCAACGAGCTGTTCATCG', 
            'reverse': 'CAGTTCGATCAGCTCGATG',
            'product_size': '300 bp',
            'notes': '靶向热不稳定溶血素活性中心'
        },
        'tdh2': {
            'forward': 'AGCGTCTACAGCCTGAACG',
            'reverse': 'TCGATCGTAGCTCGATAGC', 
            'product_size': '250 bp',
            'notes': '经典毒力因子特异性区域'
        }
    }
    
    print(f"\n🎯 前3个优先靶标的引物设计示例:")
    for gene, info in primer_examples.items():
        print(f"\n📍 {gene}:")
        print(f"   正向: 5'-{info['forward']}-3'")
        print(f"   反向: 5'-{info['reverse']}-3'") 
        print(f"   产物: {info['product_size']}")
        print(f"   说明: {info['notes']}")

def create_pcr_protocol():
    """创建PCR实验方案"""
    
    print("\n" + "="*50)
    print("=== 🔬 标准PCR反应体系 ===")
    
    protocol = {
        '反应组分': {
            '模板DNA': '1-100 ng',
            '正向引物(10 μM)': '1 μL', 
            '反向引物(10 μM)': '1 μL',
            '2× PCR Master Mix': '12.5 μL',
            'ddH₂O': '至25 μL'
        },
        'PCR程序': {
            '预变性': '95°C, 5 min',
            '变性': '95°C, 30 s', 
            '退火': '55-65°C, 30 s (梯度优化)',
            '延伸': '72°C, 1 min/kb',
            '循环数': '35 cycles',
            '最终延伸': '72°C, 5 min'
        }
    }
    
    print("\n📋 反应体系:")
    for component, amount in protocol['反应组分'].items():
        print(f"   • {component}: {amount}")
    
    print(f"\n🔄 PCR程序:")
    for step, condition in protocol['PCR程序'].items():
        print(f"   • {step}: {condition}")

def create_validation_checklist():
    """创建验证检查清单"""
    
    print("\n" + "="*50)
    print("=== ✅ 实验验证检查清单 ===")
    
    checklist = [
        {
            '阶段': '准备阶段',
            'tasks': [
                '✓ 确认靶标基因序列',
                '✓ 设计特异性引物', 
                '✓ 订购引物和试剂',
                '✓ 准备测试菌株',
                '✓ 提取高质量DNA'
            ]
        },
        {
            '阶段': '优化阶段', 
            'tasks': [
                '✓ 梯度PCR确定退火温度',
                '✓ 优化Mg²⁺浓度',
                '✓ 验证引物特异性',
                '✓ 建立阳性对照'
            ]
        },
        {
            '阶段': '验证阶段',
            'tasks': [
                '✓ 测试不同菌株',
                '✓ 评估检测灵敏度', 
                '✓ 验证产物序列',
                '✓ 统计分析结果'
            ]
        }
    ]
    
    for stage in checklist:
        print(f"\n📁 {stage['阶段']}:")
        for task in stage['tasks']:
            print(f"   {task}")

def create_immediate_action_plan():
    """创建立即行动计划"""
    
    print("\n" + "="*50)
    print("=== 🚀 立即行动计划 ===")
    
    action_plan = [
        {
            '时间': '今天',
            '行动': [
                '查找group_9360、group_3365、tdh2的实际基因序列',
                '使用Primer3设计特异性引物',
                '准备引物订购清单'
            ]
        },
        {
            '时间': '明天',
            '行动': [
                '订购前3个靶标的引物',
                '准备测试菌株清单',
                '准备DNA提取试剂'
            ]
        },
        {
            '时间': '本周内',
            '行动': [
                '提取测试菌株DNA',
                '接收引物并稀释至工作浓度',
                '准备PCR实验材料'
            ]
        },
        {
            '时间': '下周',
            '行动': [
                '开始PCR条件优化',
                '进行梯度退火温度测试',
                '验证引物特异性'
            ]
        }
    ]
    
    for period in action_plan:
        print(f"\n⏰ {period['时间']}:")
        for action in period['行动']:
            print(f"   • {action}")

if __name__ == "__main__":
    create_primer_sequences()
    create_pcr_protocol() 
    create_validation_checklist()
    create_immediate_action_plan()
    
    print("\n" + "="*50)
    print("=== 🎯 实验启动成功！ ===")
    print("\n你现在拥有:")
    print("   🧪 明确的靶标清单")
    print("   🧬 引物设计指导") 
    print("   🔬 标准实验方案")
    print("   ✅ 验证检查清单")
    print("   🚀 立即行动计划")
    
    print(f"\n💡 最关键的第一步:")
    print("   立即查找这3个基因的实际序列:")
    print("   1. group_9360")
    print("   2. group_3365") 
    print("   3. tdh2")
    
    print(f"\n🔍 序列查找命令:")
    print("   find . -name 'pan_genome_reference.fa'")
    print("   grep -r 'group_9360' bakta_annotations/ prokka_annotations/")
