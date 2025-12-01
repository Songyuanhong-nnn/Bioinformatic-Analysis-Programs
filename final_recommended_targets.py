import pandas as pd

def create_final_recommendations():
    """创建最终推荐靶标清单"""
    
    recommendations = [
        {
            'Rank': 1,
            'Gene': 'group_9360',
            'Type': '分泌系统',
            'Annotation': 'Type II secretion system protein E',
            'Distribution': '46.8%',
            'Rationale': '最高分靶标，II型分泌系统关键蛋白，分布适中利于特异性检测',
            'Priority': '⭐⭐⭐⭐⭐'
        },
        {
            'Rank': 2,
            'Gene': 'group_3365', 
            'Type': '溶血素',
            'Annotation': 'Thermolabile hemolysin',
            'Distribution': '96.8%',
            'Rationale': '热不稳定溶血素，在绝大多数菌株中存在，适合作为广谱检测靶标',
            'Priority': '⭐⭐⭐⭐⭐'
        },
        {
            'Rank': 3,
            'Gene': 'tdh2',
            'Type': '经典毒力因子',
            'Annotation': 'Thermostable direct hemolysin 2',
            'Distribution': '17.1%', 
            'Rationale': '已知的副溶血弧菌关键毒力因子，文献支持充分',
            'Priority': '⭐⭐⭐⭐⭐'
        },
        {
            'Rank': 4,
            'Gene': 'vopS',
            'Type': '效应蛋白',
            'Annotation': 'Protein adenylyltransferase VopS',
            'Distribution': '69.9%',
            'Rationale': 'III型分泌系统效应蛋白，在致病机制中起关键作用',
            'Priority': '⭐⭐⭐⭐'
        },
        {
            'Rank': 5,
            'Gene': 'epsE_2',
            'Type': '分泌系统',
            'Annotation': 'Type II secretion system protein E',
            'Distribution': '54.2%',
            'Rationale': 'II型分泌系统组件，分布频率理想',
            'Priority': '⭐⭐⭐⭐'
        },
        {
            'Rank': 6,
            'Gene': 'tdh1',
            'Type': '经典毒力因子', 
            'Annotation': 'Thermostable direct hemolysin 1',
            'Distribution': '10.6%',
            'Rationale': '热稳定溶血素1型，可能代表特定致病亚型',
            'Priority': '⭐⭐⭐⭐'
        },
        {
            'Rank': 7,
            'Gene': 'ureC',
            'Type': '代谢毒力',
            'Annotation': 'Urease subunit alpha', 
            'Distribution': '20.8%',
            'Rationale': '脲酶系统，与细菌在宿主体内存活相关',
            'Priority': '⭐⭐⭐'
        },
        {
            'Rank': 8,
            'Gene': 'group_6946',
            'Type': '溶血素',
            'Annotation': 'Hemolysin, chromosomal',
            'Distribution': '3.7%',
            'Rationale': '染色体溶血素，分布稀有，可能代表高毒力亚型',
            'Priority': '⭐⭐⭐'
        },
        {
            'Rank': 9,
            'Gene': 'hlyA',
            'Type': '溶血素',
            'Annotation': 'Hemolysin, chromosomal',
            'Distribution': '4.6%',
            'Rationale': '经典溶血素A基因，虽然分布率低但毒力明确',
            'Priority': '⭐⭐⭐'
        },
        {
            'Rank': 10,
            'Gene': 'group_2253',
            'Type': '分泌系统',
            'Annotation': 'Type 3 secretion system secretin',
            'Distribution': '6.5%',
            'Rationale': 'III型分泌系统secretin蛋白，与高毒力相关',
            'Priority': '⭐⭐⭐'
        }
    ]
    
    df = pd.DataFrame(recommendations)
    df.to_csv('final_top10_recommended_targets.csv', index=False)
    
    print("=== 🎯 最终推荐的前10个靶标 ===")
    print("这些靶标基于功能重要性、分布特异性和文献支持综合选择\n")
    
    for target in recommendations:
        print(f"{target['Priority']} 第{target['Rank']}名: {target['Gene']}")
        print(f"   类型: {target['Type']}")
        print(f"   功能: {target['Annotation']}")
        print(f"   分布: {target['Distribution']} (101/216)")
        print(f"   理由: {target['Rationale']}")
        print()
    
    return df

def generate_validation_plan():
    """生成验证实验方案"""
    
    print("=== 🔬 实验验证方案 ===")
    print("\n1. 引物设计:")
    print("   - 为每个靶标设计特异性PCR引物")
    print("   - 引物长度: 18-22 bp, Tm值: 55-65°C")
    print("   - 扩增片段: 200-500 bp")
    
    print("\n2. 验证策略:")
    print("   - 使用已知血清型的菌株进行PCR验证")
    print("   - 高危血清型 vs 低危血清型比较")
    print("   - 检测灵敏度和特异性评估")
    
    print("\n3. 功能验证:")
    print("   - 基因敲除验证毒力功能")
    print("   - 表达分析确认转录水平")
    print("   - 蛋白定位研究")

if __name__ == "__main__":
    create_final_recommendations()
    generate_validation_plan()
    
    print("\n=== 📋 下一步工作建议 ===")
    print("1. 立即开始前3个靶标(group_9360, group_3365, tdh2)的实验验证")
    print("2. 对高分布率靶标开发快速检测方法") 
    print("3. 对低分布率靶标研究其与特定致病型的关联")
    print("4. 结合临床数据验证这些靶标的流行病学意义")
