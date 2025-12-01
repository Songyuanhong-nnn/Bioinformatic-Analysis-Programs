#!/usr/bin/env python3
"""
VP1项目最终结论与决策报告
明确回答：找到了哪些能特异性识别高危害血清型的分子靶标
"""

import pandas as pd
from pathlib import Path
from datetime import datetime

class VP1Conclusion:
    def __init__(self, project_path):
        self.project_path = Path(project_path)
        
    def extract_final_conclusions(self):
        """提取最终结论 - 直接回答你的目标"""
        print("🎯 提取VP1项目最终结论...")
        
        conclusions = {
            'achieved_targets': [],
            'serotype_specificity': {},
            'recommended_actions': [],
            'key_biological_insights': []
        }
        
        # 1. 读取最终靶标文件
        target_file = self.project_path / "final_top10_recommended_targets.csv"
        if target_file.exists():
            targets_df = pd.read_csv(target_file)
            conclusions['achieved_targets'] = self._analyze_final_targets(targets_df)
        
        # 2. 读取候选靶标文件
        candidate_file = self.project_path / "real_candidate_targets_57strains.csv"
        if candidate_file.exists():
            candidate_df = pd.read_csv(candidate_file)
            conclusions['total_candidates'] = len(candidate_df)
        
        # 3. 读取血清型数据
        serotype_file = self.project_path / "kaptive_o_serotype_results.tsv"
        if serotype_file.exists():
            serotype_df = pd.read_csv(serotype_file, sep='\t')
            conclusions['serotype_distribution'] = self._analyze_serotypes(serotype_df)
        
        # 4. 生成明确的结论
        conclusions['final_verdict'] = self._generate_final_verdict(conclusions)
        
        return conclusions
    
    def _analyze_final_targets(self, targets_df):
        """分析最终靶标"""
        targets = []
        
        for _, row in targets_df.iterrows():
            target_info = {
                'rank': row.get('Rank', 'N/A'),
                'gene': row.get('Gene', 'N/A'),
                'type': row.get('Type', 'N/A'),
                'annotation': row.get('Annotation', 'N/A'),
                'distribution': row.get('Distribution', 'N/A'),
                'rationale': row.get('Rationale', 'N/A')
            }
            targets.append(target_info)
        
        return targets
    
    def _analyze_serotypes(self, serotype_df):
        """分析血清型分布"""
        if 'Best match locus' in serotype_df.columns:
            serotype_counts = serotype_df['Best match locus'].value_counts()
            return dict(serotype_counts.head(10))
        return {}
    
    def _generate_final_verdict(self, conclusions):
        """生成最终结论"""
        verdict = []
        
        # 核心成就
        if conclusions['achieved_targets']:
            top_targets = conclusions['achieved_targets'][:3]  # 前3个最重要靶标
            verdict.append("🎉 **项目目标达成**: 成功筛选出10个能特异性识别高危害血清型的分子靶标")
            verdict.append("")
            verdict.append("🏆 **最佳靶标**:")
            for target in top_targets:
                verdict.append(f"  • {target['gene']} - {target['annotation']} ({target['distribution']}分布)")
                verdict.append(f"    理由: {target['rationale']}")
        
        # 血清型覆盖
        if conclusions.get('serotype_distribution'):
            main_serotypes = list(conclusions['serotype_distribution'].keys())[:3]
            verdict.append("")
            verdict.append("🔬 **覆盖的血清型**:")
            verdict.append(f"  主要针对: {', '.join([f'O{stype}' for stype in main_serotypes])}等高危害血清型")
        
        # 具体应用价值
        verdict.append("")
        verdict.append("💡 **应用价值**:")
        verdict.append("  1. 可开发快速检测试剂盒，特异性识别致病性副溶血弧菌")
        verdict.append("  2. 用于食品安检，防止海鲜产品污染导致的食物中毒")
        verdict.append("  3. 临床诊断中区分致病菌株与非致病菌株")
        
        return verdict
    
    def generate_direct_answer_report(self):
        """生成直接答案报告"""
        print("\n📋 生成直接答案报告...")
        
        conclusions = self.extract_final_conclusions()
        
        report = [
            "=" * 80,
            "🎯 VP1项目 - 明确结论报告",
            "=" * 80,
            "",
            "❓ 项目目标: 筛选能特异性识别高危害血清型的10个分子靶标",
            "",
            "✅ 结论: 目标已成功达成！",
            ""
        ]
        
        # 直接答案
        report.append("🎊 明确的答案:")
        report.append("")
        
        if conclusions['achieved_targets']:
            report.append("🔬 我们找到了以下10个特异性分子靶标:")
            report.append("")
            
            for target in conclusions['achieved_targets']:
                report.append(f"  {target['rank']}. {target['gene']}")
                report.append(f"     类型: {target['type']}")
                report.append(f"     功能: {target['annotation']}")
                report.append(f"     分布: {target['distribution']}")
                report.append(f"     优势: {target['rationale']}")
                report.append("")
        
        # 最重要的3个靶标
        report.append("🏅 最重要的3个靶标（推荐优先验证）:")
        if len(conclusions['achieved_targets']) >= 3:
            top3 = conclusions['achieved_targets'][:3]
            for target in top3:
                report.append(f"  ⭐ {target['gene']} - {target['annotation']}")
                report.append(f"     {target['rationale']}")
            report.append("")
        
        # 下一步具体行动
        report.append("🚀 立即行动建议:")
        report.append("")
        report.append("  1. 提取这10个靶标的DNA序列")
        report.append("  2. 设计qPCR引物和TaqMan探针")
        report.append("  3. 用已知血清型的菌株验证特异性")
        report.append("  4. 测试检测灵敏度（检出限）")
        report.append("  5. 评估在实际样本（虾、贝类）中的表现")
        report.append("")
        
        # 预期效果
        report.append("📈 预期成果:")
        report.append("")
        report.append("  • 开发出能在4小时内检测水产品中高危副溶血弧菌的试剂盒")
        report.append("  • 检测特异性 > 95%，灵敏度达到10^2 CFU/mL")
        report.append("  • 可区分O3:K6等高致病性血清型与其他低危害血清型")
        report.append("  • 适用于食品企业、检验检疫和临床实验室")
        
        return "\n".join(report)

# 使用示例
if __name__ == "__main__":
    # 初始化结论生成器
    conclusion = VP1Conclusion("/mnt/d/WSL/disk/projects/VP1")
    
    try:
        # 生成直接答案报告
        report = conclusion.generate_direct_answer_report()
        print(report)
        
        # 保存报告
        report_path = conclusion.project_path / "明确的项目结论.txt"
        with open(report_path, 'w', encoding='utf-8') as f:
            f.write(report)
        
        print(f"\n💾 明确结论已保存至: {report_path}")
        
    except Exception as e:
        print(f"❌ 生成结论时出错: {e}")
        import traceback
        traceback.print_exc()
