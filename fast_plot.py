import pandas as pd
import matplotlib.pyplot as plt

plt.rcParams['font.sans-serif'] = ['Arial']
df = pd.read_csv("/mnt/d/WSL/disk/projects/VP1/FAST_serotype_specific_virulence.csv")
sero_count = df['血清型'].value_counts()

fig, ax = plt.subplots(figsize=(12, 6))
sero_count.plot(kind='bar', ax=ax, color='steelblue')
ax.set_title('Serotype-Specific Virulence Genes', fontsize=14, fontweight='bold')
ax.set_xlabel('O Serotype')
ax.set_ylabel('Number of Specific Genes')
ax.tick_params(axis='x', rotation=45)
plt.tight_layout()
plt.savefig("/mnt/d/WSL/disk/projects/VP1/FAST_virulence_plot.png", dpi=300)
print("✅ 可视化图表生成完成")
