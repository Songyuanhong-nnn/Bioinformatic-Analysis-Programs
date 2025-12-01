datasets summary genome taxon "vibrio parahaemolyticus" \
--assembly-level complete \
| jq -r '
["Assembly", "访问号", "菌株", "血清型", "总长度(bp)", "GC含量(%)", "总基因数", "组装水平", "完整性(%)", "污染率(%)", "注释更新日期"],
(.reports[] | [
  .assembly_info.assembly_name,
  .accession,
  .organism.infraspecific_names.strain // "未知",
  (.assembly_info.biosample.attributes[]? | select(.name == "serotype").value) // "未知",
  .assembly_stats.total_sequence_length,
  .assembly_stats.gc_percent,
  .annotation_info.stats.gene_counts.total // "未知",
  .assembly_info.assembly_level,
  .checkm_info.completeness // "未知",
  .checkm_info.contamination // "未知",
  .annotation_info.release_date // "未知"
]) | @csv
' > vibrio_assembly_detailed.csv
