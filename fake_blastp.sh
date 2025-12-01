#!/bin/bash
# 伪装 blastp 版本为 2.21.0（高于 Prokka 要求的 2.2）
if [ "$1" = "-version" ]; then
  echo "blastp: 2.21.0+"
  echo " Package: blast 2.21.0, build Jul  1 2025 08:59:18"
else
  /opt/ncbi-blast-2.17.0+/bin/blastp "$@"
fi

