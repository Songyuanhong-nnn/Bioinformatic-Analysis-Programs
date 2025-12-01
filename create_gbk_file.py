#!/usr/bin/env python3

# Create a GBK encoded test file
content = "这是第一行，包含测试内容\n这是第二行，没有特定内容\n这是第三行，再次包含测试内容\n这是第四行，普通文本\n这是第五行，包含TEST内容（大写）\n"

with open('test_gbk_encoded.txt', 'w', encoding='gbk') as f:
    f.write(content)

print("GBK encoded file created successfully: test_gbk_encoded.txt")
