#!/usr/bin/env python3

import pandas as pd

# Create a simple DataFrame
data = {
    'Name': ['Alice', 'Bob', 'Charlie', 'David', 'Eve'],
    'Age': [25, 30, 35, 40, 45],
    'City': ['New York', 'London', 'Paris', 'Tokyo', 'Sydney'],
    'Email': ['alice@example.com', 'bob@example.com', 'charlie@example.com', 'david@example.com', 'eve@example.com'],
    'o3:k6': ['value1', 'value2', 'value3', 'value4', 'value5']
}

df = pd.DataFrame(data)

# Save to Excel file
df.to_excel('test_excel.xlsx', index=False)

print("Test Excel file created successfully: test_excel.xlsx")
