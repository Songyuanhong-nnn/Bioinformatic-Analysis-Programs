#!/usr/bin/env python3
import argparse
import sys
import os

def main():
    parser = argparse.ArgumentParser(
        description='Filter lines in a file that contain a specific string',
        epilog='Example usage:\n  %(prog)s input.txt "search term"\n  %(prog)s input.txt "test" --case-sensitive\n  %(prog)s input.txt "keyword" -o output.txt\n  %(prog)s input.txt "内容" --encoding gbk\n  %(prog)s input.xlsx "search term"\n  %(prog)s input.csv "search term"\n  %(prog)s input.txt "search" --replace "search:replace"',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('input_file', help='Path to the input file (text, CSV, Excel)')
    parser.add_argument('search_content', help='String to search for in each line/cell')
    parser.add_argument('-o', '--output', help='Path to write filtered results (default: print to console)')
    # Simplified case handling: always case-insensitive
    # Removed complex case options to avoid confusion
    parser.add_argument('--case-sensitive', action='store_true', help='Force case-sensitive search (optional)')
    parser.add_argument('--encoding', default='utf-8', help='File encoding (default: utf-8, try gbk for Chinese files)')
    parser.add_argument('--replace', help='Replace content in matched lines, format: "old:new"')
    
    args = parser.parse_args()
    
    file_ext = os.path.splitext(args.input_file)[1].lower()
    lines = []
    
    # Handle different file types
    if file_ext in ['.xlsx', '.xls']:
        # Excel file handling
        try:
            import pandas as pd
            # Read Excel file
            df = pd.read_excel(args.input_file)
            # Convert DataFrame to list of strings (each line is a row)
            # First add the header row
            header_line = '\t'.join(df.columns.tolist()) + '\n'
            lines.append(header_line)
            # Then add data rows
            for index, row in df.iterrows():
                # Convert row to string, joining all cells with tab separator
                line = '\t'.join([str(cell) for cell in row.values]) + '\n'
                lines.append(line)
        except ImportError:
            print(f"Error: Failed to import pandas. Please install it with 'pip install pandas openpyxl xlrd'.", file=sys.stderr)
            sys.exit(1)
        except Exception as e:
            print(f"Error reading Excel file '{args.input_file}': {e}", file=sys.stderr)
            sys.exit(1)
    elif file_ext == '.csv':
        # CSV file handling
        try:
            import pandas as pd
            # Read CSV file with specified encoding
            df = pd.read_csv(args.input_file, encoding=args.encoding)
            # Convert DataFrame to list of strings
            # First add the header row
            header_line = '\t'.join(df.columns.tolist()) + '\n'
            lines.append(header_line)
            # Then add data rows
            for index, row in df.iterrows():
                line = '\t'.join([str(cell) for cell in row.values]) + '\n'
                lines.append(line)
        except ImportError:
            print(f"Error: Failed to import pandas. Please install it with 'pip install pandas'.", file=sys.stderr)
            sys.exit(1)
        except Exception as e:
            print(f"Error reading CSV file '{args.input_file}': {e}", file=sys.stderr)
            sys.exit(1)
    else:
        # Text file handling
        # Check if file is likely a binary file (excluding supported formats)
        binary_extensions = ['.doc', '.docx', '.pdf', '.zip', '.rar', '.7z', '.exe', '.dll', '.png', '.jpg', '.jpeg', '.gif', '.bmp', '.mp3', '.mp4', '.avi']
        
        if file_ext in binary_extensions:
            print(f"Error: File '{args.input_file}' is a binary file ({file_ext}), not supported.", file=sys.stderr)
            print(f"Supported formats: text files (.txt, .log, etc.), CSV (.csv), Excel (.xlsx, .xls)", file=sys.stderr)
            sys.exit(1)
        
        try:
            with open(args.input_file, 'r', encoding=args.encoding) as f:
                lines = f.readlines()
        except FileNotFoundError:
            print(f"Error: File '{args.input_file}' not found.", file=sys.stderr)
            sys.exit(1)
        except UnicodeDecodeError:
            print(f"Error: Failed to decode file '{args.input_file}' using encoding '{args.encoding}'.", file=sys.stderr)
            print(f"Try specifying a different encoding with --encoding option, e.g., --encoding gbk for Chinese files.", file=sys.stderr)
            sys.exit(1)
        except IOError as e:
            print(f"Error reading file '{args.input_file}': {e}", file=sys.stderr)
            sys.exit(1)
    
    # Filter lines and apply replacement if specified
    total_lines = len(lines)
    filtered_lines = []
    
    # Parse replace argument if provided
    replace_old = None
    replace_new = None
    if args.replace:
        if ':' in args.replace:
            replace_old, replace_new = args.replace.split(':', 1)
        else:
            print(f"Error: Invalid replace format. Use 'old:new' format.", file=sys.stderr)
            sys.exit(1)
    
    # Simplified filtering logic
    for line in lines:
        matched = False
        if args.case_sensitive:
            # Case-sensitive search
            if args.search_content in line:
                matched = True
        else:
            # Case-insensitive search (default)
            if args.search_content.lower() in line.lower():
                matched = True
        
        if matched:
            # Apply replacement if specified
            if replace_old:
                if args.case_sensitive:
                    # Case-sensitive replacement
                    line = line.replace(replace_old, replace_new)
                else:
                    # Case-insensitive replacement
                    import re
                    line = re.sub(replace_old, replace_new, line, flags=re.IGNORECASE)
            filtered_lines.append(line)
    
    matched_lines = len(filtered_lines)
    
    # Determine output file path
    output_path = args.output
    if not output_path:
        # Generate default output file in the same directory as input file
        input_dir = os.path.dirname(os.path.abspath(args.input_file))
        input_basename = os.path.basename(args.input_file)
        input_name, input_ext = os.path.splitext(input_basename)
        
        # For non-text files (Excel, CSV), use .txt extension for output
        # For text files, keep the original extension
        output_ext = input_ext
        # List of non-text file extensions that should use .txt output
        non_text_extensions = ['.xlsx', '.xls', '.csv']
        if output_ext.lower() in non_text_extensions:
            output_ext = '.txt'
        
        # Generate output filename: inputname_filtered_searchcontent.ext
        # Replace spaces in search content with underscores for filename safety
        safe_search_content = args.search_content.replace(' ', '_')
        # Limit search content length in filename to avoid issues
        safe_search_content = safe_search_content[:50] if len(safe_search_content) > 50 else safe_search_content
        
        # Generate base filename
        base_filename = f"{input_name}_filtered_{safe_search_content}"
        output_filename = f"{base_filename}{output_ext}"
        output_path = os.path.join(input_dir, output_filename)
        
        # Ensure filename is unique: add counter if file already exists
        counter = 1
        while os.path.exists(output_path):
            output_filename = f"{base_filename}_{counter}{output_ext}"
            output_path = os.path.join(input_dir, output_filename)
            counter += 1
    
    # Write results to file
    try:
        with open(output_path, 'w', encoding=args.encoding) as f:
            f.writelines(filtered_lines)
        # Output simplified feedback
        print(f"\n保存位置: {os.path.abspath(output_path)}", file=sys.stderr)
        print(f"匹配行数: {matched_lines}", file=sys.stderr)
    except IOError as e:
        print(f"Error writing to file '{output_path}': {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == '__main__':
    main()
