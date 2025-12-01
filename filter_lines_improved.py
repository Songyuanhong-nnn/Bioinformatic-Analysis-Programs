#!/usr/bin/env python3
import argparse
import sys

def main():
    parser = argparse.ArgumentParser(
        description='Filter lines in a file that contain a specific string',
        epilog='Example usage:\n  %(prog)s input.txt "search term"\n  %(prog)s input.txt "test" -i\n  %(prog)s input.txt "keyword" -o output.txt',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('input_file', help='Path to the input text file')
    parser.add_argument('search_content', help='String to search for in each line')
    parser.add_argument('-o', '--output', help='Path to write filtered results (default: print to console)')
    parser.add_argument('-i', '--ignore-case', action='store_true', help='Perform case-insensitive search')
    
    args = parser.parse_args()
    
    try:
        with open(args.input_file, 'r', encoding='utf-8') as f:
            lines = f.readlines()
    except FileNotFoundError:
        print(f"Error: File '{args.input_file}' not found.", file=sys.stderr)
        sys.exit(1)
    except IOError as e:
        print(f"Error reading file '{args.input_file}': {e}", file=sys.stderr)
        sys.exit(1)
    
    filtered_lines = []
    for line in lines:
        if args.ignore_case:
            if args.search_content.lower() in line.lower():
                filtered_lines.append(line)
        else:
            if args.search_content in line:
                filtered_lines.append(line)
    
    if args.output:
        try:
            with open(args.output, 'w', encoding='utf-8') as f:
                f.writelines(filtered_lines)
        except IOError as e:
            print(f"Error writing to file '{args.output}': {e}", file=sys.stderr)
            sys.exit(1)
    else:
        sys.stdout.writelines(filtered_lines)

if __name__ == '__main__':
    main()
