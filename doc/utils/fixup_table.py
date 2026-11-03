#!/usr/bin/env python3
"""
Reformat reStructuredText tables to ensure consistent column widths and correct syntax.
"""

import sys
import argparse

def reformat_rst_table(file_path):
    with open(file_path, "r") as f:
        lines = f.readlines()

    new_lines = []
    i = 0
    while i < len(lines):
        line = lines[i]
        if line.strip().startswith("+--"):
            # Start of a table, collect all table lines
            table_lines = []
            while i < len(lines) and (lines[i].strip().startswith("+--") or "|" in lines[i]):
                table_lines.append(lines[i])
                i += 1
            
            # Process table
            rows = []
            for t_line in table_lines:
                if "|" in t_line:
                    cells = [cell.strip() for cell in t_line.split("|") if cell.strip()]
                    if cells:
                        rows.append(cells)
            
            if not rows:
                new_lines.extend(table_lines)
                continue

            num_cols = max(len(row) for row in rows)
            col_widths = [0] * num_cols
            for row in rows:
                for j, cell in enumerate(row):
                    if j < num_cols:
                        col_widths[j] = max(col_widths[j], len(cell))

            separator = "+" + "+".join(["-" * (w + 2) for w in col_widths]) + "+"
            new_lines.append(separator + "\n")
            for row in rows:
                row_str = "|" + "|".join([f" {cell:<{col_widths[j]}} " for j, cell in enumerate(row)]) + "|"
                new_lines.append(row_str + "\n")
                new_lines.append(separator + "\n")
        else:
            new_lines.append(line)
            i += 1

    with open(file_path, "w") as f:
        f.writelines(new_lines)

def main():
    parser = argparse.ArgumentParser(description="Reformat RST tables.")
    parser.add_argument("files", nargs="+", help="RST files to reformat")
    args = parser.parse_args()

    for file_path in args.files:
        print(f"Reformatting table in {file_path}...")
        reformat_rst_table(file_path)

if __name__ == "__main__":
    main()
