#!/usr/bin/env python3
"""
Split a model grid file into two separate files:
1. Grid geometry file (all columns except Val)
2. Values file (only Val column(s))
"""

import argparse
from pathlib import Path


def split_grid_file(input_file, grid_output=None, values_output=None):
    """
    Split a model grid file into geometry and values files.
    
    Input format:
        Line 1: Number of cells
        Lines 2+: X1 X2 Y1 Y2 Z1 Z2 Val i j k (scalar, 10 columns)
              or: X1 X2 Y1 Y2 Z1 Z2 Val1 Val2 Val3 i j k (vector, 12 columns)
    
    Output grid file: X1 X2 Y1 Y2 Z1 Z2 i j k
    Output values file: Val (scalar) or Val1 Val2 Val3 (vector)
    """
    input_path = Path(input_file)
    
    if grid_output is None:
        grid_output = input_path.stem + "-grid" + input_path.suffix
    if values_output is None:
        values_output = input_path.stem + "-values" + input_path.suffix
    
    with open(input_file, 'r') as f_in:
        # Read header
        header = f_in.readline()
        n_cells = int(header.strip())
        
        # Peek at first data line to detect format
        first_line = f_in.readline()
        parts = first_line.split()
        n_cols = len(parts)
        
        if n_cols == 10:
            n_values = 1
        elif n_cols == 12:
            n_values = 3
        else:
            raise ValueError(
                f"Unexpected number of columns: {n_cols}. Expected 10 (scalar) or 12 (vector)."
            )
        
        print(f"Detected {'vector (3 components)' if n_values == 3 else 'scalar'} data")
    
    # Reopen and process
    with open(input_file, 'r') as f_in, \
         open(grid_output, 'w') as f_grid, \
         open(values_output, 'w') as f_val:
        
        # Write headers
        f_in.readline()  # skip header
        f_grid.write(header)
        f_val.write(header)
        
        # Process each cell line
        for line_num, line in enumerate(f_in, start=2):
            parts = line.split()
            
            if len(parts) != n_cols:
                raise ValueError(
                    f"Line {line_num}: expected {n_cols} columns, got {len(parts)}"
                )
            
            # X1 X2 Y1 Y2 Z1 Z2 are always first 6
            coords = parts[0:6]
            
            if n_values == 1:
                # Scalar: X1 X2 Y1 Y2 Z1 Z2 Val i j k
                val = parts[6]
                i, j, k = parts[7:10]
                f_val.write(f"{val}\n")
            else:
                # Vector: X1 X2 Y1 Y2 Z1 Z6 Val1 Val2 Val3 i j k
                val1, val2, val3 = parts[6:9]
                i, j, k = parts[9:12]
                f_val.write(f"{val1} {val2} {val3}\n")
            
            # Write grid: X1 X2 Y1 Y2 Z1 Z2 i j k
            f_grid.write(f"{' '.join(coords)} {i} {j} {k}\n")
    
    print(f"Created: {grid_output}")
    print(f"Created: {values_output}")
    print(f"Processed {n_cells} cells")


def main():
    parser = argparse.ArgumentParser(
        description="Split model grid file into grid and values files"
    )
    parser.add_argument("input_file", help="Input grid file")
    parser.add_argument("-g", "--grid", help="Output grid file name")
    parser.add_argument("-v", "--values", help="Output values file name")
    
    args = parser.parse_args()
    split_grid_file(args.input_file, args.grid, args.values)


if __name__ == "__main__":
    main()
