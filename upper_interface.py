#!/usr/bin/python3

import os

def round_x_value(x):
    return round(float(x), 3)

def main():
    input_file = "interface_fft.txt"  # Change this to your input file name
    output_file = "output.txt"  # Change this to your output file name

    with open(input_file, 'r') as infile:
        with open(output_file, 'w') as outfile:
            for line in infile:
                if line.strip():  # Check if the line is not empty
                    x, y, z = line.split()
                    rounded_x = round_x_value(x)
                    outfile.write(f"{rounded_x} {y} {z}\n")


    with open(output_file, "r") as file:
        data = [line.strip().split() for line in file]
        sorted_data = sorted(data, key=lambda x: float(x[0]))  # Sort based on the first column (x values)

    sorted_output_file = "sorted_output.txt"
    with open(sorted_output_file, "w") as output_file:
        for row in sorted_data:
            output_file.write('\t'.join(row) + '\n')  # Write each row to the file, separated by tabs

    with open(sorted_output_file, "r") as file:
        sorted_data = [line.strip().split() for line in file]

    max_y_values = {}

    for row in sorted_data:
        x = row[0]
        y = row[1]

        if x not in max_y_values or float(y) > float(max_y_values[x][1]):
            max_y_values[x] = row

    max_y_values_file = "max_y_values.txt"
    with open(max_y_values_file, "w") as output_file:
        for x, row in max_y_values.items():
            output_file.write('\t'.join(row) + '\n')  

if __name__ == "__main__":
    main()


os.system(f"rm -r output.txt")
os.system(f"rm -r sorted_output.txt")
os.system(f"mv max_y_values.txt output.txt")

