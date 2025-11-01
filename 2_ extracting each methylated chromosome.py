import csv

input_file_path = r"E:\KH LAB\lab work\Methylation data mining\final data mining results\SRR3269805_normal_liver\methylated chr all.txt"
output_file_path = r"E:\KH LAB\lab work\Methylation data mining\final data mining results\SRR3269805_normal_liver\chrY\methylated chr Y.txt"

# Define the delimiter used in the input file
delimiter = '\t'


try:
    with open(input_file_path, 'r') as infile:
        # Initialize the CSV reader
        reader = csv.reader(infile, delimiter=delimiter)

        # Prepare to write to the output file
        with open(output_file_path, 'w', newline='') as outfile:
            # Initialize the CSV writer
            writer = csv.writer(outfile, delimiter=delimiter)

            # Iterate over the rows in the input file
            for row in reader:
                # Check if the first column contains 'chrY'
                if len(row) > 0 and row[
                    0] == 'chrY':  # Ensure there is at least one column and the first column is 'chrY'
                    # Write the entire row to the new file
                    writer.writerow(row)

    print(f"Extracted rows written to {output_file_path}")

except FileNotFoundError:
    print(f"Error: The file '{input_file_path}' does not exist.")

except IOError:
    print(f"Error: Could not read or write to the file.")
