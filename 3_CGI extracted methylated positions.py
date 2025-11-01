file1 = r"E:\KH LAB\lab work\data_collection\cgi_files_hg38\cgi_start_end\chrY.txt" #filepath of file containing range of CGI
file2 = r"E:\KH LAB\lab work\Methylation data mining\final data mining results\SRR3269805_normal_liver\chrY\methylated chr Y.txt" #filpeath of file containing positions of all methylated CpG in chrY
output_file = r"E:\KH LAB\lab work\Methylation data mining\final data mining results\SRR3269805_normal_liver\chrY\chrY_methcgi_strt_end.txt"

# Read the ranges from file1
ranges = []
with open(file1, 'r') as f1:
    next(f1)
    for line in f1:
        parts = line.strip().split('\t')
        chrom = parts[1]
        start = int(parts[2])
        end = int(parts[3])
        ranges.append((chrom, start, end))

# Read the positions from file2 and filter them
filtered_positions = []
with open(file2, 'r') as f2:
    next(f2)
    range_index = 0
    current_range = None
    for line in f2:
        parts = line.strip().split('\t')
        #chrom = parts[0]
        start = int(parts[1])
        end = int(parts[2])


        while range_index < len(ranges) and (
                ranges[range_index][0] < chrom or (ranges[range_index][0] == chrom and ranges[range_index][2] < start)):
            range_index += 1


        if range_index == len(ranges):
            break


        if ranges[range_index][0] == chrom and ranges[range_index][1] <= start <= ranges[range_index][2]:
            if current_range != ranges[range_index]:
                if current_range is not None:
                    filtered_positions.append("")  # Add a blank line before starting a new range
                current_range = ranges[range_index]
            filtered_positions.append(line.strip())


with open(output_file, 'w') as out:
    out.write('track type="bedGraph" description="filtered CpG methylation fractions"\n')
    for position in filtered_positions:
        out.write(position + '\n')

print(f"Filtered positions saved to {output_file}")
