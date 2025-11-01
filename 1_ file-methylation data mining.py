def filter_rows(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            columns = line.split()
            if len(columns) >= 4 and columns[3] == '100':
                outfile.write(line)


input_file = r"E:\KH LAB\lab work\Methylation data mining\final data mining results\SRR3274243_normal_lung\Galaxy16-[MethylDackel on data 14_ CpG metylation levels].bedgraph"
output_file = r"E:\KH LAB\lab work\Methylation data mining\final data mining results\SRR3274243_normal_lung\methylated chr all.txt"
filter_rows(input_file, output_file)
