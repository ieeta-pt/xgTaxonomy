import csv
from collections import defaultdict
from pathlib import Path
import os 


def average_nc_per_genomic_type(filepaths):
    all_results = {}
    for filepath in filepaths:
        # Create a defaultdict to store the values for each compressor and genomic type
        data = defaultdict(lambda: defaultdict(list))
        genomictype =''
        # Open the CSV file and read the values
        with open(filepath) as f:
            reader = csv.reader(f)
            headers = next(reader)
            headers.pop(0)
            for i, row in enumerate(reader):
                if i == 0:
                    genomictype = row[0]
                for j,value in enumerate(row):
                    if j!=0:
                        compressor = headers[j-1]
                        if value == '':
                            value = 100000
                        data[compressor][genomictype].append(float(value))
        # Compute the averages and store the results in a dictionary
        results = {}
        for compressor, groups in data.items():
            results[compressor] = {}
            for genomic_type, values in groups.items():
                avg = sum(values) / len(values)
                results[compressor][genomic_type] = avg
            # Compute average for that compressor globally
            results[compressor]['global'] = sum(results[compressor].values()) / len(results[compressor].values())
        all_results[filepath] = results
    return all_results

def save_results(all_results, folderpath):
    for filepath, results in all_results.items():
        file_name = Path(filepath).stem
        with open(f"{folderpath}/{file_name}_compression_results.txt", 'w') as f:
            f.write(f'Compression results for {file_name}:\n')
            for compressor, groups in results.items():
                f.write(f'Compressor: {compressor}\n')
                for genomic_type, avg in groups.items():
                    f.write(f'    {genomic_type}: {avg}\n')

def main():
    filepaths = ['../results/genome_features.csv', '../results/proteome_features.csv']
    all_results = average_nc_per_genomic_type(filepaths)
    save_results(all_results, '../results')


if __name__ == '__main__':
    # get the current working directory
    cwd = os.getcwd()
    # check if the current working directory ends with 'src'
    if cwd.endswith('src'):
        main()
    else:
        print("Error: script must be run from the 'src' folder")

    