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
                genomictype = row[0]
                for j,value in enumerate(row):
                    if j!=0:
                        compressor = headers[j-1].replace('_comp','')
                        try:
                            result = float(value)
                        except ValueError:
                            result = 0                        
                        data[compressor][genomictype].append(result)
        
        # Compute the averages and store the results in a dictionary
        results = {}
        for compressor, groups in data.items():
            results[compressor] = {}
            for genomic_type, values in groups.items():
                avg = sum(values) / len(values)
                results[compressor][genomic_type] = avg        
        all_results[filepath] = results
    return all_results

def all_values_non_zero(d):
    for value in d.values():
        if value == 0:
            return False
    return True

def save_results(all_results, folderpath):
    for filepath, results in all_results.items():
        file_name = Path(filepath).stem
        with open(f"{folderpath}/{file_name}_compression_results.csv", mode='w') as f:
            fieldnames = ['Taxonomic Group'] + list(results['blzpack'].keys())
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            for genomic_type, values in results.items():
                row = {'Taxonomic Group': genomic_type}
                for compressor, avg in values.items():
                    row[compressor] = avg
                if (all_values_non_zero(row)):
                    writer.writerow(row)
            
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

    