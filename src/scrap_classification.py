import re
import csv
import os

def parse_file(filepath):
    # Create a dictionary to store the results
    results = {}
    with open(filepath) as f:
        lines = f.readlines()
        for i, line in enumerate(lines):
            # Search for the line that starts with "Classification report, using:"
            if line.startswith("Classification report, using:"):
                # Extract the method name from the line
                method = line.split(":")[1].strip()
                results[method] = {}
                for j in range(i+3, i+8):
                    # Split the line by whitespace and extract the class name, precision, recall and f1-score
                    parts = re.split(r"\s+", lines[j].strip())
                    class_name = parts[0]
                    f1_score = float(parts[3])
                    results[method][class_name] = f1_score
    return results

def save_results(results, filepath_prefix):
    # Create a new dictionary to store the results for genome and proteome separately
    results_genome = {}
    results_proteome = {}
    for method, values in results.items():
        if method.endswith("_g"):
            results_genome[method] = values
        elif method.endswith("_p"):
            results_proteome[method] = values
            
    # Save the results for genome
    with open(f"{filepath_prefix}_genome.csv", mode='w') as f:
        fieldnames = ['Taxonomic Group'] + list(results_genome['blzpack_g'].keys())
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for class_name, values in results_genome.items():
            row = {'Taxonomic Group': class_name}
            for method, f1_score in values.items():
                row[method] = f1_score
            writer.writerow(row)
            
    # Save the results for proteome
    with open(f"{filepath_prefix}_proteome.csv", mode='w') as f:
        fieldnames = ['Taxonomic Group'] + list(results_proteome['blzpack_p'].keys())
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for class_name, values in results_proteome.items():
            row = {'Taxonomic Group': class_name}
            for method, f1_score in values.items():
                row[method] = f1_score
            writer.writerow(row)

def main():
    filepath = '../results/classification_report_each_feature.txt'
    results = parse_file(filepath)
    save_results(results, '../results/classification_results_table')

if __name__ == '__main__':
    # get the current working directory
    cwd = os.getcwd()
    # check if the current working directory ends with 'src'
    if cwd.endswith('src'):
        main()
    else:
        print("Error: script must be run from the 'src' folder")
