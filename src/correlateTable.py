import csv
from scipy.stats import pearsonr, spearmanr
import subprocess

def read_csv(filepath):
    data = {}
    with open(filepath) as f:
        reader = csv.DictReader(f)
        for row in reader:
            for key, value in row.items():
                key = key.lower()
                if key == 'taxonomic group':
                    continue
                
                if key not in data:
                    data[key] = []
                data[key].append(float(value))
    return data

def compute_correlation(data1, data2):
    # Compute correlation for entire table
    all_data1 = []
    all_data2 = []
    for key in data1.keys():
        all_data1 += data1[key]
        all_data2 += data2[key]
    corr, _ = pearsonr(all_data1, all_data2)
    if abs(corr) < 0.5:
            print(f'No notable Pearson correlation in table!')
    else:
        print(f'Pearson correlation for entire table: {corr}')
    corr, _ = spearmanr(all_data1, all_data2)
    if abs(corr) < 0.5:
        print(f'No notable Spearman correlation in table!')
    else:
        print(f'Spearman correlation for entire table: {corr}')
    

def main():
    subprocess.run(["python", "evaluateCompression.py"])
    subprocess.run(["python", "evaluateClassification.py"])
    print("Analysis of the Genomics")
    data1 = read_csv('../results/classification_results_table_genome.csv')
    data2 = read_csv('../results/genome_features_compression_results.csv')
    compute_correlation(data1, data2)
    print("Analysis of the Proteomics")
    data1 = read_csv('../results/classification_results_table_proteome.csv')
    data2 = read_csv('../results/proteome_features_compression_results.csv')
    compute_correlation(data1, data2)

if __name__ == '__main__':
    main()
