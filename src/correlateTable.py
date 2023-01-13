import csv
from scipy.stats import pearsonr, spearmanr

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
    print(f'Pearson correlation for entire table: {corr}')
    corr, _ = spearmanr(all_data1, all_data2)
    print(f'Spearman correlation for entire table: {corr}')
    
    # Compute correlation for each genomic type
    for key in data1.keys():
        corr, _ = pearsonr(data1[key], data2[key])
        if abs(corr) < 0.5:
            print(f'No notable Pearson correlation for {key}!')
        else:
            print(f'Pearson correlation for {key}: {corr}')
        corr, _ = spearmanr(data1[key], data2[key])
        if abs(corr) < 0.5:
            print(f'No notable Spearman correlation for {key}!')
        else:
            print(f'Spearman correlation for {key}: {corr}')

def main():
    data1 = read_csv('../results/classification_results_table_proteome.csv')
    data2 = read_csv('../results/proteome_features_compression_results.csv')
    compute_correlation(data1, data2)

if __name__ == '__main__':
    main()
