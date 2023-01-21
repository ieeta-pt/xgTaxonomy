import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os

def analyze_csv(file_path:str):
    # Read csv file and store it in a DataFrame
    data = pd.read_csv(file_path)

    # Create a list of compressor columns
    compressor_columns = data.columns[data.columns.str.contains('comp')]

    # Create a dictionary to store the average compression time and standard deviation for each compressor
    compressor_stats = {}

    # Iterate through each compressor column
    for compressor in compressor_columns:
        compressor_data = data[compressor]
        avg_compression = compressor_data.mean()
        std_dev = compressor_data.std()
        compressor_stats[compressor] = (avg_compression, std_dev)

    # Create a list of time columns
    time_columns = data.columns[data.columns.str.contains('time')]

    # Create a dictionary to store the average time taken for each compressor
    time_stats = {}

    # Iterate through each time column
    for time_col in time_columns:
        time_data = data[time_col]
        avg_time = time_data.mean()
        time_stats[time_col] = avg_time
    return compressor_stats, time_stats



def plot_results(compressor_stats, time_stats):
    stats_df = pd.DataFrame(columns=['Compressor', 'NC', 'std_dev', 'Time'])
    for compressor, stats in compressor_stats.items():
        avg_compression, std_dev = stats
        compressor_name=compressor.replace("_comp","")
        avg_time = time_stats[compressor_name + "_time"]
        stats_df = stats_df.append({'compressor':compressor_name,'avg_compression':avg_compression,'std_dev':std_dev,'avg_time':avg_time}, ignore_index=True)
    stats_df.sort_values('avg_compression', inplace=True, ascending=False)
    stats_df['compressor'] = stats_df['compressor'].astype('category')
    plt.figure(figsize=(15,5))
    ax = sns.barplot(x="compressor", y="avg_compression", xerr=stats_df['std_dev'], data=stats_df)
    ax.set_xlabel("Compressor", fontweight='bold')
    ax.set_ylabel("Avg Compression", fontweight='bold')
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, fontweight='bold')
    for tick in ax.get_xticklabels():
        tick.set_fontsize(8)
    ax2 = ax.twinx()
    ax2.set_ylabel("Avg Time (s)", fontweight='bold')
    sns.lineplot(x="compressor", y="avg_time", data=stats_df, color="blue", ax=ax2)
    ax2.set_yscale('log')
    for tick in ax2.get_yticklabels():
        tick.set_fontsize(8)
    plt.savefig("../results/compression_time_graph.pdf")
    print(stats_df.sort_values("compressor"))
    plt.show()






def main():
    file_path = "../results/baseline_incomplete.csv"
    compressor_stats, time_stats = analyze_csv(file_path)
    print("compressor_stats")
    print(compressor_stats)
    print("time_stats")
    print(time_stats)
    plot_results(compressor_stats, time_stats)


if __name__ == '__main__':
    # get the current working directory
    cwd = os.getcwd()
    # check if the current working directory ends with 'src'
    if cwd.endswith('src'):
        main()
    else:
        print("Error: script must be run from the 'src' folder") 
