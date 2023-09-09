import os
from Bio import SeqIO

def get_dataset_summary(dataset_directory):
    domains = ['viral', 'bacteria', 'archaea', 'fungi', 'plant', 'protozoa']
    summary = {}

    for domain in domains:
        domain_path = os.path.join(dataset_directory, domain)
        
        # Ensure the domain directory exists
        if os.path.exists(domain_path):
            sequence_lengths = []
            number_of_sequences = 0

            # Traverse each sequence file within the domain directory
            for seq_file in os.listdir(domain_path):
                file_path = os.path.join(domain_path, seq_file)
                
                # Use Biopython's SeqIO to parse sequence files
                for record in SeqIO.parse(file_path, "fasta"):
                    sequence_lengths.append(len(record.seq))
                    number_of_sequences += 1

            # Store the results in the summary dictionary
            summary[domain] = {
                'average_length': sum(sequence_lengths) / len(sequence_lengths) if sequence_lengths else 0,
                'min_length': min(sequence_lengths) if sequence_lengths else 0,
                'max_length': max(sequence_lengths) if sequence_lengths else 0,
                'number_of_sequences': number_of_sequences
            }

    return summary

if __name__ == "__main__":
    dataset_directory = 'path_to_your_dataset'  # Replace with your dataset directory path
    summary = get_dataset_summary(dataset_directory)
    for domain, stats in summary.items():
        print(f"Domain: {domain}")
        print(f"  Average Length: {stats['average_length']:.2f}")
        print(f"  Min Length: {stats['min_length']}")
        print(f"  Max Length: {stats['max_length']}")
        print(f"  Number of Sequences: {stats['number_of_sequences']}\n")
