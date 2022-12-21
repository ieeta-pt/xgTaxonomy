#!/usr/bin/env python

import warnings
from xgboost import XGBClassifier 
from sklearn.metrics import accuracy_score, classification_report, f1_score 
from sklearn.model_selection import train_test_split
from sklearn.svm import LinearSVC
from sklearn.feature_selection import SelectFromModel
from sklearn.ensemble import ExtraTreesClassifier
from itertools import combinations
import csv
import sys
import numpy as np
from config import genomeFeaturesFilePath, proteomeFeaturesFilePath, numIterations
import argparse

compressor = {
    1:"blzpack_g",
    2:"bsc_g",
    3:"bzip2_g", 
    4:"GeCo3_g",
    5:"gzip_g",
    6:"JARVIS_g",
    7:"lizard_g",
    8:"lz4_g",
    9:"lzop_g",
    10:"mbgc_g",
    11:"MFCompress_g",
    12:"naf_g",
    13:"NUHT_g",
    14:"snzip_g",
    15:"zip_g",
    16:"xz_g",
    17:"zstd_g",
    18:"blzpack_p",
    19:"bsc_p",
    20:"bzip2_p",
    21:"gzip_p",
    22:"lizard_p",
    23:"lz4_p",
    24:"lzop_p", 
    25:"snzip_p", 
    26:"zip_p",
    27:"xz_p",
    28:"zstd_p" 
}





def warn(*args, **kwargs):
    pass
warnings.warn = warn

def concatenate_csv(args):
    # Open both CSV files in read mode
    with open(args.genome_filename, 'r') as file1, open(args.proteome_filename, 'r') as file2:
        # Read the contents of both files into separate variables
        reader1 = csv.reader(file1)
        reader2 = csv.reader(file2)
        combined_rows = []
        # Iterate over the rows of both lists
        for row1, row2 in zip(reader1, reader2):
            # Check if the first elements of the rows are the same
            if row1[0] == row2[0]:
                # Add the contents of the row to the combined rows list, removing the first element from row2
                combined_rows.append(row1 + row2[1:])
            else:
                break
        return combined_rows
       

def flatten_columns(args, columns):
    genome_columns, proteome_columns = columns[0], columns[1:]
    # Get the number of columns in the genome and proteome files
    with open(args.genome_filename, 'r') as file:
        reader = csv.reader(file)
        num_genome_columns = len(next(reader))

    # Flatten the columns
    return genome_columns + [i + num_genome_columns for i in proteome_columns]


def ReadData(args, columns):
    domains = {
        "viral": 0, 
        "bacteria": 1,
        "archaea": 2, 
        "fungi": 3,
        "protozoa": 4,
    }

    X_test, y_test = [], []

    # Call the concatenate_csv function and store the returned list of rows
    combined_rows = concatenate_csv(args)

    # Discard the first row (header)
    combined_rows = combined_rows[1:]

    # Select the desired columns from the combined rows list and convert them to floats
    X_test = [[float(row[i]) for i in columns] for row in combined_rows]

    # Get the labels from the first column of the combined rows
    y_test = [domains[row[0]] for row in combined_rows]
    
    return np.array(X_test), np.array(y_test)




def Classify(args, columns):
    domains = ["Viral", "Bacteria", "Archaea", "Fungi", "Protozoa"]
    accuracy_XGB = []
    f1score_XGB = []

    # Flatten the columns list
    columns = flatten_columns(args,columns)

    data, labels = ReadData(args, columns)

    if args.features_selection:
        clf = ExtraTreesClassifier(n_estimators=50)
        clf = clf.fit(data, labels)
        model = SelectFromModel(clf, prefit=True)
        print(data.shape)
        data = model.transform(data)
        print(data.shape)

    for a in range(numIterations):
        X_train, X_test, y_train, y_test = train_test_split(data, labels, test_size=0.20, stratify=labels, random_state=a)
        model = XGBClassifier(max_depth=12, learning_rate=0.89, n_estimators=500, eval_metric='mlogloss')
        model.fit(X_train, y_train)
        y_pred = model.predict(X_test)
        predictions = [round(value) for value in y_pred]
        accuracy_XGB.append(accuracy_score(y_test, predictions))
        f1score_XGB.append(f1_score(y_test, y_pred, average='weighted'))
        if args.classification_report:
            print("Classification report, using:", end=" ")
            print(" ".join([compressor[x] for x in columns]))
            print(classification_report(y_test, y_pred, target_names=domains, digits=4))
            break

    if args.accuracy:
        print("Accuracy of XGB, using:", end=" ")
        print(" ".join([compressor[x] for x in columns]))
        print(sum(accuracy_XGB)/len(accuracy_XGB))

    elif args.f1_score:
        print("F1 score of XGB, using:", end=" ")
        print(" ".join([compressor[x] for x in columns]))
        print(sum(f1score_XGB)/len(f1score_XGB))
    
    elif args.both:
        print("Accuracy of XGB, using:", end=" ")
        print(" ".join([compressor[x] for x in columns]))
        print(sum(accuracy_XGB)/len(accuracy_XGB))
        print("F1 score of XGB, using:", end=" ")
        print(" ".join([compressor[x] for x in columns]))
        print(sum(f1score_XGB)/len(f1score_XGB))
    print()
        

def help(show=False):
    parser = argparse.ArgumentParser(description="")
    helper = parser.add_argument_group('System settings', 'System parameters to run the classifier in the different modes')
    helper.add_argument('-g', '--genomeFilename', dest='genome_filename', \
                        type=str, default=genomeFeaturesFilePath, \
                        help=f'The system settings file (default: {genomeFeaturesFilePath})')
    helper.add_argument('-p', '--proteomeFilename', dest='proteome_filename', \
                        type=str, default=proteomeFeaturesFilePath, \
                        help=f'The system settings file (default: {proteomeFeaturesFilePath})')     
    helper.add_argument('-f1', '--f1-score', default=False, action='store_true', \
                            help='This flag produces the classificarion report using the F1-score (default: False)')
    helper.add_argument('-a', '--accuracy', default=False, action='store_true', \
                            help='This flag produces the classificarion report using the Accuracy metric (default: False)')
    helper.add_argument('-b', '--both', default=False, action='store_true', \
                            help='This flag produces the classificarion report using the both metrics (default: False)') 
    helper.add_argument('-fs', '--features-selection', default=False, action='store_true', \
                            help='This flag performs feature selection (default: False)') 
    helper.add_argument('-ac', '--all-columns', default=False, action='store_true', \
                            help='This flag classifies using all features (default: False)') 
    helper.add_argument('-cr', '--classification-report', default=False, action='store_true', \
                            help='This flag generates the classification report (default: False)') 
    helper.add_argument('-bf', '--brute-force', default=False, action='store_true', \
                            help='This flag performs brute force classification of all possible combination of features (default: False)') 
    if show:
        parser.print_help()
    return parser.parse_args()
    

if __name__ == "__main__":
    
    args = help()
    if args.accuracy or args.f1_score or args.both or args.classification_report:
        if args.all_columns:
            Classify(args, [[0,1,2,3,4],[0,1,2,3,4]])
        
        elif args.brute_force:
            print(2)
            all_comb_list=[]
            for x in range(1,29,1):
                com_list = list(combinations(range(1,29), x+1))
                [Classify(args,list(ele)) for ele in com_list]
        else:
            for column in range(1,29):
                Classify(args, [[column]])
    else:
        help(True)