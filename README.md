# xgTaxonomy

<H2><b>Cross-reference of Genomic Taxonomy</b></H2>

### About

xgTaxonomy is a new method for metagenomic classification that utilizes data compression algorithms, known as compressors, to classify genomic sequences. Our two-step evaluation process shows that this approach outperforms existing methods in terms of accuracy and reliability. Additionally, combining features from multiple compressors improves classification accuracy by 26,22%. This method offers a promising strategy for improving the accuracy and reliability of metagenomic classification and provides insights into the statistical and algorithmic nature of genomic data.

### Team

* Jorge M. Silva<sup id="a1">[1](#f1)</sup>
* João R. Almeida<sup id="a1">[1](#f1)</sup><sup id="a2">[2](#f2)</sup>

1. <small id="f1"> DETI/IEETA, LASI, University of Aveiro, Aveiro, Portugal </small> [↩](#a1)
2. <small id="f2"> University of A Coruña, A Coruña, Spain </small> [↩](#a2)

### Getting Started

#### Prerequisites

- Git
- Docker and Docker-compose (if using the Docker option)

#### Download Project

Get xgTaxonomy project using:

```bash
git clone https://github.com/bioinformatics-ua/xgTaxonomy.git
cd xgTaxonomy/
```

#### Using Docker

To perform installation correctly, docker and docker compose must be installed in the system (see <https://docs.docker.com/engine/install/ubuntu/>).

Then, follow these instructions:

```sh
git clone https://github.com/bioinformatics-ua/xgTaxonomy.git
cd xgTaxonomy
docker-compose build
docker-compose up -d && docker exec -it xgTaxonomy bash && docker-compose down
```

#### Install Compressors

Give run Install Compressors for Benchmark:

``` bash
bash install_compressors.sh;
```

### Result Replication

To run the pipeline and obtain all the Reports in the folder reports, use the following commands.

#### Download sequences I

For obtaining random sequences for baseline test performance run:

``` bash
cd src/
python3 getSampleSequences.py 
```

#### Baseline test

For baseline compression test run:

``` bash
cd src/
python3 compress_baseline.py
```

#### Download sequences II

For obtaining random sequences for taxonomic classification run:

``` bash
cd src/
python3 getDatabaseSequences.py 
```

### Classifiers

#### F1-score and accuracy for each compressor

```bash
cd src/
python3 classifier.py -b > ../results/f1score_accuracy_single.txt
```

#### Classification report for each compressor

```bash
cd src/
python3 classifier.py -cr > ../results/classification_reports_single.txt
```

#### Classification f1-score and accuracy for all genomic features

```bash
cd src/
python3 classifier.py -ag -b > ../results/f1_score_accuracy_all_genome_features.txt
python3 classifier.py -ag -cr > ../results/classification_report_all_genome_features.txt 
```

#### Classification f1-score and accuracy for all proteomic features

```bash
cd src/
python3 classifier.py -ap -b > ../results/f1_score_accuracy_all_proteome_features.txt
python3 classifier.py -ap -cr > ../results/classification_report_all_proteome_features.txt
```

#### Classification report using all compression features

```bash
cd src/
python3 classifier.py -cr -ac > ../results/classification_report_all_columns.txt
```

#### F1-score and accuracy using all compression features

```bash
cd src/
python3 classifier.py -ac -b > ../results/f1score_accuracy_all_columns.txt
```

#### Feature selection for f1-score and accuracy

```bash
cd src/
python3 classifier.py -fs -ac -b > ../results/feature_selection.txt
```

#### Classification f1-score and accuracy for all possible feature combinations (brute force)

```bash
cd src/
python3 classifier.py -bf -b > ../results/f1score_accuracy_all_combinations.txt
```

#### Classification report for all compressors (brute force)

```bash
cd src/
python3 classifier.py -bf -cr > ../results/classification_report_all_combinations.txt
```

#### Test Correlation

```bash
cd src/
python3 correlateTable.py
```

### Cite

Please cite the following, if you use xgTaxonomy in your work:

```bib
in progress
```

### Issues

Please let us know if there are any
[issues](https://github.com/bioinformatics-ua/COMPACT/issues).

### License
xgTaxonomy is released under the MIT License.