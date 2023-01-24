from os import listdir, path, makedirs, walk, system, remove
from os.path import isfile, join
import random
from pprint import pprint
import sys
import shutil
import wget
import requests
from bs4 import BeautifulSoup
from config import summariesPath, dstPathOfDatabaseSequences, locationOfDatabases, dstPathOfDatabaseAminoacids, tmpPrt, orfPath

#ensuring reproducibility
random.seed(0)

numberOfEntries = {
    "viral" : 1500,
    "bacteria" : 1500,
    "archaea" : 1210,
    "fungi" : 429,
    "plant" : 144,
    "protozoa" : 95
}
multiColumnPos = 19


def main():
    selectedSequences = getSequences()
    _initialize()
    downloadSequences(selectedSequences["from_url"], dstPathOfDatabaseSequences)
    splitSequences(selectedSequences, dstPathOfDatabaseSequences)
    to_amino_seq()


def _initialize():
    if not path.exists(dstPathOfDatabaseSequences):
        makedirs(dstPathOfDatabaseSequences)
    if not path.exists(tmpPrt):
        makedirs(tmpPrt)
    if not path.exists(dstPathOfDatabaseAminoacids):
        makedirs(dstPathOfDatabaseAminoacids)



def splitSequences(selectedSequences, dstPathOfDatabaseSequences):
    allFiles = [f for f in listdir(dstPathOfDatabaseSequences) if isfile(join(dstPathOfDatabaseSequences, f))]
    for x in selectedSequences:
        for domain in selectedSequences[x]:
            dst = join(dstPathOfDatabaseSequences, domain)
            if not path.exists(dst):
                makedirs(dst)
            for entry in selectedSequences[x][domain]:
                fileName = entry
                if "/" in entry:
                    fileName = entry.split("/")[-1]
                fileToMove = list(filter(lambda x: fileName in x, allFiles))[0]
                shutil.move(join(dstPathOfDatabaseSequences, fileToMove), dst)
    print("Done!")


def getSequences():
    selectedSequences = {
        "from_url":{},
        "from_db":{}
    }
    onlyfiles = [f for f in listdir(summariesPath) if isfile(join(summariesPath, f)) and join(summariesPath, f).endswith("summary.txt")]
    for fileName in onlyfiles:
        with open(join(summariesPath, fileName), 'r') as fp:
            domain = fileName.split("_")[0]
            content = fp.readlines()
            lenOfFile = len(content)
            if lenOfFile>numberOfEntries[domain]:
                listOfValues = sorted(random.sample(range(0, lenOfFile), numberOfEntries[domain]))
            else:
                listOfValues =  sorted(range(0, lenOfFile))
            selectedSequences["from_url"][domain] = _getEntries(listOfValues, content)
    return selectedSequences


def _getEntries(listOfValues, content, singleColumn=False):
    readsOfInterest = []
    for n in listOfValues:
        if singleColumn:
            readsOfInterest.append(content[n].replace("\n", ""))
        else:
            readsOfInterest.append(content[n].split("\t")[multiColumnPos])
    return readsOfInterest


def downloadSequences(selectedSequences, dst):
    for domain in selectedSequences:
        for entry in selectedSequences[domain]:
            r = requests.get(entry)
            soup = BeautifulSoup(r.text, 'html.parser')
            for link in soup.find_all('a'):
                if "genomic.fna.gz" in link.get('href') and "from_genomic" not in link.get('href'):
                    fLink = join(entry, link.get('href'))
                    print(f"Downloading {fLink} to {dst}")
                    wget.download(fLink, out=str(dst))

def to_amino_seq():
    for (root, _, file) in walk(dstPathOfDatabaseSequences):
        for f in file:
            load_path = join(root,f)
            shutil.copy(load_path, tmpPrt)
            decompressedPath = join(tmpPrt, f).replace(".gz", "")
            system(f'gzip -d {join(tmpPrt, f)}')
            system(f'{orfPath}/orfm {decompressedPath} > {tmpPrt}/PROTEIN')
            domain = root.split('/')[-1]
            savepath = join(dstPathOfDatabaseAminoacids, domain)
            if not path.exists(savepath):
                makedirs(savepath)
            save_name = join(savepath,f.replace(".gz", "").replace(".fna", ".amino"))
            system(f' cat {tmpPrt}/PROTEIN | grep -v ">" | tr -d -c "ACDEFGHIKLMNPQRSTVWY" > {save_name}')
            system(f'gzip {save_name}')
            remove(decompressedPath)
   

if __name__ == "__main__":
    main()

    