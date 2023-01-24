from os import listdir, path, makedirs, walk, system, remove
from os.path import isfile, join
import random
import shutil
import wget
import requests
from bs4 import BeautifulSoup
from config import summariesPath, dstPathOfOriginalSequences, dstPathOfOriginalAminoacids, tmpPrt, orfPath

#ensuring reproducibility
random.seed(0)

numberOfEntries = {
    "viral" : 10,
    "bacteria" : 10,
    "archaea" : 10,
    "fungi" : 10,
    "plant" : 10,
    "protozoa" : 10
}
multiColumnPos = 19


def main():
    selectedSequences = getSequences()
    _initialize()
    downloadSequences(selectedSequences["from_url"], dstPathOfOriginalSequences)
    splitSequences(selectedSequences, dstPathOfOriginalSequences)
    to_amino_seq()


def _initialize():
    if not path.exists(dstPathOfOriginalSequences):
        makedirs(dstPathOfOriginalSequences)
    if not path.exists(tmpPrt):
        makedirs(tmpPrt)
    if not path.exists(dstPathOfOriginalAminoacids):
        makedirs(dstPathOfOriginalAminoacids)


def splitSequences(selectedSequences, dstPathOfOriginalSequences):
    allFiles = [f for f in listdir(dstPathOfOriginalSequences) if isfile(join(dstPathOfOriginalSequences, f))]
    for x in selectedSequences:
        for domain in selectedSequences[x]:
            dst = join(dstPathOfOriginalSequences, domain)
            if not path.exists(dst):
                makedirs(dst)
            for entry in selectedSequences[x][domain]:
                fileName = entry
                if "/" in entry:
                    fileName = entry.split("/")[-1]
                fileToMove = list(filter(lambda x: fileName in x, allFiles))[0]
                shutil.move(join(dstPathOfOriginalSequences, fileToMove), dst)
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
            listOfValues = sorted(random.sample(range(0, lenOfFile), numberOfEntries[domain]))
            selectedSequences["from_url"][domain] = _getRandomEntries(listOfValues, content)
    return selectedSequences


def _getRandomEntries(listOfValues, content, singleColumn=False):
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
    for (root, _, file) in walk(dstPathOfOriginalSequences):
        for f in file:
            load_path = join(root,f)
            shutil.copy(load_path, tmpPrt)
            decompressedPath = join(tmpPrt, f).replace(".gz", "")
            system(f'gzip -d {join(tmpPrt, f)}')
            system(f'{orfPath}/orfm {decompressedPath} > {tmpPrt}/PROTEIN')
            domain = root.split('/')[-1]
            savepath = join(dstPathOfOriginalAminoacids, domain)
            if not path.exists(savepath):
                makedirs(savepath)
            save_name = join(savepath,f.replace(".gz", "").replace(".fna", ".amino"))
            system(f' cat {tmpPrt}/PROTEIN | grep -v ">" | tr -d -c "ACDEFGHIKLMNPQRSTVWY" > {save_name}')
            system(f'gzip {save_name}')
            remove(decompressedPath)
   


if __name__ == "__main__":
    main()

    