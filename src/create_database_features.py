#!python
#!/usr/bin/env python

import os
import sys
from os import listdir, path, makedirs
from os.path import isfile, join, getsize
sys.path.append(os.path.abspath('../src'))

compressors_path = "../compressors/"
sequencesDir = "./database_sequences/"
mergedPath="../aux/database_sequences/"
mergedPathProteins="../aux/database_protein/"
root = "../baseline_features/"

Name = "genome_features"
proName = "proteome_features"
tmpDir = os.path.join(root, "tmp/")
tmpSeqPath = join("..","tmpSq")
UHT_COMP = "../compressors_files/Unbalanced-Huffman-Tree/dist/UHT_compress/"


def main():
    _initialize()
    csvContent = getCompressionValues()
    getProteinCompressionValues()

def writeCSVLine(line, filename):
    fileName = filename+".csv"
    f = open(join(root, fileName), 'a')
    f.write(",".join(line))
    f.write("\n")
    f.close()

def getProteinCompressionValues():
    writeCSVLine(["Domain",
                "blzpack_comp","bsc_comp","bzip2_comp", "Cmix_comp","GeCo3_comp","gzip_comp",
                "lizard_comp","lz4_comp","lzop_comp","mbgc_comp","MFCompress_comp",
                "naf_comp", "NUHT_comp","snzip_comp","UHT_comp","zip_comp","xz_comp","zstd_comp"], proName)
    for domain in listdir(mergedPathProteins):
        print(f"Computing {domain}...")
        tmpPath = join(tmpDir,domain)
        print(tmpPath)
        if not path.exists(tmpPath):
            makedirs(tmpPath)

        for fileName in listdir(join(mergedPathProteins, domain)):
            name = fileName.replace(".amino.gz", "")
            csvEntry = {
                "Domain": domain,
                "blzpack_comp": None ,"bsc_comp": None,"bzip2_comp": None,"gzip_comp" : None ,
                "lizard_comp"  : None,"lz4_comp" : None ,"lzop_comp" : None, "snzip_comp" : None , "zip_comp" : None ,
                "xz_comp" : None ,"zstd_comp" : None 
            }
                
                
            print(f"\tFile {fileName}")
            filePath = join(mergedPathProteins, domain, fileName)
            if isfile(filePath):
                os.system(f'rm {tmpPath}/*')
                decompressedPath = filePath.replace(".gz", "")
                os.system(f'gzip -k -d {filePath}')
                os.system(f'mv {decompressedPath} {tmpPath}/decompressed.amino')
                os.system(f'zcat {filePath} | grep -v ">" | tr -d -c "ACGT" > {tmpPath}/GENOME_FILE')
                original_file=join(tmpPath,"GENOME_FILE")
                original_sz = os.path.getsize(original_file)
                decompressed_sz = os.path.getsize(join(tmpPath,"decompressed.amino"))
                os.system(f'gto_fasta_from_seq  < {tmpPath}/GENOME_FILE > {tmpPath}/x.fa')
                fa_sz=os.path.getsize(join(tmpPath,"x.fa"))
                
                #SNZIP
                os.system(f'snzip -k {tmpPath}/decompressed.amino')
                compressed_file=join(tmpPath,"decompressed.amino.sz")
                cmp_sz = os.path.getsize(compressed_file)
                csvEntry["snzip_comp"] = str(cmp_sz/decompressed_sz)

                #LZOP
                os.system(f'lzop -9 {tmpPath}/decompressed.amino -o {tmpPath}/seq.lzop')
                compressed_file=join(tmpPath,"seq.lzop")
                cmp_sz = os.path.getsize(compressed_file)
                csvEntry["lzop_comp"] = str(cmp_sz/decompressed_sz)

                #LZ4
                os.system(f'{compressors_path}lz4 -9  {tmpPath}/decompressed.amino {tmpPath}/seq.lz4')
                compressed_file=join(tmpPath,"seq.lz4")
                cmp_sz = os.path.getsize(compressed_file)
                csvEntry["lz4_comp"] = str(cmp_sz/decompressed_sz)

                #LIZARD 
                os.system(f'{compressors_path}lizard -49  {tmpPath}/decompressed.amino {tmpPath}/seq.lizard')
                compressed_file=join(tmpPath,"seq.lizard")
                cmp_sz = os.path.getsize(compressed_file)
                csvEntry["lizard_comp"] = str(cmp_sz/decompressed_sz)
                
                #BLZPACK
                os.system(f'{compressors_path}blzpack -9  {tmpPath}/decompressed.amino {tmpPath}/seq.blzpack')
                compressed_file=join(tmpPath,"seq.blzpack")
                cmp_sz = os.path.getsize(compressed_file)
                csvEntry["blazpack_comp"] = str(cmp_sz/decompressed_sz)
                
                #BSC
                os.system(f'{compressors_path}bsc e  {tmpPath}/decompressed.amino {tmpPath}/seq.bsc')
                compressed_file=join(tmpPath,"seq.bsc")
                cmp_sz = os.path.getsize(compressed_file)
                csvEntry["bsc_comp"] = str(cmp_sz/decompressed_sz)

                #ZSTD
                os.system(f'zstd -9  {tmpPath}/decompressed.amino -o {tmpPath}/seq.zstd')
                compressed_file=join(tmpPath,"seq.zstd")
                cmp_sz = os.path.getsize(compressed_file)
                csvEntry["zstd_comp"] = str(cmp_sz/decompressed_sz)

                #XZ
                os.system(f'xz -k -9 {tmpPath}/decompressed.amino')
                compressed_file=join(tmpPath,"decompressed.amino.xz")
                cmp_sz = os.path.getsize(compressed_file)
                csvEntry["xz_comp"] = str(cmp_sz/decompressed_sz)

                #ZIP
                os.system(f'zip -9 {tmpPath}/seq.zip {tmpPath}/decompressed.amino')
                compressed_file=join(tmpPath,"seq.zip")
                cmp_sz = os.path.getsize(compressed_file)
                csvEntry["zip_comp"] = str(cmp_sz/decompressed_sz)

                #BZIP2
                os.system(f'bzip2 -k -9 {tmpPath}/decompressed.amino')
                compressed_file=join(tmpPath,"decompressed.amino.bz2")
                cmp_sz = os.path.getsize(compressed_file)
                csvEntry["bzip2_comp"] = str(cmp_sz/decompressed_sz)

                #GZIP
                os.system(f'gzip -c {tmpPath}/decompressed.amino > {tmpPath}/seq_gzip.gz')
                compressed_file=join(tmpPath,"seq_gzip.gz")
                cmp_sz = os.path.getsize(compressed_file)
                csvEntry["gzip_comp"] = str(cmp_sz/decompressed_sz)


                os.system(f"rm {tmpPath}/*")
                writeCSVLine(csvEntry.values(), proName)

                os.system(f"rm {tmpPath}/*")
                writeCSVLine(csvEntry.values(), proName)


def getCompressionValues():
    
    writeCSVLine(["Domain",
                "blzpack_comp","bsc_comp","bzip2_comp", "Cmix_comp","GeCo3_comp","gzip_comp",
                "JARVIS_comp","lizard_comp","lz4_comp","lzop_comp","mbgc_comp","MFCompress_comp",
                "naf_comp", "NUHT_comp","snzip_comp","UHT_comp","zip_comp","xz_comp","zstd_comp"], Name)

    for domain in listdir(mergedPath):
        print(f"Computing {domain}...")
        tmpPath = join(tmpDir,domain)
        print(tmpPath)
        if not path.exists(tmpPath):
            makedirs(tmpPath)
        print("domain")
        for fileName in listdir(join(mergedPath, domain)):
            name = fileName.replace(".fna.gz", "")
            csvEntry = {
                "Domain": domain,
                "blzpack_comp" : None,"bsc_comp" : None,"bzip2_comp" : None, "Cmix_comp" : None,"GeCo3_comp" : None,"gzip_comp" : None,
                "JARVIS_comp" : None,"lizard_comp" : None,"lz4_comp" : None,"lzop_comp" : None,"mbgc_comp" : None,"MFCompress_comp" : None,
                "naf_comp" : None,"NUHT_comp" : None,"snzip_comp" : None, "UHT_comp" : None,"zip_comp" : None,
                "xz_comp" : None,"zstd_comp" : None
            }
                
                
            print(f"\tFile {fileName}")
            filePath = join(mergedPath, domain, fileName)
            if isfile(filePath):
                os.system(f'rm {tmpPath}/*')
                decompressedPath = filePath.replace(".gz", "")
                os.system(f'gzip -k -d {filePath}')
                os.system(f'mv {decompressedPath} {tmpPath}/decompressed.fna')
                os.system(f'zcat {filePath} | grep -v ">" | tr -d -c "ACGT" > {tmpPath}/GENOME_FILE')
                original_file=join(tmpPath,"GENOME_FILE")
                original_sz = os.path.getsize(original_file)
                decompressed_sz = os.path.getsize(join(tmpPath,"decompressed.fna"))
                os.system(f'gto_fasta_from_seq  < {tmpPath}/GENOME_FILE > {tmpPath}/x.fa')
                fa_sz=os.path.getsize(join(tmpPath,"x.fa"))
                
                #SNZIP
                os.system(f'snzip -k {tmpPath}/decompressed.fna')
                compressed_file=join(tmpPath,"decompressed.fna.sz")
                cmp_sz = os.path.getsize(compressed_file)
                csvEntry["snzip_comp"] = str(cmp_sz/decompressed_sz)
     

                #NAF
                os.system(f'ennaf --level 22  {tmpPath}/decompressed.fna -o {tmpPath}/seq.naf --temp-dir DIR')
                compressed_file=join(tmpPath,"seq.naf")
                cmp_sz = os.path.getsize(compressed_file)
                csvEntry["naf_comp"] = str(cmp_sz/decompressed_sz)

                #MBGC
                os.system(f'mbgc -c 3 -i {tmpPath}/decompressed.fna {tmpPath}/seq.mbgc')
                compressed_file=join(tmpPath,"seq.mbgc")
                cmp_sz = os.path.getsize(compressed_file)
                csvEntry["mbgc_comp"] = str(cmp_sz/decompressed_sz)

                #LZOP
                os.system(f'lzop -9 {tmpPath}/decompressed.fna -o {tmpPath}/seq.lzop')
                compressed_file=join(tmpPath,"seq.lzop")
                cmp_sz = os.path.getsize(compressed_file)
                csvEntry["lzop_comp"] = str(cmp_sz/decompressed_sz)

                #LZ4
                os.system(f'{compressors_path}lz4 -9  {tmpPath}/decompressed.fna {tmpPath}/seq.lz4')
                compressed_file=join(tmpPath,"seq.lz4")
                cmp_sz = os.path.getsize(compressed_file)
                csvEntry["lz4_comp"] = str(cmp_sz/decompressed_sz)

                #LIZARD 
                os.system(f'{compressors_path}lizard -49  {tmpPath}/decompressed.fna {tmpPath}/seq.lizard')
                compressed_file=join(tmpPath,"seq.lizard")
                cmp_sz = os.path.getsize(compressed_file)
                csvEntry["lizard_comp"] = str(cmp_sz/decompressed_sz)
                
                #BLZPACK
                os.system(f'{compressors_path}blzpack -9  {tmpPath}/decompressed.fna {tmpPath}/seq.blzpack')
                compressed_file=join(tmpPath,"seq.blzpack")
                cmp_sz = os.path.getsize(compressed_file)
                csvEntry["blazpack_comp"] = str(cmp_sz/decompressed_sz)
                
                #BSC
                os.system(f'{compressors_path}bsc e  {tmpPath}/decompressed.fna {tmpPath}/seq.bsc')
                compressed_file=join(tmpPath,"seq.bsc")
                cmp_sz = os.path.getsize(compressed_file)
                csvEntry["bsc_comp"] = str(cmp_sz/decompressed_sz)

                #ZSTD
                os.system(f'zstd -9  {tmpPath}/decompressed.fna -o {tmpPath}/seq.zstd')
                compressed_file=join(tmpPath,"seq.zstd")
                cmp_sz = os.path.getsize(compressed_file)
                csvEntry["zstd_comp"] = str(cmp_sz/decompressed_sz)

                #XZ
                os.system(f'xz -k -9 {tmpPath}/decompressed.fna')
                compressed_file=join(tmpPath,"decompressed.fna.xz")
                cmp_sz = os.path.getsize(compressed_file)
                csvEntry["xz_comp"] = str(cmp_sz/decompressed_sz)

                #ZIP
                os.system(f'zip -9 {tmpPath}/seq.zip {tmpPath}/decompressed.fna')
                compressed_file=join(tmpPath,"seq.zip")
                cmp_sz = os.path.getsize(compressed_file)
                csvEntry["zip_comp"] = str(cmp_sz/decompressed_sz)

                #BZIP2
                os.system(f'bzip2 -k -9 {tmpPath}/decompressed.fna')
                compressed_file=join(tmpPath,"decompressed.fna.bz2")
                cmp_sz = os.path.getsize(compressed_file)
                csvEntry["bzip2_comp"] = str(cmp_sz/decompressed_sz)

                #GZIP
                os.system(f'gzip -c {tmpPath}/decompressed.fna > {tmpPath}/seq_gzip.gz')
                compressed_file=join(tmpPath,"seq_gzip.gz")
                cmp_sz = os.path.getsize(compressed_file)
                csvEntry["gzip_comp"] = str(cmp_sz/decompressed_sz)

                #NUHT
                os.system(f'{compressors_path}NUHT_Compress  {tmpPath}/decompressed.fna')
                compressed_file=join(tmpPath,"decompress.nuht")
                cmp_sz = os.path.getsize(compressed_file)
                csvEntry["NUHT_comp"] = str(cmp_sz/decompressed_sz)

                #UHT
                os.system(f'{UHT_COMP}UHT_compress {tmpPath}/decompressed.fna')
                compressed_file=join(tmpPath,"decompress_1.uht")
                cmp_sz = os.path.getsize(compressed_file)
                csvEntry["UHT_comp"] = str(cmp_sz/decompressed_sz)

                #MFCompress
                os.system(f"{compressors_path}MFCompressC -3 {tmpPath}/x.fa")
                compressed_file=join(tmpPath,"x.fa.mfc")
                cmp_sz = os.path.getsize(compressed_file)
                csvEntry["MFCompress_comp"] = str(cmp_sz/fa_sz)

                #CMIX
                os.system(f"{compressors_path}cmix -c {tmpPath}/GENOME_FILE {tmpPath}/cmix.seq")
                compressed_file=join(tmpPath,"cmix.seq")
                cmp_sz = os.path.getsize(compressed_file)
                csvEntry["Cmix_comp"] = str(cmp_sz/decompressed_sz)

                #JARVIS
                os.system(f"JARVIS -l 3 {tmpPath}/GENOME_FILE ")
                compressed_file=join(tmpPath,"GENOME_FILE.jc")
                cmp_sz = os.path.getsize(compressed_file)
                csvEntry["JARVIS_comp"] = str(cmp_sz/original_sz)

                #GeCo3
                os.system(f"GeCo3 -v -l 3 {tmpPath}/GENOME_FILE | sed '1,6d' | sed '2d' > {tmpPath}/RESULTS_GECO")
                os.system(f"cat {tmpPath}/RESULTS_GECO")
                with open(join(tmpPath, "RESULTS_GECO"), 'r') as fp:
                    try:
                        compressed_file=join(tmpPath,"GENOME_FILE.co")
                        cmp_geco = os.path.getsize(compressed_file)
                        csvEntry["GeCo3_comp"] = str(cmp_geco/original_sz)
                    except:
                        print(f"WARNING: DNA of {name} was not computed! Maybe lack of free space in RAM!")
                        csvEntry["GeCo3_comp"] = "1"
                os.system(f"rm {tmpPath}/*")
                writeCSVLine(csvEntry.values(), Name)

                os.system(f"rm {tmpPath}/*")
                writeCSVLine(csvEntry.values(), Name)

def _initialize():
    if not path.exists(root):
        makedirs(root)
    if not path.exists(tmpDir):
        makedirs(tmpDir)
    if not path.exists(tmpSeqPath):
        makedirs(tmpSeqPath)
    if not path.exists(mergedPath):
        makedirs(mergedPath)

if __name__ == "__main__":
    if "/" in sys.argv[0]:
        print("ERROR: Please run this script inside of src/! There are relative paths defined in this code that need to be respected!")
    else:
        main()