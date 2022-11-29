#!python
#!/usr/bin/env python

'''
  Usage: python3 compress_baseline.py
'''

import os
import sys
from os import listdir, path, makedirs
from os.path import isfile, join, getsize
sys.path.append(os.path.abspath('../src'))
from pprint import pprint
from collections import Counter
import itertools
import time

compressors_path = "../compressors/"
sequencesDir = "./original_sequences/"
mergedPath="../aux/original_sequences/"
root = "../baseline_features/"
UHT_COMP = "../compressors_files/Unbalanced-Huffman-Tree/dist/UHT_compress/"



Name = "compression_time"
percentages = [0,1,2,4,6,8,10]
tmpDir = os.path.join(root, "tmp/")
tmpSeqPath = join("..","tmpSq")

def main():
    _initialize()
    csvContent = getCompressionValues()

def writeCSVLine(line):
    fileName = Name+".csv"
    f = open(join(root, fileName), 'a')
    f.write(",".join(line))
    f.write("\n")
    f.close()

def getCompressionValues():

    writeCSVLine(["Domain","blzpack_comp","bsc_comp","bzip2_comp","Cmix_comp","GeCo3_comp","gzip_comp",
            "JARVIS_comp","lizard_comp","lz4_comp","lzop_comp","mbgc_comp","MFCompress_comp",
            "naf_comp","NUHT_comp", "snzip_comp", "UHT_comp","zip_comp","xz_comp","zstd_comp",
            "blzpack_time","bsc_time","bzip2_time","Cmix_time","GeCo3_time","gzip_time","JARVIS_time",
            "lizard_time", "lz4_time", "lzop_time", "mbgc_time","MFCompress_time",
            "naf_time","NUHT_time", "snzip_time", "UHT_time", "zip_time", "xz_time", "zstd_time"])
   
    for domain in listdir(mergedPath):
        print(f"Computing {domain}...")
        tmpPath = join(tmpDir,domain)
        print(tmpPath)
        if not path.exists(tmpPath):
            makedirs(tmpPath)

        for fileName in listdir(join(mergedPath, domain)):
            name = fileName.replace(".fna.gz", "")
            csvEntry = {
                "Domain": domain,
                "blzpack_comp": None,"bsc_comp": None,"bzip2_comp": None, "Cmix_comp": None,"GeCo3_comp": None,"gzip_comp": None,
                "JARVIS_comp": None,"lizard_comp": None,"lz4_comp": None,"lzop_comp": None,"mbgc_comp": None,"MFCompress_comp": None,
                "naf_comp": None,"NUHT_comp": None,"snzip_comp": None, "UHT_comp": None,"zip_comp": None,
                "xz_comp": None,"zstd_comp": None,
                "blzpack_time": None,"bsc_time": None,"bzip2_time": None, "Cmix_time": None,"GeCo3_time": None,"gzip_time": None,
                "JARVIS_time": None,"lizard_time": None, "lz4_time": None, "lzop_time": None, "mbgc_time": None,"MFCompress_time": None,
                "naf_time": None,"NUHT_time": None, "snzip_time": None, "UHT_time": None, "zip_time": None,
                "xz_time": None, "zstd_time": None
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

                #NAF
                start_naf = time.time()
                os.system(f'ennaf --level 22  {tmpPath}/decompressed.fna -o {tmpPath}/seq.naf --temp-dir DIR')
                end_naf = time.time()
                compressed_file=join(tmpPath,"seq.naf")
                cmp_sz = os.path.getsize(compressed_file)
                elapsed_time_naf=end_naf-start_naf
                csvEntry["naf_comp"] = str(cmp_sz/decompressed_sz)
                csvEntry["naf_time"] = str(elapsed_time_naf)             

                #MBGC
                start_mbgc = time.time()
                os.system(f'mbgc -c 3 -i {tmpPath}/decompressed.fna {tmpPath}/seq.mbgc')
                end_mbgc = time.time()
                compressed_file=join(tmpPath,"seq.mbgc")
                cmp_sz = os.path.getsize(compressed_file)
                elapsed_time_mbgc=end_mbgc-start_mbgc
                csvEntry["mbgc_comp"] = str(cmp_sz/decompressed_sz)
                csvEntry["mbgc_time"] = str(elapsed_time_mbgc)

                #LZOP
                start_lzop = time.time()
                os.system(f'lzop -9 {tmpPath}/decompressed.fna -o {tmpPath}/seq.lzop')
                end_lzop = time.time()
                compressed_file=join(tmpPath,"seq.lzop")
                cmp_sz = os.path.getsize(compressed_file)
                elapsed_time_lzop=end_lzop-start_lzop
                csvEntry["lzop_comp"] = str(cmp_sz/decompressed_sz)
                csvEntry["lzop_time"] = str(elapsed_time_lzop)

                #LZ4
                start_lz4 = time.time()
                os.system(f'{compressors_path}lz4 -9  {tmpPath}/decompressed.fna {tmpPath}/seq.lz4')
                end_lz4 = time.time()
                compressed_file=join(tmpPath,"seq.lz4")
                cmp_sz = os.path.getsize(compressed_file)
                elapsed_time_lz4=end_lz4-start_lz4
                csvEntry["lz4_comp"] = str(cmp_sz/decompressed_sz)
                csvEntry["lz4_time"] = str(elapsed_time_lz4)

                #LIZARD 
                start_lizard = time.time()
                os.system(f'{compressors_path}lizard -49  {tmpPath}/decompressed.fna {tmpPath}/seq.lizard')
                end_lizard = time.time()
                compressed_file=join(tmpPath,"seq.lizard")
                cmp_sz = os.path.getsize(compressed_file)
                elapsed_time_lizard=end_lizard-start_lizard
                csvEntry["lizard_comp"] = str(cmp_sz/decompressed_sz)
                csvEntry["lizard_time"] = str(elapsed_time_lizard)
                
                #BLZPACK
                start_blazpack = time.time()
                os.system(f'{compressors_path}blzpack -9  {tmpPath}/decompressed.fna {tmpPath}/seq.blzpack')
                end_blazpack = time.time()
                compressed_file=join(tmpPath,"seq.blzpack")
                cmp_sz = os.path.getsize(compressed_file)
                elapsed_time_blazpack=end_blazpack-start_blazpack
                csvEntry["blazpack_comp"] = str(cmp_sz/decompressed_sz)
                csvEntry["blazpack_time"] = str(elapsed_time_blazpack)
                
                #BSC
                start_bsc = time.time()
                os.system(f'{compressors_path}bsc e  {tmpPath}/decompressed.fna {tmpPath}/seq.bsc')
                end_bsc = time.time()
                compressed_file=join(tmpPath,"seq.bsc")
                cmp_sz = os.path.getsize(compressed_file)
                elapsed_time_bsc=end_bsc-start_bsc
                csvEntry["bsc_comp"] = str(cmp_sz/decompressed_sz)
                csvEntry["bsc_time"] = str(elapsed_time_bsc)
        
                #ZSTD
                start_zstd = time.time()
                os.system(f'zstd -9  {tmpPath}/decompressed.fna -o {tmpPath}/seq.zstd')
                end_zstd = time.time()
                compressed_file=join(tmpPath,"seq.zstd")
                cmp_sz = os.path.getsize(compressed_file)
                elapsed_time_zstd=end_zstd-start_zstd
                csvEntry["zstd_comp"] = str(cmp_sz/decompressed_sz)
                csvEntry["zstd_time"] = str(elapsed_time_zstd)

                #XZ
                start_xz = time.time()
                os.system(f'xz -k -9 {tmpPath}/decompressed.fna')
                end_xz = time.time()
                compressed_file=join(tmpPath,"decompressed.fna.xz")
                cmp_sz = os.path.getsize(compressed_file)
                elapsed_time_xz=end_xz-start_xz
                csvEntry["xz_comp"] = str(cmp_sz/decompressed_sz)
                csvEntry["xz_time"] = str(elapsed_time_xz)


                #ZIP
                start_zip = time.time()
                os.system(f'zip -9 {tmpPath}/seq.zip {tmpPath}/decompressed.fna')
                end_zip = time.time()
                compressed_file=join(tmpPath,"seq.zip")
                cmp_sz = os.path.getsize(compressed_file)
                elapsed_time_zip=end_zip-start_zip
                csvEntry["zip_comp"] = str(cmp_sz/decompressed_sz)
                csvEntry["zip_time"] = str(elapsed_time_zip)

                #BZIP2
                start_bzip2 = time.time()
                os.system(f'bzip2 -k -9 {tmpPath}/decompressed.fna')
                end_bzip2 = time.time()
                compressed_file=join(tmpPath,"decompressed.fna.bz2")
                cmp_sz = os.path.getsize(compressed_file)
                elapsed_time_bzip2=end_bzip2-start_bzip2
                csvEntry["bzip2_comp"] = str(cmp_sz/decompressed_sz)
                csvEntry["bzip2_time"] = str(elapsed_time_bzip2)

                #GZIP
                start_gzip = time.time()
                os.system(f'gzip -c {tmpPath}/decompressed.fna > {tmpPath}/seq_gzip.gz')
                end_gzip = time.time()
                compressed_file=join(tmpPath,"seq_gzip.gz")
                cmp_sz = os.path.getsize(compressed_file)
                elapsed_time_gzip=end_gzip-start_gzip
                csvEntry["gzip_comp"] = str(cmp_sz/decompressed_sz)
                csvEntry["gzip_time"] = str(elapsed_time_gzip)
                
                #NUHT
                start_nuht = time.time()
                os.system(f'{compressors_path}NUHT_Compress  {tmpPath}/decompressed.fna')
                end_nuht = time.time()
                compressed_file=join(tmpPath,"decompress.nuht")
                cmp_sz = os.path.getsize(compressed_file)
                elapsed_time_nuht=end_nuht-start_nuht
                csvEntry["NUHT_comp"] = str(cmp_sz/decompressed_sz)
                csvEntry["NUHT_time"] = str(elapsed_time_nuht)
                
                #UHT
                start_uht = time.time()
                os.system(f'{UHT_COMP}UHT_compress {tmpPath}/decompressed.fna')
                end_uht = time.time()
                compressed_file=join(tmpPath,"decompress_1.uht")
                cmp_sz = os.path.getsize(compressed_file)
                elapsed_time_uht=end_uht-start_uht
                csvEntry["UHT_comp"] = str(cmp_sz/decompressed_sz)
                csvEntry["UHT_time"] = str(elapsed_time_uht)

                #MFCompress
                start_mfc = time.time()
                os.system(f"{compressors_path}MFCompressC -3 {tmpPath}/x.fa")
                end_mfc = time.time()
                compressed_file=join(tmpPath,"x.fa.mfc")
                cmp_sz = os.path.getsize(compressed_file)
                elapsed_time_mfc = end_mfc-start_mfc
                csvEntry["MFCompress_comp"] = str(cmp_sz/fa_sz)
                csvEntry["MFCompress_time"] = str(elapsed_time_mfc)

                #CMIX
                start_cmix = time.time()
                os.system(f"{compressors_path}cmix -c {tmpPath}/GENOME_FILE {tmpPath}/cmix.seq")
                end_cmix = time.time()
                compressed_file=join(tmpPath,"cmix.seq")
                cmp_sz = os.path.getsize(compressed_file)
                elapsed_time_cmix=end_cmix-start_cmix
                csvEntry["Cmix_comp"] = str(cmp_sz/decompressed_sz)
                csvEntry["Cmix_time"] = str(elapsed_time_cmix)

                #JARVIS
                start_jarvis = time.time()
                os.system(f"JARVIS -l 3 {tmpPath}/GENOME_FILE ")
                end_jarvis = time.time()
                compressed_file=join(tmpPath,"GENOME_FILE.jc")
                cmp_sz = os.path.getsize(compressed_file)
                elapsed_time_jarvis=end_jarvis-start_jarvis
                csvEntry["JARVIS_comp"] = str(cmp_sz/original_sz)
                csvEntry["JARVIS_time"] = str(elapsed_time_jarvis)
                
                #GeCo3
                start_geco = time.time()
                os.system(f"GeCo3 -v -l 3 {tmpPath}/GENOME_FILE | sed '1,6d' | sed '2d' > {tmpPath}/RESULTS_GECO")
                end_geco = time.time()
                elapsed_time_geco=end_geco-start_geco
                os.system(f"cat {tmpPath}/RESULTS_GECO")
                with open(join(tmpPath, "RESULTS_GECO"), 'r') as fp:
                    try:
                        compressed_file=join(tmpPath,"GENOME_FILE.co")
                        cmp_geco = os.path.getsize(compressed_file)
                        csvEntry["GeCo3_comp"] = str(cmp_geco/original_sz)
                        csvEntry["GeCo3_time"] = str(elapsed_time_geco)
                    except:
                        print(f"WARNING: DNA of {name} was not computed! Maybe lack of free space in RAM!")
                        csvEntry["GeCo3_comp"] = "1"
                os.system(f"rm {tmpPath}/*")
                writeCSVLine(csvEntry.values())



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