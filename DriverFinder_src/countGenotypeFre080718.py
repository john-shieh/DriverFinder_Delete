#Author: John Shieh

import os
import os.path
import sys
import shutil
import time
import pandas as pd

def count_Genotype(parent, fileSelected, fullPath):
    samfile_in = fullPath
    samfile_out = fullPath.replace(".sam", "[noH].sam")
    command0 = "samtools view -Sq 1 " + samfile_in + ">" + samfile_out #remove header
    os.system(command0)

    #read lines
    plines = 0
    f_in = open(samfile_out, 'r')
    line = f_in.readline().replace("\n", "")
    genotypeList = []
    while (line):
        plines = plines + 1
        items = line.split("\t")
        r_rname = items[2]
        r_mapq = int(items[4])
        line=f_in.readline().replace("\n", "")
        ind= str(r_rname).index('-')
        genotype = str(r_rname)[ind+1:]
        genotypeList.append([r_rname, genotype, r_mapq])
    genomedf = pd.DataFrame(genotypeList)
    genomedf.columns = ['bigG', 'littleG', 'mapQ']
    #genomedf.to_csv('E:/DoubleSample/genotypes.csv')
    #print(genomedf.head())

    #Genotypes (littleG)
    genomedfcountl = (genomedf['littleG'].value_counts()).reset_index()
    total2 = genomedfcountl['littleG'].sum()
    mainGl = genomedfcountl.loc[0,'littleG']
    #Sub-genotypes (big G)
    genomedfcountb = (genomedf['bigG'].value_counts()).reset_index()
    total = genomedfcountb['bigG'].sum()
    mainG = genomedfcountb.loc[0,'bigG']
    if (int(total2) == int(mainGl)):
        print(genomedfcountl.head())
        print("Total Reads: ", total2)
        print("Percentage for first genotype: 100%" + "\n")

        print(genomedfcountb.head())
        print("Total Reads: ", total)
        print("Percentage for first genotype: 100%")
        firstBigG = genomedfcountb.loc[0, 'index']
        pred = firstBigG
    else:
        refPerl = str((mainGl/total2)*100) + "%"
        secondGl = genomedfcountl.loc[1,'littleG']
        refPerl2 = str((secondGl/total2)*100) + "%"
        print(genomedfcountl.head())
        print("Total Reads: ", total2)
        print("Percentage for first genotype: " + str(refPerl))
        print("Percentage for second genotype: " + str(refPerl2) + "\n")
        topG = str(genomedfcountl.loc[0,'index'])
        
        #16 genotypes (bigG)
        refPer = str((mainG/total)*100) + "%"
        print(genomedfcountb.head())
        print("Total Reads: ", total)
        print("Percentage for top genotype (16): " + str(refPer) + "\n")
        #genomedfcountb.to_csv('big.csv')
        firstBigG = genomedfcountb.loc[0, 'index']
        secondBigG = genomedfcountb.loc[1, 'index']
        indB = str(firstBigG).index('-')
        firstLittleG = str(firstBigG)[indB+1:]
        pred = ""
        if topG == firstLittleG:
            pred = firstBigG
        else:
            pred = secondBigG
        print ("Predicted Genotype: ", pred)
    
    #check
    if (int(total2) > int(mainGl)):
        secondGl = genomedfcountl.loc[1,'littleG']
        refPerl2 = str((secondGl/total2)*100) + "%"
    else:
        refPerl2 = "0%"

    #Finish program
    parent.logInfo("Predicted genotype: " + str(pred))
    if (parent != None):
        parent.enableOptions()
 
""" testing
    print(genomedfcountl.head())
    print("Total Reads: ", total2) 
    print("Percentage for first genotype: ", refPerl)
    print("Percentage for second genotype: ", refPerl2)

    genomedfcountl.to_csv('E:/DoubleSample/little.csv')

if len(nd_index)>0:
    output = pd.DataFrame(uniqueDICT)
    gene_result=(output['Gene'].value_counts()).reset_index()
    gene_result.columns=['Gene', 'Count']
    parent.logInfo(str(gene_result))
    outpath = inputFilesFolder + foldername + '/' + filename + 'DupRemoved.csv'
    output.to_csv(outpath, index = False)
    viewseq.getGetSequenceFile(outpath, filename + 'DupRemoved')

#filepath = "E:/IJ1916/IJ19.R_split16__viralAlign.sam"
#count_Genotype(filepath)
"""

def process(parent, fileSelected=None):
    def logInfo(message):
        if (parent != None):
            parent.logInfo(message)
        else:
            print(message)
        return

    def showStatus(message):
        if (parent != None):
            parent.showStatus(message)
        else:
            print(message)
        return
    parent.logInfo("Finding genotype...")
    if (fileSelected == None):
        print ('Please provide file name.')
        return False
    else:
        # get all the folders/files 
        #fileFullPath = os.listdir(fileSelected)
        fileFullPath = fileSelected
        count_Genotype(parent, fileSelected, fileFullPath)
""" testing of counting genotype
        for item_name in folders:
            full_path_name = filesFolder+"/"+item_name
            if os.path.isfile(full_path_name): # it is a file
                # skip non sample file folder
                continue
            else: # it is a sample file split folder
                count_Genotype(parent, filesFolder, full_path_name, item_name)
"""
if __name__ == '__main__':
    if not os.path.isfile("samtools.exe"):
        print("samtools not installed.")
    else:
        folder_name = "E:/SingleSam"
        process(folder_name)

