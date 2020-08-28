import os
import os.path
import sys
import shutil
import time
import pandas as pd

#direc 0 = copyend
#direc 1 = copybeg

def change_Ref(parent, fileSelected, fullPath, direc, size):
    fasta_in = fullPath
    if direc == 0:
        fasta_out = fullPath.replace(".fa","[copyend].fa")
    if direc == 1:
        fasta_out = fullPath.replace(".fa","[copybeg].fa")
    if direc == 2:
        fasta_out = fullPath.replace(".fa","[cut1600].fa")
    f_in = open(fasta_in, 'r')
    line = f_in.readline().replace("\n", "")
    newLine = ""
    while (line):
        if not line[0] == ">":
            newLine = newLine + line
        else:
            header = line
        
        line = f_in.readline().replace("\n", "")
    #print(newLine)
    f_in.close()
    f_out = open(fasta_out,'w')
    if direc == 0:
        newLine = newLine[-(size):len(newLine)] + newLine
    if direc == 1:
        newLine = newLine + newLine[0:(size)]
    if direc == 2:
        newLine = newLine[1600:len(newLine)] + newLine[0:1600]
    f_out.writelines(header + '\n')
    for i in range(0,len(newLine),80):
        if (i+80)<len(newLine):
            f_out.writelines(newLine[i:i+80] + '\n')
        else:
            f_out.writelines(newLine[i:len(newLine)] + '\n')
    f_out.close()

    parent.logInfo("Created the reference: " + str(fasta_out))
    parent.enableOptions()

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
    parent.logInfo("Changing reference...")
    if (fileSelected == None):
        print ('Please provide folder name.')
        return False
    else:
        # get all the folders/files 
        #fileFullPath = os.listdir(fileSelected)
        fileFullPath = fileSelected
        change_Ref(parent, fileSelected, fileFullPath, 0, 600)
        
#filepath = "E:/ChimericSeq_1_14_5/ViralReference16/X02763/X02763.fa"
#change_Ref(filepath, 2, 600)
