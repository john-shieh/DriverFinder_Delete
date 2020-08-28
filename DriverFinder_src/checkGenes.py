# -*- coding: utf-8 -*-

import os, csv
from tkinter import filedialog, ttk
import platform
# from tkinter import *
import tkinter as tk
from tkinter import font as tkfont
import pandas as pd
from openpyxl import load_workbook

version="1.0.0"
hexversion=0x010000

class CheckGenes():
    platformType = platform.system()
    app_version = "1.0"
    app_title = "Check Genes"
    database_ok = False

    def __init__(self, parent):
        self.p = parent
        self.title = self.app_title + " " + self.app_version
        self.logInfo("Driver database loading...")
        result = self.loadGeneDatabase()
        if (result):
            self.logInfo("Driver database loading complete!")

    def quit_Application(self):
        root.destroy()
        return 1

    def logInfo(self, message):
        self.p.logInfo(message)
        return
    
    def showData(self, message):
        self.p.showData(message)
        return
    
    def loadGeneDatabase(self):
        
        file = 'DriverGeneDatabase.xlsx'
        if not os.path.isfile(file):
            self.logInfo("Error: database file "+file+" not found!")
            return False
        df = pd.read_excel(file, 'Driver Genes')
        self.GeneList=[]
        self.GeneListAll=[]
        for i in df.index:
            gene=df['Gene'][i]
            """
            try:
                risk=str(int(df['Risk'][i]))
            except:
                risk=''
            """
            risk='' 
            try:
                support_reads=int(df['Supporting Reads'][i])
            except:
                support_reads=0
            #print(gene, support_reads, risk)
            self.GeneList.append(gene)
            self.GeneListAll.append([gene, support_reads, risk])
        self.database_ok = True
        return True
        
    def evaluateRisk(self, filename):
        if not self.database_ok:
            self.logInfo("Error: no database loaded for risk evaluation!")
            return
        self.checkRiskWindow = tk.Toplevel()
        self.checkRiskWindow.wm_attributes("-topmost", 1)
        self.checkRiskWindow.title("Check Genes - "+filename)
        self.frame1 = tk.Frame(self.checkRiskWindow)
        self.frame1.grid(row=0, column=0, sticky='nw', columnspan=1)
        disp_font=tkfont.Font(family="Monaco", size=10)
        self.label1 = tk.Label(self.frame1, text='found        gene', font = disp_font)
        self.label1.grid(row=0, column=0, sticky='nw', columnspan=12)
        """
        self.label2 = tk.Label(self.frame1, text='----- ---- -------------------------')
        self.label2.grid(row=2, column=0, sticky='nw', columnspan=12)
        """
        # gene
        boxHeight = 20
        self.frame3 = tk.Frame(self.checkRiskWindow)
        self.frame3.grid(row=2, column=0, sticky='nw', columnspan=1)
        self.seqBox=tk.Text(self.frame3,wrap=tk.NONE, height=boxHeight,width=65,relief=tk.SUNKEN, font = disp_font)
        self.seqBox.grid(row=0,column=0,sticky='nw')
        self.scrollbar1=ttk.Scrollbar(self.frame3,command=self.seqBox.yview)
        self.seqBox.config(yscrollcommand=self.scrollbar1.set)
        self.scrollbar1.grid(row=0,column=1,sticky='ns')
        self.loadFileAndCheck(filename)
        return

    def changeLabel(self, label, message):
        label.configure(state=tk.NORMAL)
        label.delete(1.0, tk.END)
        label.insert(tk.INSERT, message)
        label.configure(state=tk.DISABLED)

    def loadFileAndCheck(self, filename):
        #self.GeneList=[]
        self.patientGene=[]
        self.seqBox.delete(1.0, tk.END)
        self.seqBox.update()       
        print("load file: " + filename)
        # 0-9:  'ReadName','ReadLength','Sequence','Chromosome','Gene','InsideGene','DistanceToGene','GeneDirection','Focus','GeneObj',
        # 10-19:'FocusObj','HostLocalCords','Hlength','HostRefCords','HostOrientation','ViralAccession','ViralLocalCords','Vlength','ViralRefCords','ViralOrientation',
        # 20-29:'HTM','HTMAdjusted','VTM','VTMAdjusted','Overlap','OverlapTM','OverlapTMAdjusted','Inserted','Microhomology','HostMapFlag',
        # 30-39:'ViralMapFlag','HMapQ','VMapQ','FastqSource','Index']
        with open(filename, 'r') as csvfile:
            reader = csv.reader(csvfile, delimiter=',')
            read_index = 0
            ndx = -1
            for row in reader:
                Gene = row[4]
                ndx += 1
                if (ndx == 0):
                    continue  # skip column title

                if not(Gene in self.patientGene):
                    self.patientGene.append(Gene)
                    if (Gene in self.GeneList):
                        ndx=self.GeneList.index(Gene)
                        risk=self.GeneListAll[ndx][2]
                        if risk=='':
                            disp_info = " >>>>        "+Gene
                        else:
                            disp_info = " >>>>  !!!!  "+Gene
                    else:
                        disp_info = "             "+Gene
                        
                    self.seqBox.insert(tk.END, disp_info+"\n")
        self.seqBox.update()
        return
