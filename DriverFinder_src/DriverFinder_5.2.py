# -*- coding: utf-8 -*-

import platform
import os
from tkinter import filedialog,ttk
from tkinter import messagebox
import tkinter as tk
from tkinter import font as tkfont
import GetSummaryAndDrivers
#import GetSummaryAndDrivers_B
import CombineSam
import countGenotypeFre080718
import changeRef

platformTest=platform.system()
pysam_ok = False
if platformTest == "Windows":
    #print("Detect HBV Mutation does not support in Windows!!!")
    pysam_ok = True
    import FindHBVmutations_John_0717 as FindHBVmutations
else:
    try:
        import pysam
        pysam_ok = True
        import FindHBVmutations_John_0717 as FindHBVmutations
    except:
        print("Detect HBV Mutation does not support due to pysam module not installed!!!")
        
import viewAlignedSequence
import checkGenes
import threading

version="5.2.0"
hexversion=0x050100


class MainGUI(tk.Frame):
    global version
    platformType=platform.system()
    app_version=version
    app_title = "DriverFinder"
    bttnColor="blue"
    folderSelected=""
    fileSelected=""
    HBVSelected=""
    rootFrame=tk.Frame
    
    # custom sort items
    def __init__(self, master=None):
        tk.Frame.__init__(self, master)
        self.title = self.app_title+" "+self.app_version
        self.master.title(self.title)
        self.create_widgets()
        self.master=master
        """
        if self.platformType == "Windows":
            self.logInfo("\n------------------------------------------------------")
            self.logInfo("Detect HBV Mutations only works under Mac environment!")
            self.logInfo("------------------------------------------------------\n")
        else:
            if (pysam_ok == False):
                self.logInfo("\n-------------------------------------------------------------------")
                self.logInfo("Detect HBV Mutations is disabled due to pysam module not installed!")
                self.logInfo("-------------------------------------------------------------------\n")
        """    
        #root.protocol('WM_DELETE_WINDOW', self.quit_Application)

    def quit_Application(self):
        root.destroy()
        return 1

    def selectFilesFolder(self):
        temp=filedialog.askdirectory(title='Select Folder')
        if (temp != ""):
            self.folderSelected = temp
            self.logInfo("Folder selected: "+self.folderSelected)
        return
    
    def selectFile(self):
        if self.platformType == "Windows":
            #temp=filedialog.askopenfilename(title='Select a ChimericSeq Output File',filetypes=(("Output Files", "*.csv"),("All files", "*.*")))
            temp=filedialog.askopenfilename(title='Select a ChimericSeq Output File')
        else:
            temp=filedialog.askopenfilename(title='Select a ChimericSeq Output File')
        if (temp != ""):
            self.fileSelected = temp
            self.logInfo("File selected: "+self.fileSelected)
        return

    def selectHBV(self):
        if self.platformType == "Windows":
            temp=filedialog.askopenfilename(title='Select HBV Reference File',filetypes=(("HBV Reference", "*.fa"),("All files", "*.*")))
        else:
            temp=filedialog.askopenfilename(title='Select HBV Reference File')
        if (temp != ""):
            self.HBVSelected = temp
            self.logInfo("HBV Reference File selected: "+self.HBVSelected)
        return
    
    
    def buildDrivers(self):
        if GetSummaryAndDrivers.checkRequirements(self) == False:
            self.logInfo("Stop getting summary and drivers")
            return
        try:
            self.logInfo("\nStart to get summary and drivers.\n")
            self.disableOptions()
            t = threading.Thread(target = GetSummaryAndDrivers.process, args = (self, self.folderSelected,))
            t.start()
        except:
            self.logInfo("Can't start thread to run GetSummaryAndDrivers.")
            self.enableOptions()
        return

    def buildDrivers_B(self):
        if GetSummaryAndDrivers_B.checkRequirements(self) == False:
            self.logInfo("Stop getting summary and drivers")
            return
        try:
            self.logInfo("\nStart to get summary and drivers.\n")
            self.disableOptions()
            t = threading.Thread(target = GetSummaryAndDrivers_B.process, args = (self, self.folderSelected,))
            t.start()
        except:
            self.logInfo("Can't start thread to run GetSummaryAndDrivers.")
            self.enableOptions()
        return

    def detectMutations(self):
        folderSelected = None
        HBVSelected = None
        if (self.folderSelected != ""):
            print("mutation files folder selected = ", self.folderSelected)
            folderSelected = self.folderSelected
        else:
            print("mutation files folder is default 'viralsamfiles'")
        if (self.HBVSelected != ""):
            HBVSelected = self.HBVSelected
        if FindHBVmutations.checkRequirements(self, folderSelected) == False:
            self.logInfo("Stop HBV mutations detection")
            return
        try:
            self.disableOptions()
            t = threading.Thread(target = FindHBVmutations.process, \
                                 args = (self, folderSelected, HBVSelected, ))
            t.start()
        except:
            self.logInfo("Can't start thread to run HBVmutations.")
            self.enableOptions()
        return

    def combineSamFiles(self):
        folderSelected = None
        if (self.folderSelected != ""):
            print("Folder selected with SAM files to combine = ", self.folderSelected)
            folderSelected = self.folderSelected
        try:
            self.disableOptions()
            t = threading.Thread(target = CombineSam.process, \
                                 args = (self, folderSelected))
            t.start()
        except:
            self.logInfo("Can't start thread to run CombineSam.")
            self.enableOptions()
        return

    def findGenotype(self):
        if (self.fileSelected != ""):
            print("File selected for genotyping: ", self.fileSelected)
            fileSelected = self.fileSelected
        try:
            self.disableOptions()
            t = threading.Thread(target = countGenotypeFre080718.process, \
                                 args = (self, fileSelected))
            t.start()
        except:
            self.logInfo("Can't start thread to run countGenotypeFre")
            self.enableOptions()
        return

    def changeReference(self):
        if(self.fileSelected != ""):
            print("Reference selected: ", self.fileSelected)
            fileSelected = self.fileSelected
        try:
            self.disableOptions()
            t = threading.Thread(target = changeRef.process, \
                                 args = (self, fileSelected))
            t.start()
        except:
            self.logInfo("Can't start thread to run changeRef")
            self.enableOptions()
    def viewReads(self):
        if self.fileSelected != "":
            self.disableOptions()
            vv=viewAlignedSequence.ViewAlignedSequence(self)
            self.logInfo("Load file")
            vv.loadFile(self.fileSelected)
            self.logInfo("Display reads")
            vv.alignSort2()
            self.enableOptions()
        else:
            self.logInfo("Error: no file selected!")
        return

    def evaluateRisk(self):
        if self.fileSelected != "":
            self.disableOptions()
            er=checkGenes.CheckGenes(self)
            er.evaluateRisk(self.fileSelected)
            self.enableOptions()
        else:
            self.logInfo("Error: no file selected!")
        return

    def showData(self, data):
        """
        self.dataText.configure(state=tk.NORMAL)
        self.dataText.insert(tk.END,data+'\n')
        self.dataText.configure(state=tk.DISABLED)
        self.dataText.see(tk.END)
        """
        return

    def logInfo(self, message):
        self.logText.configure(state=tk.NORMAL)
        self.logText.insert(tk.END,message+'\n')
        self.logText.configure(state=tk.DISABLED)
        self.logText.see(tk.END)
        return

    def showStatus(self,message):
        self.statusText.configure(state=tk.NORMAL)
        self.statusText.delete(1.0,tk.END)
        self.statusText.insert(tk.INSERT,message)
        self.statusText.configure(state=tk.DISABLED)

    def disableOptions(self):
        self.SummaryBttn.configure(state=tk.DISABLED)
        self.SummaryBttn_B.configure(state=tk.DISABLED)
        self.HBVBttn.configure(state=tk.DISABLED)
        self.ViewBttn.configure(state=tk.DISABLED)
        self.checkRiskBttn.configure(state=tk.DISABLED)
        self.combineSamBttn.configure(state=tk.DISABLED)
        self.findGenotypeBttn.configure(state=tk.DISABLED)
        return

    def enableOptions(self):
        global pysam_ok
        self.SummaryBttn.configure(state=tk.NORMAL)
        self.SummaryBttn_B.configure(state=tk.NORMAL)
        #if (self.platformType == "Windows") or (pysam_ok == False):
        if (pysam_ok == False):
            self.HBVBttn.configure(state=tk.DISABLED)
        else:
            self.HBVBttn.configure(state=tk.NORMAL)
        self.ViewBttn.configure(state=tk.NORMAL)
        self.checkRiskBttn.configure(state=tk.NORMAL)
        self.combineSamBttn.configure(state=tk.NORMAL)
        self.findGenotypeBttn.configure(state=tk.NORMAL)
        return
        
        
    def create_widgets(self):
        global pysam_ok
        self.frame1=tk.Frame(self.master)
        self.frame1.grid(row=0,column=0,sticky='nw',columnspan=10)         
        root.style=ttk.Style()
        root.style.configure('TButton', bg='blue')
        self.folderBttn=ttk.Button(self.frame1,text = ' Select Input Files Folder ', 
                                width=80, command=self.selectFilesFolder)
        self.folderBttn.grid(row=0,column=0,sticky='',ipadx=3,ipady=2)
        self.fileBttn=ttk.Button(self.frame1,text = ' Select Input File ',
                                width=80, command=self.selectFile)
        self.fileBttn.grid(row=1,column=0,sticky='',ipadx=3,ipady=2)
        self.HBVBttn=ttk.Button(self.frame1,text = ' Select HBV Reference File ',
                                width=80, command=self.selectHBV)
        self.HBVBttn.grid(row=2,column=0,sticky='',ipadx=3,ipady=2)
        self.frame2=tk.Frame(self.master)
        self.frame2.grid(row=1,column=0,sticky='nw',columnspan=1)         
        self.SummaryBttn=ttk.Button(self.frame2,text = 'Build\nSummary/Identify\nDrivers ',
                                width=20, command=self.buildDrivers)
        self.SummaryBttn.grid(row=0,column=0,sticky='',ipadx=3,ipady=2)
        """
        self.SummaryBttn_B=ttk.Button(self.frame2,text = 'Build\nSummary/Identify\nDrivers v2',
                                width=20, command=self.buildDrivers_B)
        self.SummaryBttn_B.grid(row=1,column=0,sticky='',ipadx=3,ipady=2)
        """
        self.sp1=ttk.Separator(self.frame2, orient=tk.HORIZONTAL)
        self.sp1.grid(row=2, column=0)
        self.HBVBttn=ttk.Button(self.frame2,text = ' Detect HBV \n Mutations ', 
                                width=20, command=self.detectMutations)
        self.HBVBttn.grid(row=3,column=0,sticky='',ipadx=3,ipady=2)
        self.sp2=ttk.Separator(self.frame2, orient=tk.HORIZONTAL)
        self.sp2.grid(row=4, column=0)
        # if (self.platformType == "Windows") or (pysam_ok == False):
        # print("FSS create_widgets skip Windows check for HBV mutation.....")
        if (pysam_ok == False):
            self.HBVBttn.configure(state=tk.DISABLED)
        self.ViewBttn=ttk.Button(self.frame2,text = ' Visualize \n Reads ', 
                                width=20, command=self.viewReads)
        self.ViewBttn.grid(row=5,column=0,sticky='',ipadx=3,ipady=2)
        self.sp3=ttk.Separator(self.frame2, orient=tk.HORIZONTAL)
        self.sp3.grid(row=5, column=0)
        self.checkRiskBttn=ttk.Button(self.frame2,text = ' Check \n Genes ',
                                width=20, command=self.evaluateRisk)
        self.checkRiskBttn.grid(row=7,column=0,sticky='',ipadx=3,ipady=2)
        
        #Combine Sam Files Button
        self.combineSamBttn=ttk.Button(self.frame2,text = ' Combine SAM Files ',
                                width=20, command=self.combineSamFiles)
        self.combineSamBttn.grid(row=8,column=0,sticky='',ipadx=3,ipady=2)

        #Find/Count genotype button
        self.findGenotypeBttn=ttk.Button(self.frame2, text = ' Find Genotype ',
                                width=20, command=self.findGenotype)
        self.findGenotypeBttn.grid(row=9, column=0, sticky='', ipadx=3, ipady=2)

        #Change the reference
        self.changeRefBttn=ttk.Button(self.frame2, text = ' Create Circular ',
                                width = 20, command=self.changeReference)
        self.changeRefBttn.grid(row=10, column=0, sticky='', ipadx=3, ipady=2)
        
        self.frame3=tk.Frame(self.master)
        self.frame3.grid(row=1,column=1,sticky='nw',columnspan=1)         
        self.logText= tk.Text(self.frame3,wrap=tk.NONE, height=25,width=75,relief=tk.SUNKEN)
        self.logText.grid(row=0,column=0)         
        self.scrollbar1=ttk.Scrollbar(self.frame3,command=self.logText.yview)
        self.scrollbar2=ttk.Scrollbar(self.frame3,orient=tk.HORIZONTAL,command=self.logText.xview)
        self.logText.config(xscrollcommand=self.scrollbar2.set)
        self.logText.config(yscrollcommand=self.scrollbar1.set)
        self.scrollbar1.grid(row=0,column=1,sticky='ns')
        self.scrollbar2.grid(row=1,column=0,sticky='ew')
        
        self.logText= tk.Text(self.frame3,wrap=tk.NONE, height=10,width=75,relief=tk.SUNKEN)
        self.logText.grid(row=3,column=0)         
        self.scrollbar1=ttk.Scrollbar(self.frame3,command=self.logText.yview)
        self.scrollbar2=ttk.Scrollbar(self.frame3,orient=tk.HORIZONTAL,command=self.logText.xview)
        self.logText.config(xscrollcommand=self.scrollbar2.set)
        self.logText.config(yscrollcommand=self.scrollbar1.set)
        self.scrollbar1.grid(row=3,column=1,sticky='ns')
        self.scrollbar2.grid(row=4,column=0,sticky='ew')
        
        self.frame4=tk.Frame(self.master)
        self.frame4.grid(row=2,column=1,sticky='nw',columnspan=2)         
        self.statusText= tk.Text(self.frame4,wrap=tk.NONE, height=1,width=75,relief=tk.SUNKEN)
        self.statusText.grid(row=0,column=0)         
        return

root = tk.Tk()
GUI = MainGUI(root)
GUI.mainloop()

