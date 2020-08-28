# -*- coding: utf-8 -*-

import subprocess,os,math,re,sys,zipfile,datetime,multiprocessing,webbrowser,threading,traceback,time,csv
import platform
import difflib
#from tkinter import *
from tkinter import filedialog,ttk
from tkinter import messagebox
import tkinter as tk
from threading import Thread,Lock
from Bio.Seq import Seq
from bisect import bisect
from Bio import SeqIO
from tkinter import font as tkfont

class viewSequence():
    platformType=platform.system()
    NA_value=999999999
    maxSeqLength=210 # FSS changed from 185
    maxWindowLength=210
    maxDataLabelLength=20
    seqLabelLength=0
    seqLabel=""
    app_version="1.11"
    app_title = "JBS ChimericSeq Viewer"
    rows=[]
    list_data=[]
    display_list_data=[]
    maxLength=0
    infoWindow=0
    sortBttnColor="light gray"
    selectedSortBttnColor="light green"
    customSortbttn_id=4
    savebttn_id=5
    undobttn_id=6
    alignViewbttn_id=7

    # custom sort items
    customSortWindow = 0
    helpWindow = 0
    listSelection = -1
    listSortingItems =[]
    NONE_ITEM="None"
    VIRAL="Viral"
    HOST="Host"
    ENTIRE_SEQ="Entire Seq"
    CHROMOSOME="Chromosome"
    HOST_REF_COORDS_START="HostRefCoords Start"
    HOST_REF_COORDS_STOP="HostRefCoords Stop"
    HOST_ORIENTATION="Host Orientation"
    VIRAL_REF_COORDS_START="ViralRefCoords Start"
    VIRAL_REF_COORDS_STOP="ViralRefCoords Stop"
    VIRAL_ORIENTATION="Viral Orientation"
    GENE="Gene"
    DISTANCE_TO_GENE="Distance To Gene"
    VIRAL_FIRST="Viral First"
    HOST_FIRST="Host First"
    sortItemsOptionList=(NONE_ITEM, VIRAL, HOST, ENTIRE_SEQ, CHROMOSOME, HOST_REF_COORDS_START, HOST_REF_COORDS_STOP,
                         HOST_ORIENTATION, VIRAL_REF_COORDS_START, VIRAL_REF_COORDS_STOP, VIRAL_ORIENTATION,
                         GENE, DISTANCE_TO_GENE, VIRAL_FIRST, HOST_FIRST)
    name = ''
    filename = ''
    

    def alignSort(self):
        if len(self.list_data)==0: # data not loaded
            # messagebox.showinfo("Error", "No file selected")
            print ('Error, No File Selected!')
            return
        list_data=self.list_data
        display_type=0
        item=None
        item=self.VIRAL_REF_COORDS_START
        if item==self.VIRAL_REF_COORDS_START:
            list_data=sorted(list_data, key=lambda item:item[20]) 
        item=self.HOST_REF_COORDS_STOP
        if item==self.HOST_REF_COORDS_STOP:
            list_data=sorted(list_data, key=lambda item:item[18]) 
        item=self.CHROMOSOME
        if item==self.CHROMOSOME:
            list_data=sorted(list_data, key=lambda item:item[15]) 
            list_data=sorted(list_data, key=lambda item:item[16])
        item=self.GENE
        if item==self.GENE:
            list_data=sorted(list_data, key=lambda item:item[23]) 
        if item==self.VIRAL:
            list_data=sorted(list_data, key=lambda item:item[5])
            list_data=sorted(list_data, key=lambda item:item[6])
        if item==self.HOST:
            list_data=sorted(list_data, key=lambda item:item[3])
            list_data=sorted(list_data, key=lambda item:item[4])
        if item==self.ENTIRE_SEQ:
            list_data=sorted(list_data, key=lambda item:item[1])
            list_data=sorted(list_data, key=lambda item:item[2])
        if item==self.HOST_REF_COORDS_START:
            list_data=sorted(list_data, key=lambda item:item[17]) 
        if item==self.HOST_ORIENTATION:
            list_data=sorted(list_data, key=lambda item:item[19]) 
        if item==self.VIRAL_REF_COORDS_STOP:
            list_data=sorted(list_data, key=lambda item:item[21]) 
        if item==self.VIRAL_ORIENTATION:
            list_data=sorted(list_data, key=lambda item:item[22]) 
        if item==self.DISTANCE_TO_GENE:
            list_data=sorted(list_data, key=lambda item:item[24]) 
        if item==self.VIRAL_FIRST:
            display_type=1
            list_data=sorted(list_data, key=lambda item:item[14]) 
        if item==self.HOST_FIRST:
            display_type=2
            list_data=sorted(list_data, key=lambda item:item[14], reverse=True) 

                
        self.displayAlignedSequence(list_data, display_type)
        # self.show_Info_Window("Sorting...", 0, "Custom Sorting")
        return

    def displayAlignedSequence(self, list_data, sorting_type=0):
        filename = self.name + "Align.html"

        """
        <!DOCTYPE html>
        <html>
        <style>
        s1 {background-color: yellow;}
        s2 {color: black;}
        s3 {color: red;}
        </style>
        <body>
        <font face="courier" size="2">
        <s2>normal</s2>
        <s1><s2>host</s2></s1> 
        <s1><s3><u>overlap</u></s3></s1> 
        <s3>viral</s3>
        <s2>I am black</s2>
        </font>
        </body>
        </html>
        """
        fo=open(filename, "w")
        fo.write("<!DOCTYPE html>\n")
        fo.write("<html>\n")
        fo.write("<style>\n")
        fo.write("s1 {background-color: yellow;}\n")
        fo.write("s2 {color: black;}\n")
        fo.write("s3 {color: red;}\n")
        fo.write("</style>\n")
        fo.write("<body>\n")
        fo.write('<font face="courier" size="3">\n')
        V1=self.seqLabel.index("Red")
        V2=V1+3
        H1=self.seqLabel.index("Highlighted")
        H2=H1+11
        O1=self.seqLabel.index("Underlined")
        O2=O1+10
        label=self.seqLabel
        column_label=label[:V1]+'<s3>'+label[V1:V2]+'</s3>'+label[V2:H1]+'<s1>'+label[H1:H2]+'</s1>'+label[H2:O1]+'<s1><s3><u>'+label[O1:O2]+'</u></s3></s1>'+label[O2:]
        column_label=column_label.replace(" ", "&nbsp;")+"<br>"
        fo.write(column_label+'\n')
        line=""
        for i in range(len(self.seqLabel)):
            line=line+"="
        line=line+"<br>"
        fo.write(line+'\n')


        self.display_list_data=list_data
        # self.seqBox.configure(state=tk.NORMAL)
        # self.seqBox.delete(1.0, tk.END)
        blanks=""
        for i in range(10):
            blanks=blanks+"                    "
        ub=0
        group_count=0
        group_symbols="*!+"
        group_symbol=group_symbols[0]
        group_symbol_index=0
        symbol_applied=False
        smallet_host_start=9999999999
        non_mapped_extra=0
        for i in range(len(list_data)):
            #print("current ", i, list_data[i][15] , list_data[i][23])
            if i > ub:
                smallet_host_start=list_data[i][17]
                group_count=0
                non_mapped_extra=0
                if list_data[i][14] == True:
                    #host first
                    non_mapped_extra=list_data[i][7] # host start
                else:
                    #host not first
                    non_mapped_extra=list_data[i][2] - list_data[i][8] # seq length - host stop
                for j in range(i+1, len(list_data)):
                    # move upper bound to next if same Chromosome, Gene, host ref stop
                    if (list_data[i][15] == list_data[j][15]) and \
                       (list_data[i][23] == list_data[j][23]) and \
                       (list_data[i][18] == list_data[j][18]):
                        ub=j
                        group_count=group_count+1
                        if list_data[j][17] < smallet_host_start:
                            smallet_host_start = list_data[j][17]
                        if list_data[j][14] == True:
                            #host first
                            non_mapped_extra_temp=list_data[j][7] # host start 
                        else:
                            #host not first
                            non_mapped_extra_temp=list_data[j][2] - list_data[j][8] # seq length - host stop
                        if non_mapped_extra_temp > non_mapped_extra:
                            non_mapped_extra = non_mapped_extra_temp
                    else:
                        break
                if (group_count > 0):
                    group_symbol_index+=1
                    group_symbol=group_symbols[group_symbol_index % len(group_symbols)]
                else:
                    non_mapped_extra=0
                    
            disp_seq=list_data[i][1]
            ndx=i+1
            extra_offset=0
            if list_data[i][17] > smallet_host_start:
                extra_offset= (list_data[i][17] - smallet_host_start)
            offset=16+extra_offset

            # check for flip over -----------------------------
            o_start = int(list_data[i][12])
            o_stop = int(list_data[i][13])
            h_start = int(list_data[i][7])
            h_stop = int(list_data[i][8])
            v_start = int(list_data[i][9])
            v_stop = int(list_data[i][10])
            if list_data[i][14] == False:  #not host first
                read_length = int(list_data[i][2])
                if (o_start != 0) or (o_stop != 0):
                    o_start = read_length - o_start
                    o_stop = read_length - o_stop
                    temp=o_start
                    o_start=o_stop
                    o_stop=temp
                h_start = read_length - h_start
                h_stop = read_length - h_stop                
                temp=h_start
                h_start=h_stop
                h_stop=temp
                v_stop = read_length - v_stop
                v_start = read_length - v_start
                temp=v_start
                v_start=v_stop
                v_stop=temp

                sequence_complement=str(Seq(disp_seq).reverse_complement())
                disp_seq = sequence_complement

            add_blanks_non_mapped = 0
            if non_mapped_extra > 0:
                add_blanks_non_mapped = (non_mapped_extra - h_start)

            H_start = str(h_start+offset+add_blanks_non_mapped)
            H_stop=str(h_stop+offset+add_blanks_non_mapped)
            V_start = str(v_start+offset+add_blanks_non_mapped)
            V_stop=str(v_stop+offset+add_blanks_non_mapped)
            O_start = str(o_start+offset+add_blanks_non_mapped)
            O_stop=str(o_stop+offset+add_blanks_non_mapped)

            # if a sequence is too long, it will be truncated.
            # so, we have to make sure that color will be truncated as well
            if (int(H_stop)>self.maxSeqLength):
                H_stop = str(self.maxSeqLength)
            if (int(V_stop)>self.maxSeqLength):
                V_stop = str(self.maxSeqLength)
            if (int(O_stop)>self.maxSeqLength):
                O_stop = str(self.maxSeqLength)
            overlap=list_data[i][11]
            rec_num=str(ndx)+"      "
            rec_num=rec_num[:6]+"["
            rec_num1=str(list_data[i][0])+"      "
            if list_data[i][14] == True:
                rec_num=rec_num+rec_num1[:6]+"]_"
            else:
                rec_num=rec_num+rec_num1[:6]+"]v"
            if (group_count>0):
                if (i==len(list_data) - 1): #last record
                    if (list_data[i][20] == list_data[i-1][20]):
                        rec_num=rec_num+group_symbol
                    else:
                        rec_num=rec_num+" "
                else: # not the last record
                    # same host stop, but different viral start
                    if (symbol_applied == False):
                        # same host stop & same viral start
                        if (list_data[i][20] == list_data[i+1][20]):
                            rec_num=rec_num+group_symbol
                            symbol_applied = True
                        else:
                            rec_num=rec_num+" "
                    else:
                        # symbol has been applied
                        if (list_data[i][20] == list_data[i-1][20]): # still same viral start
                            rec_num=rec_num+group_symbol
                        else:   # different viral start   
                            if (list_data[i][20] == list_data[i+1][20]): 
                                group_symbol_index+=1
                                group_symbol=group_symbols[group_symbol_index % len(group_symbols)]
                                rec_num=rec_num+group_symbol
                            else:
                                symbol_applied = False
                                group_symbol_index+=1
                                group_symbol=group_symbols[group_symbol_index % len(group_symbols)]
                                rec_num=rec_num+" "
            else:
                rec_num=rec_num+" "
            #disp_info=rec_num+disp_seq+"  "+list_data[i][15][:self.maxDataLabelLength]
            temp_info=rec_num
            for ii in range(0, extra_offset++add_blanks_non_mapped):
                temp_info = temp_info+" "
            temp_info=temp_info+disp_seq+blanks
            # Chromosome
            chromosome_info=list_data[i][15]+"                    "
            disp_info=temp_info[:(self.maxSeqLength-2)]+"    "+chromosome_info[:20]
            # Host Ref Cords & orientation
            H_ref_cords=" ["+("        "+str(list_data[i][17]))[-10:]+", "+("        "+str(list_data[i][18]))[-10:]+"] "
            disp_info=disp_info+H_ref_cords
            H_orientation = "  "+(list_data[i][19]+"     ")[:5]
            disp_info=disp_info+"      "+H_orientation+"        "
            # Viral Ref Cords & orientation
            V_ref_cords=" ["+("        "+str(list_data[i][20]))[-7:]+", "+("        "+str(list_data[i][21]))[-7:]+"] "
            disp_info=disp_info+V_ref_cords
            V_orientation = "  "+(list_data[i][22]+"     ")[:5]
            disp_info=disp_info+"         "+V_orientation+"    "
            # Gene
            Gene=(list_data[i][23]+"               ")[:17]
            disp_info=disp_info+"    "+Gene+"  "
            if list_data[i][24]!=self.NA_value:
                DistanceToGene=str(list_data[i][24])
            else:
                DistanceToGene="N/A"
            disp_info=disp_info+"    "+DistanceToGene+"     "
            # make sure disp_info has same legth as the seqLabel, so the bottom scroll bar can move
            # both seqLabel & disp_info at same time
            for i in range(len(disp_info), self.seqLabelLength):
                disp_info=disp_info+" "
            # add info into display 
            # self.seqBox.insert(tk.END, disp_info+"\n")
            # self.seqBox.tag_add('Host',str(ndx)+"."+H_start,str(ndx)+"."+H_stop)
            # self.seqBox.tag_add('Virus',str(ndx)+"."+V_start,str(ndx)+"."+V_stop)
            # if (O_start != O_stop):
            #     self.seqBox.tag_add('Overlap',str(ndx)+"."+O_start,str(ndx)+"."+O_stop)

            result_HTML=self.convertToHTML(disp_info, 0, H_start, H_stop, V_start, V_stop, O_start, O_stop)
            fo.write(result_HTML+"<br>\n")

        # update display
        # self.seqBox.update()

        # finish HTML file
        fo.write("</font>\n")
        fo.write("</body>\n")
        fo.write("</html>\n")
        fo.close()

        # self.set_Sort_Button_Color(self.alignViewbttn_id)

        self.filename = self.filename + "Align.html"
        return



    def saveToHTML(self, list_data, sorting_type=0, filename = 'ChimericSeqDisp.html'):
        """
        <!DOCTYPE html>
        <html>
        <style>
        s1 {background-color: yellow;}
        s2 {color: black;}
        s3 {color: red;}
        </style>
        <body>
        <font face="courier" size="2">
        <s2>normal</s2>
        <s1><s2>host</s2></s1> 
        <s1><s3><u>overlap</u></s3></s1> 
        <s3>viral</s3>
        <s2>I am black</s2>
        </font>
        </body>
        </html>
        """
        fo=open(filename, "w")
        fo.write("<!DOCTYPE html>\n")
        fo.write("<html>\n")
        fo.write("<style>\n")
        fo.write("s1 {background-color: yellow;}\n")
        fo.write("s2 {color: black;}\n")
        fo.write("s3 {color: red;}\n")
        fo.write("</style>\n")
        fo.write("<body>\n")
        fo.write('<font face="courier" size="3">\n')
        V1=self.seqLabel.index("Red")
        V2=V1+3
        H1=self.seqLabel.index("Highlighted")
        H2=H1+11
        O1=self.seqLabel.index("Underlined")
        O2=O1+10
        label=self.seqLabel
        column_label=label[:V1]+'<s3>'+label[V1:V2]+'</s3>'+label[V2:H1]+'<s1>'+label[H1:H2]+'</s1>'+label[H2:O1]+'<s1><s3><u>'+label[O1:O2]+'</u></s3></s1>'+label[O2:]
        column_label=column_label.replace(" ", "&nbsp;")+"<br>"
        fo.write(column_label+'\n')
        line=""
        for i in range(len(self.seqLabel)):
            line=line+"="
        line=line+"<br>"
        fo.write(line+'\n')
        blanks=""
        for i in range(10):
            blanks=blanks+"                    "
        for i in range(len(list_data)):
            disp_seq=list_data[i][1]
            ndx=i+1
            offset=16
            if (sorting_type == 1): #viral sorting
                if int(list_data[i][7]) <= int(list_data[i][9]): #host before viral
                    extra_offset=self.maxLength-len(disp_seq)
                    disp_seq = blanks[:extra_offset]+disp_seq
                    offset = offset+extra_offset
            if (sorting_type == 2): #host sorting
                if int(list_data[i][7]) >= int(list_data[i][9]): #host after viral
                    extra_offset=self.maxLength-len(disp_seq)
                    disp_seq = blanks[:extra_offset]+disp_seq
                    offset = offset+extra_offset
            H_start = str(list_data[i][7]+offset)
            H_stop=str(list_data[i][8]+offset)
            V_start = str(list_data[i][9]+offset)
            V_stop=str(list_data[i][10]+offset)
            O_start = str(list_data[i][12]+offset)
            O_stop=str(list_data[i][13]+offset)
            # if a sequence is too long, it will be truncated.
            # so, we have to make sure that color will be truncated as well
            if (int(H_stop)>self.maxSeqLength):
                H_stop = str(self.maxSeqLength)
            if (int(V_stop)>self.maxSeqLength):
                V_stop = str(self.maxSeqLength)
            if (int(O_stop)>self.maxSeqLength):
                O_stop = str(self.maxSeqLength)
            overlap=list_data[i][11]
            rec_num=str(ndx)+"      "
            rec_num=rec_num[:6]+"["
            rec_num1=str(list_data[i][0])+"      "
            rec_num=rec_num+rec_num1[:6]+"]_ "
            #disp_info=rec_num+disp_seq+"  "+list_data[i][15][:self.maxDataLabelLength]
            temp_info=rec_num+disp_seq+blanks
            # Chromosome
            chromosome_info=list_data[i][15]+"                    "
            disp_info=temp_info[:(self.maxSeqLength-2)]+"    "+chromosome_info[:20]
            # Host Ref Cords & orientation
            H_ref_cords=" ["+("        "+str(list_data[i][17]))[-10:]+", "+("        "+str(list_data[i][18]))[-10:]+"] "
            disp_info=disp_info+H_ref_cords
            H_orientation = "  "+(list_data[i][19]+"     ")[:5]
            disp_info=disp_info+"      "+H_orientation+"        "
            # Viral Ref Cords & orientation
            V_ref_cords=" ["+("        "+str(list_data[i][20]))[-7:]+", "+("        "+str(list_data[i][21]))[-7:]+"] "
            disp_info=disp_info+V_ref_cords
            V_orientation = "  "+(list_data[i][22]+"     ")[:5]
            disp_info=disp_info+"         "+V_orientation+"    "
            # Gene
            Gene=(list_data[i][23]+"               ")[:17]
            Gene=Gene.replace("-", "_")
            disp_info=disp_info+"    "+Gene+"  "
            if list_data[i][24]!=self.NA_value:
                DistanceToGene=str(list_data[i][24])
            else:
                DistanceToGene="N/A"
            disp_info=disp_info+"    "+DistanceToGene+"     "
            # make sure disp_info has same legth as the seqLabel, so the bottom scroll bar can move
            # both seqLabel & disp_info at same time
            for i in range(len(disp_info), self.seqLabelLength):
                disp_info=disp_info+" "
            # add info into display 
            result_HTML=self.convertToHTML(disp_info, 0, H_start, H_stop, V_start, V_stop, O_start, O_stop)
            fo.write(result_HTML+"<br>\n")

        fo.write("</font>\n")
        fo.write("</body>\n")
        fo.write("</html>\n")
        fo.close()
        return
    
    def convertToHTML(self, data, offset, H_start, H_stop, V_start, V_stop, O_start, O_stop):
        """
        <style>
        s1 {background-color: yellow;}
        s2 {color: black;}
        s3 {color: red;}
        </style>
        <p>I am&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;normal</p>
        <s2>normal</s2>
        <s1><s2>host</s2></s1> 
        <s1><s3><u>overlap</u></s3></s1> 
        <s3>viral</s3>
        <s2>I am black</s2>
        """
        H1=int(H_start)-offset
        H2=int(H_stop)-offset
        V1=int(V_start)-offset
        V2=int(V_stop)-offset
        O1=int(O_start)-offset
        O2=int(O_stop)-offset
        data=data.replace("-", "_") # some viewers may treat '-' as a special character
        info=data[offset:]
        result=""
        if (O1 == O2): # no overlap -------------
            if H1 < V1:  #host first
                #      normal           host                normal             viral               normal
                result=info[:H1]+"<s1>"+info[H1:H2]+"</s1>"+info[H2:V1]+"<s3>"+info[V1:V2]+"</s3>"+info[V2:]
            else: #viral first
                #      normal           viral               normal             host                normal
                result=info[:V1]+"<s3>"+info[V1:V2]+"</s3>"+info[V2:H1]+"<s1>"+info[H1:H2]+"</s1>"+info[H2:]
        else: # over lap found --------------------
            if H1 < V1:  #host first with overlap
                #      normal           host                              overlap                             viral               normal
                result=info[:H1]+"<s1>"+info[H1:V1]+"</s1>"+"<s1><s3><u>"+info[V1:H2]+"</u></s3></s1>"+"<s3>"+info[H2:V2]+"</s3>"+info[V2:]
            else: # viral first with overlap
                #      normal           viral                             overlap                             host                normal
                result=info[:V1]+"<s3>"+info[V1:H1]+"</s3>"+"<s1><s3><u>"+info[H1:V2]+"</u></s3></s1>"+"<s1>"+info[V2:H2]+"</s1>"+info[H2:]
        result=result.replace(" ", "&nbsp;") # replace blank with HTML special blanks, otherwie only one blank available between words
        #result=result.replace("-", "&#45;") # some viewers may treat '-' as a special character
        return result
  

    def getGetSequenceFile(self, filename, name):
        label='index   row#      Sequence: Red=Viral Highlighted=Host Underlined=Overlap'
        sl=len(label)
        for i in range(sl, self.maxSeqLength):
            label=label+' '
        label1="|   Chromosome        |      HostRefCords      | Host Orientation |      ViralRefCords     | ViralOrientation |        Gene       |   DistToGene    "
        seq_Label=label+label1
        self.seqLabelLength=len(seq_Label)
        self.seqLabel=seq_Label

        self.name = filename[:-4]
        self.filename = name
        temp=filename
        if len(temp)>0:
            # self.changeLabel(self.fileLabel,temp)
            self.loadFile(temp)
        else:
            print("no file selected")

        return self.filename

    def loadFile(self, filename):
        self.rows=[]
        self.list_data=[]
        self.display_list_data=[]
        self.maxLength=0
        #0-9:  'ReadName','ReadLength','Sequence','Chromosome','Gene','InsideGene','DistanceToGene','GeneDirection','Focus','GeneObj',
        #10-19:'FocusObj','HostLocalCords','Hlength','HostRefCords','HostOrientation','ViralAccession','ViralLocalCords','Vlength','ViralRefCords','ViralOrientation',
        #20-29:'HTM','HTMAdjusted','VTM','VTMAdjusted','Overlap','OverlapTM','OverlapTMAdjusted','Inserted','Microhomology','HostMapFlag',
        #30-39:'ViralMapFlag','HMapQ','VMapQ','FastqSource','Index']
        """
        try:
        """
        with open (filename,'r') as csvfile:
            reader = csv.reader(csvfile, delimiter=',')
            read_index = 0
            ndx=-1
            for row in reader:
                self.rows.append(row)
                ReadLength=row[1]
                seq=row[2]
                Hlength=row[12]
                HostLocalCords=row[11]
                HostRefCords=row[13]
                HostOrientation=row[14]
                Vlength=row[17]
                ViralLocalCords=row[16]
                ViralRefCords=row[18]
                ViralOrientation=row[19]
                Overlap=row[24]
                Chromosome=row[3]
                Gene=row[4]
                ndx+=1
                if (ndx==0):
                    continue #skip column title
                
                # distance data was saved in float string format or "N/A", so convert it into float and then int
                # use exceptto catch "N/A"
                try:
                    DistanceToGene=int(float(row[6]))
                except:
                    DistanceToGene=self.NA_value
                
                # cords format [xxxxx,xxxxx]
                # orientation format xxxxxx/xxxxx
                H_ref_start, H_ref_stop = HostRefCords.replace("[", "").replace("]", "").replace(" ", "").split(",")
                H_orientation_1, H_orientation_2 = HostOrientation.replace(" ","").split("/")
                H_orientation_2 = H_orientation_1[0]+"_"+H_orientation_2[0]
                V_ref_start, V_ref_stop = ViralRefCords.replace("[", "").replace("]", "").replace(" ", "").split(",")
                V_orientation_1, V_orientation_2 = ViralOrientation.replace(" ","").split("/") 
                V_orientation_2 = V_orientation_1[0]+"_"+V_orientation_2[0]
                if int(ReadLength) > self.maxLength:
                    self.maxLength = int(ReadLength)
                #print(ReadLength+" "+HostLocalCords+" "+ViralLocalCords)
                # cords format [xxxxx,xxxxx]
                H_start, H_stop=HostLocalCords.replace("[", "").replace("]", "").replace(" ", "").split(",")
                H_start = str(int(H_start)-1)
                H_stop = str(int(H_stop))
                
                V_start, V_stop=ViralLocalCords.replace("[", "").replace("]", "").replace(" ", "").split(",")
                V_start = str(int(V_start)-1)
                V_stop = str(int(V_stop))

                overlap=int(Overlap)
                if (overlap > 0):
                    v_start = int(V_start)
                    h_start = int(H_start)
                    if (v_start > h_start):
                        O_start = V_start
                        O_stop = str(v_start + overlap)
                    else:
                        O_start = H_start
                        O_stop = str(h_start + overlap)
                else:
                    O_start = O_stop = str(0)
                    
                disp_seq=seq
                H_seq=disp_seq[int(H_start):int(H_stop)]
                V_seq=disp_seq[int(V_start):int(V_stop)]
                if (int(V_start) <= int(H_start)):
                    H_first = False
                else:
                    H_first = True
                """    
                # check for flip over -----------------------------
                if (H_first == False):
                    read_length = int(ReadLength)
                    o_start = int(O_start)
                    o_stop = int(O_stop)
                    if (o_start != 0) or (o_stop != 0):
                        o_start = read_length - o_start
                        o_stop = read_length - o_stop
                        O_start = str(o_stop)
                        O_stop = str(o_start)
                        
                    h_start = int(H_start)
                    h_stop = int(H_stop)
                    h_start = read_length - h_start
                    h_stop = read_length - h_stop
                    H_start = str(h_stop)
                    H_stop = str(h_start)
                    
                    v_start = int(V_start)
                    v_stop = int(V_stop)
                    v_stop = read_length - v_stop
                    v_start = read_length - v_start
                    V_start = str(v_stop)
                    V_stop = str(v_start)

                    sequence_complement=str(Seq(seq).reverse_complement())
                    disp_seq = sequence_complement
                    #print(H_start, H_stop, V_start, V_stop, O_start, O_stop)
                """    
                temp=Chromosome[:2]
                temp=temp.replace("_", "")
                try:
                    Chromosome_num=int(temp)
                except:
                    temp=temp.upper()
                    if (temp==""):
                        Chromosome_num=99                        
                    elif (temp[0]=='M'):
                        Chromosome_num=23
                    elif (temp[0]=='X'):
                        Chromosome_num=24
                    elif (temp[0]=='Y'):
                        Chromosome_num=25
                    else:
                        Chromosome_num=26
                        
                self.list_data.append([ndx+1, disp_seq, len(disp_seq), H_seq,int(H_stop)-int(H_start),V_seq,int(V_stop)-int(V_start),
                                  int(H_start),int(H_stop),int(V_start),int(V_stop), int(Overlap), int(O_start), int(O_stop),
                                  H_first,str(Chromosome), int(Chromosome_num), int(H_ref_start), int(H_ref_stop),
                                  str(H_orientation_2), int(V_ref_start), int(V_ref_stop),str(V_orientation_2), str(Gene),
                                  int(DistanceToGene)])
                #print(seq)
            
        #self.seqBox.configure(state=DISABLED)
        """    
        except:
            print("Error: can not open file: "+filename)
        """

        self.saveToHTML(self.list_data, 0)
        self.alignSort()
        # self.displaySequence(self.list_data)
        # self.reset_Sort_Button_Color()
        return
    
# if __name__ == '__main__':        
#     root = tk.Tk()
#     app = Application(master=root)
#     app.mainloop()

