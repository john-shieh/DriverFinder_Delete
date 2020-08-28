# -*- coding: utf-8 -*-

import csv, os
import platform
from tkinter import filedialog,ttk
from tkinter import messagebox
import tkinter as tk
from Bio.Seq import Seq
from Bio import SeqIO

from tkinter import font as tkfont

class ViewAlignedSequence():
    platformType=platform.system()
    list_data_viral=[]
    app_version="1.1"
    app_title = "Viral Aligned Viewer"
    NA_value=999999999
    seqBoxFontSizeDefault = 8
    seqBoxFontSize = seqBoxFontSizeDefault
    file_aligned = ""

    helpWindow = 0
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
    
    def __init__(self, parent):
        verion=" 1.12"
        self.title = self.app_title+" "+self.app_version
        self.p=parent
        
    def quit_Application(self):
        root.destroy()
        return 1

    def logInfo(self, message):
        self.p.logInfo(message)
        return
    
    def alignSort2(self):
        #[0]: index (row # in CSV file)
        #[1]: sequence
        #[2]: sequence length
        #[3]: host sequence
        #[4]: host length
        #[5]: viral sequence
        #[6]: viral length
        #[7]: host start
        #[8]: host stop
        #[9]: viral start
        #[10]: viral stop
        #[11]: overlap length
        #[12]: overlap start
        #[13]: overlap stop
        #[14]: host first
        #[15]: Chromosome
        #[16]  
        #[17]: host ref cords start
        #[18]: host ref cords stop
        #[19]: host orientation
        #[20]: viral ref cords start
        #[21]: viral ref cords stop
        #[22]: viral orientation
        #[23]: gene
        #[24]: distance to gene
        #sortItemsOptionList=(NONE_ITEM, VIRAL, HOST, ENTIRE_SEQ, CHROMOSOME, HOST_REF_COORDS_START, HOST_REF_COORDS_STOP,
        #                HOST_ORIENTATION, VIRAL_REF_COORDS_START, VIRAL_REF_COORDS_STOP, VIRAL_ORIENTATION,
        #                GENE, DISTANCE_TO_GENE, VIRAL_FIRST, HOST_FIRST)
        self.logInfo("Sorting...")

        list_data=self.list_data
        display_type=0
        item=None
        item=self.HOST_REF_COORDS_STOP
        if item==self.HOST_REF_COORDS_STOP:
            list_data=sorted(list_data, key=lambda item:item[18]) 
        item=self.VIRAL_REF_COORDS_START
        if item==self.VIRAL_REF_COORDS_START:
            list_data=sorted(list_data, key=lambda item:item[20]) 
        if item==self.CHROMOSOME:
            list_data=sorted(list_data, key=lambda item:item[15]) 
            list_data=sorted(list_data, key=lambda item:item[16])
        #item=self.GENE
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
                
        #self.displayAlignedSequence2(list_data, display_type)
        self.list_data_viral = list_data
        self.displayAlignedViral(list_data, display_type)
        self.logInfo("Sorting complete!")
        return

    def setSeqBoxSize(self, size):
        w, h = self.WH[size]
        self.seqBox_v.configure(width=w, height=h)
        self.seqBox_v.update()
        self.labelText_v_pos.configure(width=w)
        self.labelText_v_pos.update()
        
    def increaseSeqBox_V_FontSize(self, event):
        if self.seqBoxFontSize < self.seqBoxFontSizeDefault:
            self.seqBoxFontSize += 1
            disp_font2=tkfont.Font(family="Monaco", size=self.seqBoxFontSize)
            self.seqBox_v.configure(font = disp_font2)
            self.labelText_v_pos.configure(font = disp_font2)
            self.setSeqBoxSize(self.seqBoxFontSize)
            #self.displayAlignedViralSequence(list_data)
        return

    def reduceSeqBox_V_FontSize(self, event):
        if self.seqBoxFontSize > 1:
            self.seqBoxFontSize -= 1
            disp_font2=tkfont.Font(family="Monaco", size=self.seqBoxFontSize)
            self.seqBox_v.configure(font = disp_font2)
            self.labelText_v_pos.configure(font = disp_font2)
            self.setSeqBoxSize(self.seqBoxFontSize)
            #self.displayAlignedViralSequence(list_data)
        return

    def cancelSeqBox_V_FontSize(self, event):
        self.seqBoxFontSize = self.seqBoxFontSizeDefault
        disp_font2=tkfont.Font(family="Monaco", size=self.seqBoxFontSize)
        self.seqBox_v.configure(font = disp_font2)
        self.labelText_v_pos.configure(font = disp_font2)
        self.setSeqBoxSize(self.seqBoxFontSize)
        #self.displayAlignedViralSequence(list_data)
        return


    def xview_v(self, *args):
        self.labelText_v_pos.xview(*args)
        self.seqBox_v.xview(*args)
        return

    def showChimericReads(self, event):
        list_data = self.list_data_viral
        self.displayAlignedViralSequence(list_data, 1)
        return
    
    def showViralReads(self, event):
        list_data = self.list_data_viral
        self.displayAlignedViralSequence(list_data, 2)
        return
    
    def showMajorJuctionReads(self, event):
        list_data = self.list_data_viral
        self.displayAlignedViralSequence(list_data, 3)
        return

    def displayAlignedViral(self, list_data, sorting_type=0):
        self.WH=[[200,40],[1390,175],[695,175],[695,135],[465,90],[348,72],[280,70],[280,45],[200,40]]
        self.alignViralWindow = tk.Toplevel()
        self.alignViralWindow.wm_attributes("-topmost", 1)
        self.alignViralWindow.title("Align Viral - "+self.file_aligned)
        disp_font=tkfont.Font(family="Monaco", size=self.seqBoxFontSize)
        self.disp_font = disp_font
        self.alignViralWindow_frame1 = tk.Frame(self.alignViralWindow)
        self.alignViralWindow_frame1.grid(row=0,column=0,sticky='nw',columnspan=1)
        #label
        self.labelText_v= tk.Text(self.alignViralWindow_frame1,wrap=tk.NONE, height=1,width=200,relief=tk.SUNKEN, font = disp_font)
        self.labelText_v.tag_configure("Virus", foreground="red")
        self.labelText_v.tag_configure('Host',background='yellow')
        self.labelText_v.tag_configure('Overlap',underline=True)
        seq_Label_v='                  Sequence: Red=Viral Highlighted=Host Underlined=Overlap'

        self.changeLabel(self.labelText_v,seq_Label_v)
        self.labelText_v.tag_add('Virus','1.28','1.31')
        self.labelText_v.tag_add('Host','1.38','1.49')
        self.labelText_v.tag_add('Virus','1.55','1.65')
        self.labelText_v.tag_add('Host','1.55','1.65')
        self.labelText_v.tag_add('Overlap','1.55','1.65')
        self.labelText_v.configure(state=tk.DISABLED)
        self.labelText_v.grid(row=0,column=0,sticky='nw')
        self.seqBox_v_width=3404
        self.seqBox_v_width_max=3510
        #self.seqBoxFontSize = 1
        disp_font=tkfont.Font(family="Monaco", size=self.seqBoxFontSize)
        self.labelText_v_pos= tk.Text(self.alignViralWindow_frame1,wrap=tk.NONE, height=1,width=200,relief=tk.SUNKEN, font = disp_font)
        self.labelText_v_pos.grid(row=1,column=0,sticky='nw')
        disp_info=""
        for i in range(self.seqBox_v_width):
            disp_info=disp_info+"."
        for i in range(0, self.seqBox_v_width, 50):
            num=str(i)
            disp_info=disp_info[:i]+num+disp_info[i+len(num):]
        for i in range(self.seqBox_v_width, self.seqBox_v_width_max):
            disp_info="."+disp_info
        self.changeLabel(self.labelText_v_pos,disp_info)

        #seq box
        self.seqBox_v=tk.Text(self.alignViralWindow_frame1,wrap=tk.NONE, height=40,width=200,relief=tk.SUNKEN, font = disp_font)
        #self.seqBox_v=tk.Text(self.alignViralWindow_frame1,wrap=tk.NONE, height=40,width=150,state=tk.DISABLED,relief=tk.SUNKEN, font = disp_font)
        self.seqBox_v.grid(row=2,column=0,sticky='nw')
        self.seqBox_v.bind('<F6>',self.increaseSeqBox_V_FontSize)
        self.seqBox_v.bind('<F5>',self.reduceSeqBox_V_FontSize)
        #self.seqBox_v.bind('<Alt-Right>',self.cancelSeqBox_V_FontSize)
        self.seqBox_v.bind('<F1>',self.showChimericReads)
        self.seqBox_v.bind('<F2>',self.showViralReads)
        self.seqBox_v.bind('<F3>',self.showMajorJuctionReads)
        self.alignViralWindow.bind('<F6>',self.increaseSeqBox_V_FontSize)
        self.alignViralWindow.bind('<F5>',self.reduceSeqBox_V_FontSize)
        #self.alignViralWindow.bind('<Alt-Right>',self.cancelSeqBox_V_FontSize)
        self.alignViralWindow.bind('<F1>',self.showChimericReads)
        self.alignViralWindow.bind('<F2>',self.showViralReads)
        self.alignViralWindow.bind('<F3>',self.showMajorJuctionReads)
        self.seqBox_v.tag_configure("Virus", foreground="red")
        self.seqBox_v.tag_configure('Host',background='yellow')
        self.seqBox_v.tag_configure('Overlap',underline=True)
        self.scrollbar1_v=ttk.Scrollbar(self.alignViralWindow_frame1,command=self.seqBox_v.yview)
        self.scrollbar2_v=ttk.Scrollbar(self.alignViralWindow_frame1,orient=tk.HORIZONTAL)
        self.seqBox_v.config(xscrollcommand=self.scrollbar2_v.set)
        self.seqBox_v.config(yscrollcommand=self.scrollbar1_v.set)
        self.scrollbar2_v.configure(command=self.xview_v)
        self.scrollbar2_v.grid(row=3,column=0,sticky='ew')
        self.scrollbar1_v.grid(row=2,column=1,sticky='ns')
        
        #self.customSortWindow.protocol('WM_DELETE_WINDOW', self.customSort_Close)
        self.setSeqBoxSize(self.seqBoxFontSize)
        self.displayAlignedViralSequence(list_data)
        
        return

    def displayAlignedViralSequence(self, list_data, disp_type = 1):
        self.seqBox_v.configure(state=tk.NORMAL)
        self.seqBox_v.delete(1.0, tk.END)
        #[0]: index (row # in CSV file)
        #[1]: sequence
        #[2]: sequence length
        #[3]: host sequence
        #[4]: host length
        #[5]: viral sequence
        #[6]: viral length
        #[7]: host start
        #[8]: host stop
        #[9]: viral start
        #[10]: viral stop
        #[11]: overlap length
        #[12]: overlap start
        #[13]: overlap stop
        #[14]: host first
        #[15]: Chromosome
        #[16]  
        #[17]: host ref cords start
        #[18]: host ref cords stop
        #[19]: host orientation
        #[20]: viral ref cords start
        #[21]: viral ref cords stop
        #[22]: viral orientation
        #[23]: gene
        #[24]: distance to gene
        seqBox_v_start_offset = self.seqBox_v_width_max - self.seqBox_v_width
        blank_line=""
        for i in range(self.seqBox_v_width_max):
            blank_line=blank_line+" "
        disp_info=""
        blanks=""
        for i in range(10):
            blanks=blanks+"                    "
        non_mapped_extra=0
        line=0
        v_start_min=0
        prev_host_ref_stop=-1
        prev_viral_ref_start=-1
        next_host_ref_stop=-1
        next_viral_ref_start=-1
        for i in range(len(list_data)):
            junction_found = False
            curr_host_ref_stop = int(list_data[i][18])
            curr_viral_ref_start = int(list_data[i][20])
            if i < (len(list_data) - 1):
                next_host_ref_stop = int(list_data[i+1][18])
                next_viral_ref_start = int(list_data[i+1][20])
            else:
                next_host_ref_stop=-1
                next_viral_ref_start=-1
            line+=1
            whole_seq=list_data[i][1]
            host_seq=list_data[i][3]
            viral_seq=list_data[i][5]
            #print(viral_seq)

            disp_seq=list_data[i][5]
            ndx=line
            o_start = int(list_data[i][12])
            o_stop = int(list_data[i][13])
            h_start = int(list_data[i][7])
            h_stop = int(list_data[i][8])
            v_start = int(list_data[i][9])
            v_stop = int(list_data[i][10])
            # check for flip over -----------------------------
            if ((curr_host_ref_stop == prev_host_ref_stop) and (curr_viral_ref_start == prev_viral_ref_start)) or \
               ((curr_host_ref_stop == next_host_ref_stop) and (curr_viral_ref_start == next_viral_ref_start)):
                #print("juction found ndx = ", i)
                junction_found = True
                if list_data[i][14] == False:  #not host first
                    read_length = int(list_data[i][2])
                    if (o_start != 0) or (o_stop != 0):
                        o_start_new = read_length - o_stop
                        o_stop_new = read_length - o_start
                        o_start=o_start_new
                        o_stop=o_stop_new
                    h_start_new = read_length - h_stop
                    h_stop_new = read_length - h_start               
                    h_start=h_start_new
                    h_stop=h_stop_new
                    v_start_new = read_length - v_stop
                    v_stop_new = read_length - v_start
                    v_start=v_start_new
                    v_stop=v_stop_new

                    sequence_complement=str(Seq(whole_seq).reverse_complement())
                    whole_seq = sequence_complement
                    sequence_complement=str(Seq(viral_seq).reverse_complement())
                    viral_seq = sequence_complement

            # if a sequence is too long, it will be truncated.
            # so, we have to make sure that color will be truncated as well
            rec_num=str(ndx)+"      "
            rec_num=rec_num[:6]+"["
            rec_num1=str(list_data[i][0])+"      "
            rec_num=rec_num+rec_num1[:6]+"]_"
            rec_num=rec_num+" "
            disp_info=blank_line
            disp_info=rec_num+disp_info[len(rec_num):]

            #[20]: viral ref cords start
            #[21]: viral ref cords stop
            #disp_info=temp_info+disp_seq
            v_r_start=int(list_data[i][20]) + seqBox_v_start_offset
            v_r_stop=int(list_data[i][21]) + seqBox_v_start_offset
            if (disp_type == 2):
                disp_info=disp_info[:v_r_start]+viral_seq+disp_info[v_r_stop:]

            display_viral=True
            if (disp_type == 1) or (disp_type == 3):
                if (disp_type == 3) and (not junction_found):
                    display_viral=False
                seq_offset = 0 - v_start
                seq_start=v_r_start + seq_offset
                seq_stop=seq_start+(int(list_data[i][2]))
                if display_viral:
                    disp_info=disp_info[:seq_start]+whole_seq+disp_info[seq_stop:]
                offset = h_start - v_start
                offset2 = o_start - v_start
                h_start=v_r_start + offset
                h_stop=h_start+(int(list_data[i][4]))
                o_start=v_r_start + offset2
                o_stop=o_start+(int(list_data[i][11]))
                #print(ndx, " --- ", h_start, h_stop, v_start, v_stop, o_start, o_stop)
                #print(ndx, "     ", list_data[i][7], list_data[i][8],list_data[i][9], list_data[i][10])
            # add info into display 
            self.seqBox_v.insert(tk.END, disp_info+"\n")
            if (disp_type == 1) or (disp_type == 3):
                H_start = str(h_start)
                H_stop = str(h_stop)
                if display_viral:
                    self.seqBox_v.tag_add('Host',str(ndx)+"."+H_start,str(ndx)+"."+H_stop)
                    if (o_start != o_stop):
                        O_start = str(o_start)
                        O_stop = str(o_stop)
                        self.seqBox_v.tag_add('Overlap',str(ndx)+"."+O_start,str(ndx)+"."+O_stop)
            V_start = str(v_r_start)
            V_stop = str(v_r_stop+1)
            if display_viral:
                self.seqBox_v.tag_add('Virus',str(ndx)+"."+V_start,str(ndx)+"."+V_stop)

            prev_host_ref_stop = int(list_data[i][18])
            prev_viral_ref_start = int(list_data[i][20])

        # update display
        self.seqBox_v.update()

        return
    
    def changeLabel(self,label,message):
        label.configure(state=tk.NORMAL)
        label.delete(1.0,tk.END)
        label.insert(tk.INSERT,message)
        label.configure(state=tk.DISABLED)
        

    def getGetSequenceFile(self):
        if self.platformType == "Windows":
            temp=filedialog.askopenfilename(title='Select Sequence File',filetypes=(("Sequence Files", "*.csv"),("All files", "*.*")))
        else:
            temp=filedialog.askopenfilename(title='Select Sequence File')
        if len(temp)>0:
            self.loadFile(temp)
        else:
            print("no file selected")
        return

    def loadFile(self, filename):
        self.rows=[]
        self.list_data=[]
        self.display_list_data=[]
        self.maxLength=0
        self.file_aligned = filename
        print("load file: "+filename)
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
                #[0]: index (row # in CSV file)
                #[1]: sequence
                #[2]: sequence length
                #[3]: host sequence
                #[4]: host length
                #[5]: viral sequence
                #[6]: viral length
                #[7]: host start
                #[8]: host stop
                #[9]: viral start
                #[10]: viral stop
                #[11]: overlap length
                #[12]: overlap start
                #[13]: overlap stop
                #[14]: host first
                #[15]: Chromosome
                #[16]  
                #[17]: host ref cords start
                #[18]: host ref cords stop
                #[19]: host orientation
                #[20]: viral ref cords start
                #[21]: viral ref cords stop
                #[22]: viral orientation
                #[23]: gene
                #[24]: distance to gene
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
        return

