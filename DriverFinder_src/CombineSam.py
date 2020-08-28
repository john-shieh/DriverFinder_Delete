import os
import os.path
import sys
import shutil
import time

def combine_sam_files(parent, result_full_path, sample_folder_full_path, sample_folder):
    print("======= combine sam files for sample: "+sample_folder_full_path)
    combine_bam_file = result_full_path+"/"+"[combine]"+sample_folder+".bam"
    #remove existing file
    if os.path.isfile(combine_bam_file):
        os.remove(combine_bam_file)
    split_folders = os.listdir(sample_folder_full_path)
    first_file = True
    for name in split_folders:
        split_sample_full_name = sample_folder_full_path+"/"+name
        if os.path.isfile(split_sample_full_name): # it is a file
            #skip file
            continue
        else: # it is a sample file split folder
            files = os.listdir(split_sample_full_name)
            for file in files:
                full_name = split_sample_full_name+"/"+file
                full_name = full_name.replace("/", "\\")
                if "[" in full_name:
                    continue #skip if file includes "[" as part of file name
                if "viralAlign.sam" in full_name:
                    print("   ---- process viral sam file:"+full_name)

                    """
                    """
                    # convert to bam format ----------------------------------
                    # the created bam file is named [new]file
                    sam_file=full_name
                    bam_file=split_sample_full_name+"/"+"[new]"+file.replace(".sam", ".bam")
                    bam_file=bam_file.replace("/", "\\")
                    command0 = "samtools view -bS "+sam_file+ " > " + bam_file
                    print("   create bam:       "+command0)
                    os.system(command0)

                    # remove ummapped reads
                    # the created mapped bam file is anmed [map]file
                    map_bam_file=bam_file.replace("new", "map")
                    command1 = "samtools view -b -F 4 "+bam_file+" > "+ map_bam_file
                    print("   remove unmapped:  "+command1)
                    os.system(command1)
                    os.remove(bam_file) #remove bam file

                    # sort reads from bam file -------------------------------
                    # the created sorted bam file is named [sorted]file
                    sort_bam_file=map_bam_file.replace("map", "sorted").replace(".bam", "")
                    command2 = "samtools sort " + map_bam_file + " " + sort_bam_file
                    print("   sort:             "+command2)
                    os.system(command2)
                    os.remove(map_bam_file) #remove map bam file
                    
                    # remove duplicates from bam file ------------------------
                    # the created bam file with duplicate removed is named [rmdup]file
                    sort_bam_file=sort_bam_file+".bam"
                    rmdup_bam_file=sort_bam_file.replace("sorted", "rmdup")
                    command3 = "samtools rmdup -S " + sort_bam_file + " " + rmdup_bam_file
                    print("   remove duplicate: "+command3)
                    os.system(command3)
                    os.remove(sort_bam_file) #remove temp sort bam file

                    # combine bam files
                    if first_file == True:
                        first_file = False
                        # no need to merge first file
                        shutil.copy(rmdup_bam_file, combine_bam_file)
                        print("   copy 1st bam file to combine bam")
                    else:
                        tmp_combine_bam_file = combine_bam_file.replace("combine", "tmp")
                        command4 = "samtools merge -f " + tmp_combine_bam_file + " " + combine_bam_file + " "+rmdup_bam_file
                        print("   merge:            "+command4)
                        os.system(command4)
                        shutil.copy(tmp_combine_bam_file, combine_bam_file)
                        os.remove(tmp_combine_bam_file)

                    # debug purpose
                    tmp_sam_file=rmdup_bam_file.replace(".bam", ".sam")
                    command5 = "samtools view -h " + rmdup_bam_file + " > " + tmp_sam_file
                    os.system(command5)
                    
                    os.remove(rmdup_bam_file) #remove temp rmdup bam file

    # convert the final bam to sam format ----
    # sort combine bam
    sort_combine_bam_file=combine_bam_file.replace("combine", "sorted").replace(".bam", "") 
    command6 = "samtools sort " + combine_bam_file + " " + sort_combine_bam_file
    print("   sort combine      "+command6)
    os.system(command6)

    # remove duplicates again from combine bam file
    sort_combine_bam_file = sort_combine_bam_file+".bam"
    rmdup_combine_bam_file=sort_combine_bam_file.replace("sorted", "rmdup")
    command7 = "samtools rmdup -S " + sort_combine_bam_file + " " + rmdup_combine_bam_file
    print("   remove duplicate: "+command7)
    os.system(command7)
    os.remove(sort_combine_bam_file) #remove temp sort bam file
    
    combine_sam_file=combine_bam_file.replace(".bam", ".sam")
    combine_sam_file = combine_sam_file.replace("[combine]", "")
    combine_sam_file = combine_sam_file.replace(".sam", "__viralAlign.sam")
    command8 = "samtools view -h " + rmdup_combine_bam_file + " > " + combine_sam_file
    print(" create final sam: "+command8)
    os.system(command8)
    os.remove(rmdup_combine_bam_file) #remove temp file
    os.remove(combine_bam_file) #remove temp file
    print("   combine sam file: "+combine_sam_file)
    parent.logInfo("Created the combined SAM file: " + combine_sam_file)
    if (parent != None):
        parent.enableOptions()
    
def process(parent, filesFolder=None):

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

    #start = timer()

    logInfo("Starting the combination and de-duplication process.")
    if (filesFolder == None):
        print ('Please provide folder name.')
        return False
    else:
        # get all the folders/files 
        folders = os.listdir(filesFolder)
        for item_name in folders:
            full_path_name = filesFolder+"/"+item_name
            if os.path.isfile(full_path_name): # it is a file
                # skip non sample file folder
                continue
            else: # it is a sample file split folder
                logInfo("Processing " + item_name)
                combine_sam_files(parent, filesFolder, full_path_name, item_name)

# this will not be called if import as a module
if __name__ == '__main__':
    if not os.path.isfile("samtools.exe"):
        print("samtools not installed.")
    else:
        folder_name = "E:/SingleSample"
        process(folder_name)

