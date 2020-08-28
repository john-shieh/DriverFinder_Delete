
def Proc(skip_pysam=False):
    print("check Python version")
    result = pythonVerChk()
    if result == False:
        return result
    result = checkRequiredModules(skip_pysam)
    if result == False:
        return result
    return True

# check Python version
def pythonVerChk():
    status = True
    import sys
    if sys.hexversion < 0x3050000:
        print("Error: Python 3.5 or greater needed")
        print("Please install Python 3.5 or later")
        print("Program STOP!!!")
        status = False
    return status

# check required modules
def module_version_check(version, need_version):
    status = True
    need_info = need_version.split(".")
    major = int(need_info[0])
    minor = int(need_info[1])
    ver_info=version.split(".")
    ver_major = int(ver_info[0])
    ver_minor = int(ver_info[1])
    if ver_major < major:
        status = False
    if ver_major > major:
        status = True
    if ver_major == major:
        if ver_minor >= minor:
            status = True
        else:
            status = False
    return status

def checkRequiredModules(skip_pysam):
    modulesStatus = True
    module_error=''
    version_error=''
    try:
        import numpy
        need_ver="1.13.x"
        print("numpy module, version = ", numpy.__version__)
        if module_version_check(numpy.__version__, need_ver) == False:
            print("  >>>Please update numpy module to version "+need_ver+".")
            version_error=version_error + "numpy "
    except:
        module_error = module_error + "numpy "

    try:
        import pandas
        need_ver="0.20.x"
        print("pandas module, version = ", pandas.__version__)
        if module_version_check(pandas.__version__, need_ver) == False:
            print("  >>>Please update pandas module to version "+need_ver+".")
            version_error=version_error + "pandas "
    except:
        module_error = module_error + "pandas "
        
    try:
        import Bio # for Biopython module
        need_ver = "1.67.x"
        print("Biopython module, version = ", Bio.__version__)
        if module_version_check(Bio.__version__, need_ver) == False:
            print("  >>>Please update biopython module to version "+need_ver+".")
            version_error=version_error + "biopython "
    except:
        module_error = module_error + "Biopython "
        
    try:
        import xlrd # for openpyxl
        need_ver="1.0.x"
        print("xlrd module, version = ", xlrd.__VERSION__)
        if module_version_check(xlrd.__VERSION__, need_ver) == False:
            print("  >>>Please update xlrd module to version "+need_ver+".")
            version_error=version_error + "xlrd "
    except:
        module_error = module_error + "xlrd "

    try:
        import openpyxl
        need_ver="2.4.x"
        print("openpyxl module, version = ", openpyxl.__version__)
        if module_version_check(openpyxl.__version__, need_ver) == False:
            print("  >>>Please update openpyxl module to version "+need_ver+".")
            version_error=version_error + "openpyxl "
    except:
        module_error = module_error + "openpyxl "

    try:
        import cython # for pysam
        need_ver="0.24.x"
        print("cython module, version = ", cython.__version__)
        if module_version_check(cython.__version__, need_ver) == False:
            print("  >>>Please update cython module to version "+need_ver+".")
            version_error=version_error + "cython "
    except:
        module_error = module_error + "cython "

    if (skip_pysam != True):
        try:
            import pysam
            need_ver="0.11.x"
            print("pysam module, version = ", pysam.__version__)
            if module_version_check(pysam.__version__, need_ver) == False:
                print("  >>>Please update pysam module to version "+need_ver+".")
                version_error=version_error + "pysam "
        except:
            module_error = module_error + "pysam"


    if module_error != '':
        print("\n--------------------------------------------------------------")
        print("modules not installed: "+module_error)
        print("Please install required modules first.")
        modules=module_error.split(" ")
        if (len(modules) == 1) and (modules[0]=="pysam"):
            print("\nNote: The 'Detect HBV Mutations' will be disabled if the pysam module is not installed.")
        print("--------------------------------------------------------------\n")

    if version_error != '':
        print("\n--------------------------------------------------------------")
        print("modules update needed: "+version_error)
        print("Please update those modules.")
        print("--------------------------------------------------------------\n")


    if (module_error != '') or (version_error != ''):
        modulesStatus = False
    else:
        print("All required modules loaded!!!")

    return modulesStatus 
#===========================================
testModulesStatus = Proc(skip_pysam=True)
print("test required modules = ", testModulesStatus)





