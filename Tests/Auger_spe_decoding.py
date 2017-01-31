# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a separate function for pulling out binary data from spe files and saving in usable
CSV format.
It has also been incorporated into the larger Auger batch processing project.
"""
import struct, pandas
if 'C:\\Users\\tkc\\Documents\\Python_Scripts' not in sys.path:
    sys.path.append('C:\\Users\\tkc\\Documents\\Python_Scripts')
import Auger_batch_import_functions as Auger
    
# open and process binary spe file data
def processbinarySPE(AugerFileName):      
    AugerFileName='C2010W_18Nov15_105.spe'
    with open(AugerFileName, 'rb') as file:
        filedata = file.read()
    end=filedata.find(b'EOFH')+3  
    headerdata=filedata[0:end+1] # works to cut off binary part (EOFT then \r\n)  
    header=headerdata.decode(encoding='cp437') # more generic encoding than utf-8
    # import binary string
    start=end+3 # this index gives correct binary string (cuts EOFT then \r\n) 
    bindata=filedata[start:] # read in binary part
    
    # Determine # areas, # points
    
    # test output of bytes data
    with open('survey_2areas_binary.spe', 'wb') as mybytes:
        mybytes.write(bindata)
    



# 
testbyte=mybytes[36:44]
s = struct.Struct('I I')
unpacked_bytes = s.unpack(testbyte)
print('Unpacked Values:', unpacked_bytes)

# figuring out hex value associated with known float value (use pack)
testfloat=float(580240)
packed_bytes=s.pack(testfloat)
with open('test_binary.spe', 'wb') as mybytes:
    mybytes.write(packed_bytes)


# function for determining structure of unknown binary strings
# binfile is binary byte object (usually with header string removed)
#startbyte 
mybytes=binfile

startbyte=112
numobjects=1
bytechunk=4
structstring=numobjects*'f'
for i in range(startbyte,startbyte+10):
    testbyte=binfile[i:i+numobjects*bytechunk]
    s = struct.Struct(structstring)
    unpacked_bytes = s.unpack(testbyte)
    print('Unpacked Values:', unpacked_bytes)

def bintestreader(mybytes, startbyte, bytechunk, structstring):
    Binaryparse=pandas.DataFrame(columns=['Startbyte','Byterange','Unpacked_bytes'])
    for i in range(10):
        templist={'Startbyte':'','Byterange':'','Unpacked_bytes':''}
        testbyte=binfile[startbyte+bytechunk*i:startbyte+bytechunk*(i+1)]
        s = struct.Struct(structstring)
        unpacked_bytes = s.unpack(testbyte)
        templist.update({'Startbyte'=startbyte+bytechunk*i,'Byterange'=startbyte+bytechunk*i})
        Binaryparse.loc[i]=pandas.Series(templist)