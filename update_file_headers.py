import os
import sys

pyfiles = [f for f in os.listdir("./") if ((f[-2:] == '.c') or (f[-2:] == '.h')) and (f != 'update_file_headers.py')]

appendtxt = open('fileheaders.txt', 'r').read()
print appendtxt
# relPath = "./test"
# relPath = "./src"
for (path, dirs, files) in os.walk(relPath):
    for f in files:
        if ( ((f[-2:] == '.c') or (f[-2:] == '.h')) and (f != "clapack.h") and (f!= "f2c.h") and (f != "CuTest.h") and (f != "CuTest.c")):
            print "modifying file: ", f
            ff = open(os.path.join(path,f), 'r')
            text = ff.read().split('//Code')

            replaceText = appendtxt + text[-1]
            print replaceText
            ff = open(os.path.join(path,f), 'w')
            ff.write(replaceText)

