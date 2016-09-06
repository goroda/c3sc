import os
import sys

pyfiles = [f for f in os.listdir("./") if ((f[-2:] == '.c') or (f[-2:] == '.h')) and (f != 'update_file_headers.py')]

#appendtxt = open('fileheaders.txt', 'r').read()
appendtxt = "#Code\n\n"
print appendtxt
relPath = "./"
for (path, dirs, files) in os.walk(relPath):
    for f in files:
        if ( ((f[-2:] == '.c') or (f[-2:] == '.h')) and (f != "clapack.h") and (f!= "f2c.h")):
            print "modifying file: ", f
            ff = open(os.path.join(path,f), 'r')
            text = ff.read()
            replaceText = appendtxt + text
            #print replaceText
            ff = open(os.path.join(path,f), 'w')
            ff.write(replaceText)

