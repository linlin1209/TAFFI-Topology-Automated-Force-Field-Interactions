


import sys,os

if len(sys.argv) < 3:
    print "ERROR: file and molecule names must be supplied as the first and second arguments."
    quit()
elif os.path.isfile(sys.argv[1]) is False:
    print "ERROR: could not find file {}".format(sys.argv[1])
else:
    with open(sys.argv[1]) as f:
        for lines in f:
            fields = lines.split()
            if len(fields) >= 2:
                if sys.argv[2] in fields:
                    print " ".join(fields[fields.index(sys.argv[2]) + 1:])
                    

