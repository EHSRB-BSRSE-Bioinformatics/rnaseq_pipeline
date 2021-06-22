# forked from DESeq2Metadata.py, originally © 2018 Adrian Verster

# This script creates the Metadata table that is required by DESeq_Pipe.R
# Uses a regular expression, where group(1) is the technical replicate and group(2) is the dose

import sys, os, re
import pandas as pd
import glob
from optparse import OptionParser

if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-i", "--indir")
    parser.add_option("-o", "--outfile")
    parser.add_option("-r", "--regex_use",default="(M[A-Z]+([0-9]{1})[0-9]{2})") #Of the format MFI105 where the dose is 1, or MB548 where the dose is 5
    parser.add_option("-f", "--first", default = True, action="store_false")
    (options, args) = parser.parse_args()
    indir = options.indir
    outfile = options.outfile
    regex = options.regex_use

    directories = os.listdir(indir)
    print(directories)
    if options.first:
        doses = [re.search(regex, x).group(2) for x in directories]
        animals = [re.search(regex, x).group(1) for x in directories]
    else:
        doses = [re.search(regex, x).group(1) for x in directories]
        animals = [re.search(regex, x).group(2) for x in directories]
    out = pd.DataFrame( {"condition":doses ,"combine":animals} )
    out.index = directories
    out.to_csv(outfile, sep = ",", index = True, index_label = False)