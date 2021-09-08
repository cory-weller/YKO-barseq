#!/usr/bin/env python

import sys

filename = sys.argv[1]

with open(filename, 'r') as infile:
    for line in infile:
        barcode, lineNo = line.strip().split(",")
        print("""%s,%s""" % (barcode[-20:], lineNo))