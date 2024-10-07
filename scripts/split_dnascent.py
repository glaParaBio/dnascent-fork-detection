#!/usr/bin/env python3

# Why? See https://stackoverflow.com/questions/14207708/ioerror-errno-32-broken-pipe-python
from signal import signal, SIGPIPE, SIG_DFL
import sys
import argparse
import subprocess
import tempfile
import atexit
import os
from collections import OrderedDict

signal(SIGPIPE, SIG_DFL) 

parser = argparse.ArgumentParser(description= 'Split the output of `DNAscent` subcommands `detect` or `forkSense`')
parser.add_argument('input', help= 'Input file or stdin [%(default)s]', default= '-', nargs= '?')
parser.add_argument('--prefix', '-p', help= 'Output prefix for output files', required= True)
parser.add_argument('--additional-suffix', '-s', help= 'Additional suffix for output files [%(default)s]', default= '.dna')
parser.add_argument('--num-records', '-n', help= 'Number of records per file [%(default)s]', default= 10000, type= int)
parser.add_argument('--version', '-v', action='version', version='%(prog)s 0.1.0')

if __name__ == '__main__':
    args= parser.parse_args()
    
    if args.num_records <= 0:
        raise Exception('Number of records per file must be greater than 0')

    if args.input == '-':
        fin= sys.stdin
    else:
        fin= open(args.input)
    
    fileid = 0
    fname = '%s%04d%s' % (args.prefix, fileid, args.additional_suffix)
    fout = open(fname, 'w')

    n = 0
    tot = 0
    for line in fin:
        if line.startswith('#'):
            continue
        if line.startswith('>'):
            if n == args.num_records:
                n = 0
                fout.close()
                sys.stderr.write('Finished writing %s\r' % fname)
                fileid += 1
                if fileid > 2000:
                    raise Exception('This is a safety error: Too many files in output')
                fname = '%s%04d%s' % (args.prefix, fileid, args.additional_suffix)
                fout = open(fname, 'w')
            n += 1
            tot += 1
        fout.write(line)
    fin.close()
    fout.close()
    sys.stderr.write('%s records written to %s files \n' % (tot, fileid + 1))
