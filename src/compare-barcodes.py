#!/usr/bin/env python3

from fuzzywuzzy import fuzz
from fuzzywuzzy import process
from subprocess import Popen, PIPE, check_output
import gzip
import io
from collections import Counter


splitreadsFn = '03_initial_miseq_analysis/data/processed/combined-splitreads.tab.gz'
splitreads = '03_initial_miseq_analysis/data/processed/bcs.txt'

p1 = Popen(['zcat', splitreadsFn], stdout = PIPE)
p2 = Popen(['tail', '-n', '+2'], stdin = p1.stdout, stdout=PIPE)
p3 = Popen(['awk', '-F', ',', '{print $1,$3}'], stdin=p2.stdout, stdout=PIPE)
bcs = []

with open(splitreads, 'r') as f:
    for line in f:
        s = line.split()
        bcs.append([s[0]+s[1]])

with gzip.open(splitreadsFn, 'rt') as f:
    for line in f:
        s = line.split(',')
        bcs.append([s[0],s[2]])


for line in io.TextIOWrapper(p3.stdout, encoding='utf-8'):
    bcs.append(''.join(line.split()))

while True:
    line = p3.stdout.readline()
    if not line:
        break
    print("test:", line.rstrip())

output = [x.strip().split() for x in p3.stdout.readlines()]

'zcat 03_initial_miseq_analysis/data/processed/combined-splitreads.tab.gz | tail -n +2'

