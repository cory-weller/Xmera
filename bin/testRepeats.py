#!/usr/bin/env python

import regex as re
import sys

infileName = sys.argv[1]
stringLength = int(sys.argv[2])
mismatches = int(sys.argv[3])

with open(infileName, 'r') as infile:
    seq = infile.readlines()[1:]
    seq = ''.join([x.strip() for x in seq])


print('\t'.join([str(x) for x in ["stringStart", "stringLength", "substring", "matchStart", "matchEnd", "sub", "ins", "del"]]))
totalMatches = 0
for i in range(len(seq) - stringLength):
    substring = seq[i:(i+stringLength)]
    pattern = re.compile('(%s){e<=%s}' % (substring, mismatches))
    matches = re.finditer(pattern, seq, overlapped = False)
    matches = list(matches)
    N = len(matches)
    if N == 1:
        continue
    else:
        for match in matches:
            start = match.start()
            end = match.end()
            subtituions, insertions, deletions = match.fuzzy_counts
            if abs(start - i) <= insertions:
                continue
            else:
                totalMatches += 1
                print('\t'.join([str(x) for x in [i, substring, start, end, subtituions, insertions, deletions]]))
print('\t'.join([str(x) for x in ["stringStart", "stringLength", "substring", "matchStart", "matchEnd", "sub", "ins", "del"]]))
print('%s total pairs when string length=%s and mismatches<=%s' % (totalMatches, stringLength, mismatches))