import sys

with open (sys.argv[1], 'r') as infile:
	for line in infile:
		if line[0] == '>':
			length = line.split('_')[3]
			print(length)