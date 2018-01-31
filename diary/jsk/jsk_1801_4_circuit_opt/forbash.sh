#!/bin/bash

for i in {1..667}
do
	python main.py
done

dumpfile="netnodump.txt"
echo '1' > $dumpfile
