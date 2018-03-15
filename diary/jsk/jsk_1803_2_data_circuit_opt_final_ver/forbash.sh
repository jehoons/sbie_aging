#!/bin/bash

for i in {1..10}
do
	python main.py &
done
wait
