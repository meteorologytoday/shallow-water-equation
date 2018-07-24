#!/bin/bash

./bin/makefield-gaussian.out

./bin/main.out
python3 plot.py
