#!/bin/bash

HERE=$(pwd)
cd ..
sudo make clean
sudo make install
cd $HERE
python example_ascii.py

