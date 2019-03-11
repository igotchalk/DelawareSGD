# -*- coding: utf-8 -*-
"""
Created on Wed Jan 25 11:59:41 2017

@author: ianpg

Write CSV file to tab delimited
"""
import os
os.chdir('E:\SGEMS\Python_SGEMS')
input_file = open('testfile.csv', 'r')
output_file = open('testfile_spacedelimited.txt', 'w')
input_file.readline() # skip first line 
for line in input_file:
    (a, date, time, lon, lat, country) = line.strip().split(',')
    output_file.write('\t'.join([lon, lat, country, date, time]) + '\n')
input_file.close()
output_file.close()