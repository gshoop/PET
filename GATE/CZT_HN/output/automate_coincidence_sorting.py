import os
import sys

cwd = os.getcwd()  # Get the current working directory (cwd)
files = os.listdir(cwd)  # Get all the files in that directory
file_list = []
for file in files:
    if file.endswith('Hits.dat'):
        file_list.append(file)

for file in file_list:
    command1 = './main ' + file + ' non bin2 non'
    os.system(command1)

files = os.listdir(cwd)  # Get all the files in that directory
file_list = []
for file in files:
    if file.endswith('Hits.dat.lst'):
        file_list.append(file)

for file in file_list:
    command2 = 'python3 reorgcpp.py ' + file
    os.system(command2)
