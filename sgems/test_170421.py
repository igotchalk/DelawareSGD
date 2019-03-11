# -*- coding: utf-8 -*-
"""
Created on Fri Apr 21 23:56:15 2017
Renaming script because I forgot to put the underscore in the name 
"""
ucn_dir = 'C:\\Users\\ianpg\\Dropbox\\Aurora_ARR\\MODFLOW\\FlopyModel\\Flopy_NWT\\All_concs\\'
os.chdir(ucn_dir)
all_files = os.listdir(ucn_dir)

#SELECTION CRITERIA
files = []
for i in range(len(all_files)):
    if len(all_files[i].split('_')) ==4:
        if all_files[i].split('_')[3][-8:-4] == '1498':
            files.append(all_files[i])
            

for file in files:
    splitup = file.split('_')
    if not len(splitup) == 4:
        break
    else:
        smashlength =len(splitup[3]) #length of suffix should be ==18
        simnum = splitup[3][0:smashlength-18]
        newname = splitup[0] + '_' + splitup[1] + '_' + splitup[2] + '_' + str(simnum) + '_' + splitup[3][smashlength-18:]
        print(newname) 
        os.rename(file,newname)


