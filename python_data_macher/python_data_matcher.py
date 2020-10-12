#!/usr/bin/python2
# Python Data Matcher 
# Python script to extract data from two spread sheets and produce a single output file.
# Luzie U. Wingen, JIC, 08/11/2018
#from decimal import *
#from sqlobject import *
#from datetime import datetime 
import os
import re
import sys
import string
import copy 
from collections import defaultdict

def get_input(prompt):
  'For compatibility Python 2 and 3'
  if sys.hexversion > 0x03000000:
      return input(prompt)
  else:
      return raw_input(prompt)

def transposed(lists):
   if not lists: return []
   return list(zip(*lists))

def open_file(csv_file): 
  'open file'
  infile = ''
  try : 
    infile = open(csv_file,'r')# start of data upload from file
  except : ### IOError
    print("##                                                                                ##")
    print("##  The file",csv_file,"was not found. Please move it to the current directory.   ##")
    print("##                                                                                ##")
    print("####################################################################################")
  return(infile)

def strip_quotes(item) :
  'strip quotes from list items'
  return re.sub('"','',item)

def choose_file_and_load_data(csv_file,message): 
  'Choose csv file and load spread sheet'
  prompt2 = message
  fname=str.strip(get_input(prompt2))
  if fname=='': # give an example file 
    print("##  No file name was given. An example file will be used:                         ##")
    print("##                                                                                ##")
    fname=csv_file
  print("##  Load data from: "+fname+"")
  infile = open_file(fname)
  csvL=[]
  if infile!='' :
    for line in infile :
      csvline=line.strip().split(',')
      csvL.append(csvline)
  if len(csvL)>0 :
    print("##  "+str(len(csvL))+" rows found in file.")
    print("##                                                                                ##")
###    print("####################################################################################")
  return [fname, csvL] 

def load_csv_and_add(): 
  'load csv files and merge by first column'
  print("")
  csv_file1='data/UK_winter_wheat_RL_2018-19.csv'
  csv_file2='data/UK_winter_wheat_RL_2017-18.csv'
  message = "####################################################################################\n##                                                                                ##\n##  Please type the name of the fist csv' spread sheet:                           ##\n##                                                                                ##\n##  "
  returnL=choose_file_and_load_data(csv_file1,message)
  csv1L=returnL[1]
  csv_file1=returnL[0]
  message = "####################################################################################\n##                                                                                ##\n##  Please type the name of the second 'csv' spread sheet:                        ##\n##                                                                                ##\n##  "
  returnL=choose_file_and_load_data(csv_file2,message)
  csv2L=returnL[1]
  csv_file2=returnL[0]
  if (len(csv1L)>0)&(len(csv2L)>0) :
     merge_tables(csv1L,csv2L,csv_file1)
  else :
    print("##  Please correct the errors and run the script again.                           ##")
    print("####################################################################################")

def merge_tables(csv1L,csv2L,csv_file1):
  'merge the two tables by first column entry'
  print("####################################################################################")
  print("")
  print("    Python Data Matcher will now merge the spread sheets:")
  csv1LT=transposed(csv1L)
  csv2LT=transposed(csv2L)
  featureD=defaultdict()
  for i,n in enumerate(csv2LT[0]) : ### elements in first column
    if (n not in featureD.keys()) and (i!=0) :
      featureD[n]=defaultdict() 
    elif (n not in featureD.keys()) and (i==0) :
      if csv1LT[0][0] not in featureD.keys() :
        featureD[csv1LT[0][0]]=defaultdict() 
      featureD[csv1LT[0][0]]=[csv1LT[0][0]]+csv2L[i]
    elif n in featureD.keys() :
      print("##"+n+"is present twice in the second file. The second item overwrites the first. ##")
    if n in csv1LT[0] : ### csv1L[0] identifier names (first column). Get index of identfier.
      featureD[n]=[csv1LT[0].index(n)]+csv2L[i]
    elif i!=0  :### n not in file 1 - excluding first row
      if n=='' :
        print("##  Empty key ('') in file 2.")
      else :
        featureD[n]=['']+csv2L[i]
        print("    Entry "+n+" (file 2) was not found.")
    else :### first row
      print("")
      print("    Note: the first rows are merged although they do not have the same label:")
      print(n)
      print(csv1LT[i][0])
      print("")
  print("")
  outname = re.sub('[.]*c*s*v*$','_extended.csv',csv_file1)
  csvNewL=[]
  for csvL in csv1L :
    csvline=csvL+['']*len(csv2L[0]) ### extend each row for the new features
    n=csvline[0]
    if n in featureD.keys() :
      featurenumber=len(featureD[n])
      for j in range(1,featurenumber) :
        csvline[len(csvline)-featurenumber+j]= featureD[n][j]
    csvNewL.append(csvline)
  missing=[k for k in featureD.keys() if k not in csv1LT[0]] ### items that are in file 2 but not in file 1
  print("")
  if len(missing)>0 :
    print("####################################################################################")
    print("##                                                                                ##")
    print("##  Some item(s) in the second file are not present in the first file.            ##")
    print("##                                                                                ##")
    prompt="##  Do you want these items added to the merged table?                            ##\n##                                                                                ##\n##  Type 'y' for yes (return or any other key for no):                            ##\n##  "
    answer= str.strip(get_input(prompt))
    if answer=='y': ### add items from file 2 to the merged file.
      for n in missing:
        csvline=['']*len(csv1L[0])+['']*len(csv2L[0])
        csvline[0]=n
        featurenumber=len(featureD[n])
        for j in range(1,featurenumber) :
          csvline[len(csvline)-featurenumber+j]= featureD[n][j]
        csvNewL.append(csvline)
  print("##                                                                                ##")
  print("####################################################################################")
  print("                                                                                  ##")
  print("    The merged table is written to file:                                          ##")
  print("    "+outname)
  print("                                                                                  ##")
  outfile = open(outname,'w')# start of data upload from file
  for line in csvNewL :
    if len(line) != len(csvL+['']*len(csv2L[0])): ### extend each row for the new features
###      print("##  Line length discrepancy detected. Check outcome carefully!: "+str(len(line))+" "+str(len(csv1L))+" "+csv1L[0])
      print("##  Line length discrepancy detected. Check outcome carefully!: ##")
      print(str(len(line)))
      print(str(len(csv1L)))
      print(str(csv1L[0]))
###      print("##  Line length discrepancy detected. Check outcome carefully!: "+str(len(line))+" "+str(len(csv1L))+" "+csv1L[0])
##    outfile.write('\"'+'\",\"'.join(line)+'\"\n') ### write items with quotes?
    outfile.write(",".join(line)+"\n")

### Here is the start. 
### Under Linux execute by typing in the shell: ./python_data_match.py 
### Under Windows: execute by clicking on "python_data_matcher.py" icon.
### The csv file names will be interactively asked for.
print("####################################################################################")
print("##                                                                                ##")
print("##                          Python Data Matcher                                   ##")
print("##                                                                                ##")
print("##  Joining two 'csv' spread sheets by matching the items in the first columns.   ##")
print("##                                                                                ##")
print("##  The script will ask for the names of two spread sheet files to match.         ##")
print("##  The files have to be in 'csv' (comma seperated value) format.                 ##")
print("##  Please type valid file names otherwise the script will fail.                  ##")
print("##                                                                                ##")
print("##  Python Data Matcher written by Luzie U. Wingen, JIC, 8 November 2018.         ##")
print("##                                                                                ##")
print("####################################################################################")
load_csv_and_add()
print("####################################################################################")
print("##                                                                                ##")
print("##            All work done, Python Data Matcher says Good bye!                   ##")
print("##                                                                                ##")
print("####################################################################################")
print("")
prompt = "Press return to finish."
waiting = str.strip(get_input(prompt))

