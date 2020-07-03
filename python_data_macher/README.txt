Author: Luzie U. Wingen, John Innes Centre
Date: 14/11/2019

Outline: script to match items in the first column of two different spread sheets.
Useful to match up data sets from different sources, e.g. Watseq phenotype data from different parties put together
in one big csv spread sheet.

Software: Python 3.5 or higher

Basic instructions:

The script "python_data_matcher.py" works only on "csv" files. Save spread sheets from Excel as "csv".
Install Python 3.5 or higher.

On Windows click on "python_data_matcher.py" icon.

A window will pop up and ask you for the name of your first file.
You best put the file in the same directory as the python script - less typing.

You have to type the name  and hit "return".

The script will then ask you for the name of the second file. 
You have to type the name of the second file and hit "return".

The script then reports back if it finds items in the second file not present in the first. 
If you want items added to your list, you answer that question with:
 "y" + "return".
If you don't want items added just hit "return".

"python_data_matcher.py" will then produce the joint file, using the first file name and adds 'extended' to it.

All done, you can hit "return" to leave the script.

