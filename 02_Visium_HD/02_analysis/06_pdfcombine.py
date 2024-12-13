# Script to combine GEP pdfs in a folder
import os
import sys
from PyPDF2 import PdfFileMerger
from natsort import natsorted

input1 = sys.argv[1]
input2 = sys.argv[2]

print(input1)
print(input2)

rank = "K" + str(input1)

directory = str(input2)
directory = os.path.join(directory, "hires_vectors")

outputname = directory + "/Combined_plots_" + rank + ".pdf"

pdfs = []

merge = PdfFileMerger()
rank = rank + "_"

# Iterates through the directory and merges the pdfs of a certain rank
for filename in os.listdir(directory):

    # String for each filename
    f = os.path.join(directory, filename)

    # Checks if it is a file
    if os.path.isfile(f):

        # Checks if the rank# is in the file name
        if rank in f:

            # Adds the file to the list to merge
            pdfs.append(f)

# Sorts the list
pdfs = natsorted(pdfs)
# Merges and writes the output file
for pdf in pdfs:
    merge.append(pdf)

merge.write(outputname)

