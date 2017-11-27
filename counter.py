import os
import sys
import csv

'''
count the amount of data
'''

rootdir = os.path.dirname(os.path.realpath(__file__))
current_dir = rootdir

new_data_dir = os.path.join(rootdir, 'data')
try:
    current_dir = os.path.join(rootdir, new_data_dir)
except:
    print("no data directory")

# make csv file
with open("count.csv", "w", newline='') as f:

    #print(len(os.listdir(current_dir)))
    for trial in os.listdir(current_dir):


        count = len(os.listdir(os.path.join(current_dir, trial)))

        arr = [trial, count]
        print(arr)

        # save to csv
        writer = csv.writer(f)
        writer.writerow(arr)
        sys.stdout.flush()
