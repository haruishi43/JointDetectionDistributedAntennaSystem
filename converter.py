import os
import sys
import scipy.io as sio
import numpy as np
from PIL import Image
from tqdm import *


## Params:
rootdir = os.path.dirname(os.path.realpath(__file__))
current_dir = rootdir

rb = 24
time_interval = 1

height = 3
width = 14

data_root = os.path.join(rootdir, 'old_data')
if not os.path.exists(data_root):
    print("Error: cannot find old data")
    quit()

all_data_paths = os.listdir(data_root)

## Make a data directory if not created yet
new_data_dir = os.path.join(rootdir, 'new_data')
if not os.path.exists(new_data_dir):
    print("Making a new directory: 'data'\n")
    os.makedirs(new_data_dir)


## Functions:
# normalization to change into pixel values
def normalize_nparray(array):
    if type(array).__module__ == np.__name__:
        # normalize
        array *= 255.0/array.max()
        #return rounded to int values
        return np.rint(array)
    else:
        print("error at normalize_nparray")


# saving to image using PIL
def save_as_image(channel, throughput, save_dir, index, count):
    try:
        # concatenate the array
        new_array = np.insert(channel, [0], [throughput[0], throughput[1], throughput[2]], axis=1)
        new_array = np.insert(new_array, [2], [throughput[0], throughput[1], throughput[2]], axis=1)
        new_array = np.insert(new_array, [4], [throughput[0], throughput[1], throughput[2]], axis=1)
        new_array = np.insert(new_array, [6], [throughput[0], throughput[1], throughput[2]], axis=1)
        new_array = np.insert(new_array, [8], [throughput[0], throughput[1], throughput[2]], axis=1)
        new_array = np.insert(new_array, [10], [throughput[0], throughput[1], throughput[2]], axis=1)
        new_array = np.insert(new_array, [12], [throughput[0], throughput[1], throughput[2]], axis=1)
        new_array = new_array.flatten()
        pixels = new_array.tolist()

        # init the image
        image = Image.new("L", (width, height))

        # make the image
        image.putdata(pixels)

        # save it to the path

        complete_name = os.path.join(save_dir, str(index) + "_" + str(count) + ".png")
        image.save(complete_name)

        return 0
    except:
        print("error at save_as_image")
        return 1


## Script:

for i, data_path_top in enumerate(all_data_paths):
    print(data_path_top)
    current_data_path = os.path.join(data_root, data_path_top)
    data_paths = os.listdir(current_data_path)
    print("# of data: "+str(len(data_paths)))
    sys.stdout.flush()
    count = 0
    pbar = tqdm(total=len(data_paths))
    for trial in data_paths:
        if trial.startswith('.'):
            # ignore the files that starts with '.'
            print('will not iterate through this directory')
        else:
            #print(trial)

            # go into each trial
            trial_dir = os.path.join(current_data_path, trial)

            os.chdir(trial_dir)

            # mat files
            mat_files = os.listdir(trial_dir)

            # import data from mat files:
            channel_response = np.array(sio.loadmat(mat_files[0])['channel_response'])
            combination = np.array(sio.loadmat(mat_files[1])['combination'])
            past_throughput = np.array(sio.loadmat(mat_files[2])['past_throughput'])

            # go back to root
            os.chdir(rootdir)

            comb_list = []
            for j in range(0, rb):
                # make the channel_response array for resource block
                channel = normalize_nparray(channel_response[i])
                comb = combination[i,0]

                # make a combination dir
                comb_dir = os.path.join(new_data_dir, str(comb))
                if not os.path.exists(comb_dir):
                    os.makedirs(comb_dir)

                if not comb in comb_list:
                    if len(os.listdir(comb_dir)) > 999:
                        #print("dir full")
                        continue
                    else:

                        # make past throughput array for the time interval
                        throughput = normalize_nparray(past_throughput[:])

                        # save the Image
                        saved = save_as_image(channel, throughput, comb_dir, i, count)

                        if saved != 0:
                            print("error")

                        # cleanup
                        count = count + 1
                        comb_list.append(comb)

            pbar.update(1)
            sys.stdout.flush()
    pbar.close()
