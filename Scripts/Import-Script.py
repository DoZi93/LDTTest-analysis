"""
This script is the first of the four scripts
1. Import script
2. Outlier script
3. Filtering script
4. Analysis script
established to analyse the output of the Light/Dark transition test.
It uses pandas to combine raw datafiles and the metainformation obtained
to form a single dataframe for each experiment. The script operates on the
 DanioVision chambers output and the metafiles.

It was developed at the Computational Ecology working group,
Institute for Environmental Research, Biology V, RWTH Aachen.

For questions please contact: dominik.ziaja@rwth-aachen.de
"""

import numpy as np
import csv
import pandas as pd
import re
import pathlib
import os


def compress_list(nested_list):
    unnested_list = []
    for sublist in nested_list:
        for item in sublist:
            unnested_list.append(item)
    return unnested_list


def list_files(pathlib_path):
    r = []
# find all files ending with .txt in the dirrefectory
    for file in pathlib_path.glob('**/*.txt'):
        r.append(file.as_posix())
    return r


def process_wellplate_metafile(metafile):
    # create two empty lists for both informations
    # treatments and the fish_positions
    treatments = []
    wellplate_position = []
    for row in metafile[6:]:
        # check whether the row is empty or contains something
        if row:
            # check whether '#' is at the beginning like in
            # "#2_17.5 mg/L Cadmiumchlorid"
            if re.match('[#]', row[0]):
                # append only the number and treatment (17.5 mg/L Substance)
                treatments.append(re.split('[_#]', row[0])[1:])
        # if the first element is either 1 to 2 numbers long OR contains nan
        # ("1" "11" "nan") (number up to 99 [99 treatments] is possible)
        # it will be recognized as a number indicating treatment and
        # wellplate_position
            elif re.match('[\b\d{1,2}\b|^nan]', row[0]):
                wellplate_position.append(row)

    return wellplate_position, treatments


# this function replaces the numbers in the read in meta_wellplate file
# with the actual concentration + Substances
def replace_wellplate_treatments(wellplate_position, treatments):
    # create a vector of the positions of the fish in 'U32' dataformat as in an
    # default format 'U1' it will only store 1-character long strings
    wellplate_array_flat = np.array(wellplate_position, dtype='U32').flatten()
    # iterate through the different treatments and replace the number
    # with the corresponding treatment e.g.'2' = '17.5 mg Cadmiumchlorid'
    for treatment in treatments:
        # replaces every element with the corresponding treatment.
        np.place(wellplate_array_flat,  # array where to replace
                 wellplate_array_flat == treatment[0],  # condition
                 treatment[1])  # value to insert

    return wellplate_array_flat


def read_in_file(datapath):
    x = []
    with open(datapath) as csvfile:
        data = csv.reader(csvfile, delimiter=';')
        for row in data:
            # 'if row' checks whether the row is empty or not and only appends
            # if the row contains something
            if row:
                x.append(row)

    return x


# format the data into pandas dataframe and return both, headers and values
def format_data_into_dataframe(data):
    header_raw = data[:35]
    values = pd.DataFrame(data[35:])
    # spacebars of columnnames (line 34) are replaced with '_'
    header_replaced = [i.replace(' ', '_') for i in header_raw[33]]
    # create a dictionary of old and new header
    new_columnnames = {header_old: header_new
                       for header_old, header_new in
                       zip(range(0, len(header_replaced)), header_replaced)}
    # and replace old headers in the dataframe with the new ones
    dataframe = values.rename(columns=new_columnnames)

    return dataframe, header_raw


# format the values into float type to allow calculations and boolean checks
def format_values(df):
    df['Trial_time'] = df['Trial_time'].astype(float)
    df['Distance_moved'].replace(to_replace='-', value=np.nan, inplace=True)
    df['Distance_moved'] = df['Distance_moved'].astype(float)

    return df


# Updates the header and inserts units ('mm' and 's')
def update_column_labels(df):
    old_labels = df.columns.values.tolist()
    new_labels = ['Trial_time [s]', 'Distance_moved [mm]']
    dict_labels = {old_name: new_name
                   for old_name, new_name in zip(old_labels, new_labels)}
    df = df.rename(columns=dict_labels)

    return df


# uses the light_dark information to insert 0 for light off, and 1 for light on
def insert_lighton_lightoff(dataframe, light_dark_metadata):
    # make a new empty column in the dataframe to insert the light/dark-info
    dataframe['Light_on_off'] = 0
    # iterate through the different time borders.
    for idx, (border, lightcondition) in enumerate(light_dark_metadata):

        # if last transitionpoint is reached
        # do 'smaller/equal' instead of 'between' comparison
        # size -1 because max index != size
        if idx == len(light_dark_metadata)-1:
            # set values in the Light column to light on/light off
            # based on the Trial time
            dataframe.loc[dataframe['Trial_time [s]'] >= border,
                          'Light_on_off'] = lightcondition
        # if the border is between different borders
        else:
            # ask which values are in between the two timeborders
            condition = ((dataframe['Trial_time [s]'] >= border)
                         & (dataframe['Trial_time [s]']
                         < light_dark_metadata[idx+1][0]))
            # set the values between timeborders to the lightcondition
            dataframe.loc[condition, 'Light_on_off'] = lightcondition

    # get the rows in the dataframe where the light is on/light is off
    # and replace light on with 1
    dataframe.loc[dataframe['Light_on_off'] == 'Light on', 'Light_on_off'] = 1
    # and replace light off with 0
    dataframe.loc[dataframe['Light_on_off'] == 'Light off', 'Light_on_off'] = 0

    return dataframe


# set the path where the script is located as the current working directory
script_location = os.path.dirname(os.path.abspath(__file__))
os.chdir(script_location)
# and print it out for control
print("operating in datapath {} ".format(os.getcwd()))
# get the paths of all files located in the current working directory in a list
operating_path = pathlib.Path(os.getcwd())
files_path_list = [file_path.as_posix()
                   for file_path in operating_path.glob('*')]
# generate an empty dataframe where all fish_files will be stored in
# and an empty list to append the paths of the raw files for each fish
data = pd.DataFrame()
fishmovement_file_paths = []

# iterate through all files in the folder where the script is located
for item in files_path_list:
    pathlib_item = pathlib.Path(item)
    # check if 'meta' is in a case-insensitive version of the filename
    if ('meta' in pathlib_item.name.lower()):
        # then read the file as "metafile" in
        with open(item) as csvfile:
            metafile_reader = csv.reader(csvfile, delimiter=',', quotechar='"')
            metafile = [line for line in metafile_reader]
    # check which of the three possible metafiles it is
    # and process them accordingly.
    # check case insensitive whether the metafile
    # contains the light/dark information
        if (('light' in pathlib_item.name.lower())
                or ('dark' in pathlib_item.name.lower())):
            light_dark_meta = [
                               [float(border), light]
                               for border, light in metafile[6:]]
    # or the treatment and wellplate position information
        elif ('wellplate' in pathlib_item.name.lower()):
            wellplate, treatment = process_wellplate_metafile(metafile)
            treatment = replace_wellplate_treatments(wellplate, treatment)
            #treatment = (
            #            replace_wellplate_treatments(
            #                            process_wellplate_metafile(metafile)))
    # or the general metainformations (hpf, etc.)
        elif ('expdesign' in pathlib_item.name.lower()):
            exp_design_meta = compress_list(metafile[6:])
    # if 'meta' is not in the filename, it is recognized as a raw data file and
    # appended in the fishmovement_list
    elif (('.txt' in item) and ('hardware' not in item.lower())):
        fishmovement_file_paths.append(item)

# logical check whether less than 96 files were detected, indicating
# some files might have been forgotten to be inserted in the folder
if len(fishmovement_file_paths) < 96:
    print("Are you sure, all fish files are in the folder and the names are "
          "correctly formatted? I register < 96 fishfiles")
    input("Press Enter to continue...")
# store informations about hpf and the replicate_ID in a variable
hpf = exp_design_meta[1]
replicate = exp_design_meta[3]

for idx, filepath in enumerate(fishmovement_file_paths):
    # print the filenumber which is being processed on the display
    print("processing file {} of {}.".format(
        idx+1, len(fishmovement_file_paths)))
    # read in the csv file
    temp = read_in_file(filepath)
    df, header = format_data_into_dataframe(temp)
    # subset the dataframe to df2 keeping only the columns which are necessary
    df2 = df[['Trial_time', 'Distance_moved']].copy()
    # format trial time+Distance moved into floats
    df2 = format_values(df2)
    # insert units into the column labels
    df2 = update_column_labels(df2)
    # insert the meta information about light/dark times
    df2 = insert_lighton_lightoff(df2, light_dark_meta)

    # 1. get the Individuum number the raw data file describes
    # +1 to make the range from 1 to 96 instead of 0 to 95
    individuum_number = int(header[6][1])+1
    # and set it as a new column
    df2['Individuum'] = individuum_number

    # 2. add new columns containing info about the treatment to the dataframe
    Concentration_Substance = treatment[individuum_number-1].split(' ')
    # replace ',' in the concentration with '.'
    df2['Concentration'] = Concentration_Substance[0].replace(',', '.')
    df2['Concentration_unit'] = Concentration_Substance[1]
    df2['Substance'] = Concentration_Substance[2]
    df2['hpf'] = hpf

    # set up an unique ID of the fish
    df2 = df2.assign(
        # the new column 'ID' is defined as the Individuum_numbers
        ID=(str(df2.iloc[1].Individuum)
            + '_'
            # + the first six letters of the Substance
            + ''.join(re.split('', df2.iloc[1].Substance)[:6])
            # + the concentration
            + str(df2.iloc[1].Concentration)
            + '_'
            # + the hpf
            + str(df2.iloc[1].hpf)
            + 'hpf'
            + '_'
            # + the ID of the replicate
            + str(replicate))
            )
    # then the processed information is appended to the big data-dataframe
    data = data.append(df2)
    # and the next raw file will be processed and appended the same way
print("Writing the dataframe onto the harddisk...")
# save the pandas dataframe containing all processed informations of the folder
# as a csv in the folder with the replicate name
data.to_csv(f"Behaviour_df_{replicate}.csv",
            sep=',', index=False, na_rep='nan')
