"""
This script is the second of the four scripts
1. Import script
2. Outlier script
3. Filtering script
4. Analysis script
established to analyse the output of the Light/Dark transition test.
It applies a threshold based on fish lost during the tracking in order to
eliminate these outliers from the dataset. The script operates on the output of
 the import (third) script

It was developed at the Computational Ecology working group,
Institute for Environmental Research, Biology V, RWTH Aachen.

For questions please contact: dominik.ziaja@rwth-aachen.de
"""
import pandas as pd
import pathlib
import os
import csv


# set up a multiindex containing Trial time and ID for the dataframe
def set_indices(df):
    df2 = df.copy()
    # set first level of the multiindex
    df2.set_index(['Trial_time [s]'], inplace=True, drop=True)
    # set second level of the multiindex
    df2.set_index(['ID'], append=True, inplace=True)
    # return the dataframe
    return df2


def remove_outliers(threshold, df):
    # roll over 1500 values (60seconds)
    # and calculate the sum for each Individual
    rolled_result = (df.groupby(['Individuum'])
                     ['Distance_moved [mm]'].rolling(window=1500).sum())
    # if a fish moved over the defined threshold,
    # the ID of the fish will be saved
    outliers = rolled_result.index[rolled_result
                                   > threshold].get_level_values('ID').unique()

    # The function returns a list of outliers and
    # the dataframe where outliers are removed
    return outliers, df.drop(outliers, level='ID')


def rearrange_columns(df):
    # reset the dataframe indices so the indices "Trial time" and "ID" become
    # columns again
    df.reset_index(drop=False, inplace=True)
    # get a list of all the dataframe header
    header = df.columns.values.tolist()
    # rearrange the headers
    new_order = header[:1] + header[2:] + header[1:2]
    # update the arrangement of the columns
    df = df[new_order]

    return df


threshold = 750  # > 750 mm per minute moved will be determined as outlier
                 # (1500 window width = 60 seconds)
# set the location of the script as the current working directory
script_location = os.path.dirname(os.path.abspath(__file__))
os.chdir(script_location)
# and create a pathlib path of the working directory
path = pathlib.Path(os.getcwd())
#
pathcontainer = [filepath for filepath in path.glob('**/*.csv')
                 if (('_processed' not in filepath.name)
                     & (
                        ('_R_' in filepath.name)
                        | ('_Replikat_' in filepath.name)
                        ))]

for counter, file in enumerate(pathcontainer):
    print("processing file number {} of {}"
          .format(counter+1, len(pathcontainer)))
    # read in the file
    df = pd.read_csv(file)
    # update all the indices of the dataframe for further analysis
    df2 = set_indices(df)
    # identify and remove the outliers which
    # have more than 750 mm movement within a minute
    outliers, df2 = remove_outliers(threshold, df2)
    # check if any outliers exist
    if outliers.any():
        print("removed the fishs {} due to Movement > {}mm within a minute"
              .format(outliers, threshold))
        # and write them into the file "outliers.txt"
        with open('outliers.txt', 'w', newline='') as outlierfile:
            writer = csv.writer(outlierfile, delimiter=',')
            # format outliers as list to be iterable
            for line in [outliers]:
                writer.writerow(line)

    df2 = rearrange_columns(df2)
    # write the dataframe to csv without the index (Trial time)
    print("writing {} to the harddisk".format(df2.ID[10]))
    df2.to_csv(file.stem+'_wo_outliers.csv',
               sep=',', index=False)
