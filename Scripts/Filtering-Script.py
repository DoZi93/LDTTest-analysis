"""
This script is the fourth of the four scripts
1. Import script
2. Outlier script
3. Filtering script
4. Analysis script
established to analyse the output of the Light/Dark transition test.
It applies the moving average and moving standard deviation to extract
characteristics from the dataset and eliminate the zero-values. The script
 operates on the output of the outlier (third) script.

It was developed at the Computational Ecology working group,
Institute for Environmental Research, Biology V, RWTH Aachen.

For questions please contact: dominik.ziaja@rwth-aachen.de
"""

# Script for combining several datasets in the folder and calculating the
# rolling mean/StdDev. Only Negative Control and max. Concentration are kept
# Except for Ethanol, where the max Concentration
# resulted in death, thus 2nd highest Concentration is kept
import pandas as pd
import numpy as np
import pathlib
import os


# Get a list of all csv-files with "processed_with_na" in their name
# as this is the name saved by the script beforehand
def get_file_paths(path):
    filepaths = [filepath for filepath in
                 pathlib.Path(path).glob('*.csv')]

    return filepaths


# format the dataframe into a trialtime x sample size format
# where the header consists of each fish's ID
# and each row of one timepoint
# while the values in the actual df are the Distance moved
def format_dataframe(df):
    df2 = pd.DataFrame()
    for Individuum in df['Individuum'].unique().tolist():
        subset_Ind = df[df['Individuum'] == Individuum]
        ID = subset_Ind['ID'].unique()[0]
        # needs to be defined as f array beforehand, otherwise dataframe wont
        # accept the values
        array = np.asfarray(subset_Ind['Distance_moved [mm]'])
        # create new column with the individual ID as name
        df2[ID] = array

    return df2


# only works if all fish have the equal amount of time points
# Thus, na-values in the dataset need to be dropped after this, not beforehand
def set_index_trialtime(df_all, df_singlefish):
    df_all.set_index(df_singlefish['Trial_time [s]'].unique(), inplace=True)


def add_Dataframes_together(df, df_together):
    for key, values in df.iteritems():
        # allows duplicates, however - in best case no IDs are duplicates.
        # only happens if replicate ID in the metafile is used for another
        # experiment as well
        df_together.insert(len(df_together.columns), key, values,
                           allow_duplicates=True)

    return df_together


# set the path where the script is located as the current working directory
script_location = os.path.dirname(os.path.abspath(__file__))
os.chdir(script_location)
# and print it out for control
print("operating in datapath {} ".format(os.getcwd()))
# get all csv files with "na_processed" in the scripts directory.
datafile_paths = get_file_paths(os.getcwd())

# create an empty dataframe to collect all values in it
df_all = pd.DataFrame()


for idx, filepath in enumerate(datafile_paths, 1):
    df = pd.read_csv(filepath, sep=',', header=0)

    print(f"file number {idx} of {len(datafile_paths)} is being processed.")
    # extract only the min and max concentration (NegControl, max Conc)
    # of each treatment except:
    # if EtOH is the treatment,
    # Then NegaControl and second highest Concentration are selected
    # Yet "EtOH" needs to be written that way.
    # Maybe make it case insensitive check
    if np.isin(df['Substance'].unique(), 'EtOH').any():
        df_subset = df[
                np.logical_or(
                             df['Concentration']
                             == df['Concentration'].min(skipna=True),
                             df['Concentration']
                             == np.sort(df['Concentration'].unique())[-2])
                ]
    else:
        df_subset = df[
                (df['Concentration'] == df['Concentration'].max())
                | (df['Concentration'] == df['Concentration'].min())]
    # set each ind. to a column in a transposed dataframe
    df2 = format_dataframe(df_subset)

    df_all = add_Dataframes_together(df2, df_all)

# set the trial time as index of the big dataframe
try:
    df_all.set_index(df['Trial_time [s]'].unique(), inplace=True)
# if an exception is raised, print the following error.
except:
    print("Exception arised, maybe the dataframes don't"
          "have the equal length of Time points")

# delete the big dataframe to get some RAM back
del df
# when every file is concatenated, drop na-values
df_all.dropna(inplace=True)
# define window width of the filters
rolling = df_all.rolling(center=True, window=12000)
# and apply the filters, save it transposed (fish x timepoints format)
rolling_mean = rolling.mean().T
rolling_stddev = rolling.std().T
# save the dataframe unrolled as well as rolled
print("Writing the combined dataframe without "
      "applied filters to the harddisk")
df_all.to_csv(f'Fish_behaviour_unfiltered.csv',
              index=True, header=True, sep=',')
print("Writing the moving standard deviation to the harddisk")
rolling_stddev.to_csv(f'Fish_behaviour_moving_stddev.csv',
                      index=True, header=True, sep=',')
print("Writing the moving average to the harddisk")
rolling_mean.to_csv(fr'Fish_behaviour_moving_average.csv',
                    index=True, header=True, sep=',')
