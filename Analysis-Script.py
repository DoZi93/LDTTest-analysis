"""
This script is the fourth of the four scripts
1. Import script
2. Outlier script
3. Filtering script
4. Analysis script
established to analyse the output of the Light/Dark transition test.
It applies hierarchical clustering and visualizes the results to allow for
interpretation. The script operates on the output of the filtering (third)
 script.

It was developed at the Computational Ecology working group,
Institute for Environmental Research, Biology V, RWTH Aachen.

For questions please contact: dominik.ziaja@rwth-aachen.de
"""

import matplotlib as mlb
import pathlib
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns; sns.set()
from scipy.cluster import hierarchy
import re
import matplotlib.patches as mpatches
import os


def get_file_paths(folderpath):
    # append all csv-files in the directory containing (un)rolled in their name
    path_list = [item for item in
                 pathlib.Path(folderpath).glob('*.csv')
                 if re.search('_moving_', item.name)]

    return path_list


# join the Susbtance, Concentration and hpf from the ID to remove the
# "uniqueness" of each ID while keeping the treatment
# for color-visualization later
def get_treatments(df):
    IDs = pd.Series(df.index).apply(
            lambda row: '-'.join([re.split('_', row)[1],
                                 re.split('_', row)[2]])
            )
    return np.array(IDs)


# get all negative controls once, appended in a list
def get_Negctrls_unique(IDs):
    Negative_controls = []

    for nc in IDs:
        # '0-' indicates Negative control
        if ('0-' in nc) and (nc not in Negative_controls):
            Negative_controls.append(nc)

    return np.array(Negative_controls)


# replace the unique negative controls with the label "NegCtrl"
def replace_negctrls(IDs, negative_controls):
    index = np.isin(IDs, negative_controls)
    IDs[index] = 'NegCtrl'

    return IDs


# apply the three functions defined to replace the treatment-IDs in the df
def get_treatments_and_replace(df):
    IDs = get_treatments(df)
    neg_ctrls = get_Negctrls_unique(IDs)
    IDs = replace_negctrls(IDs, neg_ctrls)

    return IDs


# ask the user for an input of the amount of clusters he would like to receive
def get_amount_cluster():
    n = input("Please type in how many Clusters you would"
              " like to have returned (positive integer!): ")
    try:
        # try to set the input as integer, raise exception otherwise
        n = int(n)
        if n < 0:
            print(f"{n} is not positive.")
            get_amount_cluster()
        # if everything was entered correct, return n
        return n
    # if no integer number was entered, call the function again
    except ValueError:
        print("You entered a non integer number.")
        get_amount_cluster()


def calculate_hierarchy_linkage(df, amount_cluster):
    df_cluster = df.copy()
    counter = 0
    # list for checking the condition whether there are still 1-fish-cluster
    bool_check_list = [True]
    # will be outputted as the fish removed
    outlier_list = []
    # while there is one cluster with only 1 fish
    while True in bool_check_list:
        # calculate the linkage method metrics and values,
        # cut them at the wished amount of cluster
        # and save them into a pd.Series-format
        link = hierarchy.linkage(df_cluster,
                                 metric='euclidean',
                                 method='complete')
        cut_tree = hierarchy.cut_tree(link, amount_cluster)
        cut_tree = np.squeeze(cut_tree)
        Cluster_series = pd.Series(cut_tree, index=df_cluster.index)
        # every step, generate the check list new
        bool_check_list = []
        # iterate through every cluster
        for cluster in np.unique(cut_tree):
            # check if only 1 individuum is in a cluster
            if np.sum(cut_tree == cluster) <= 1:
                # append the boolean check result if its true
                bool_check_list.append(True)
                # get the ID of the fish which makes up one cluster by himself
                outlier_fish = Cluster_series.index[Cluster_series == cluster]
                # and then drop that fish
                df_cluster.drop(outlier_fish, inplace=True)
                print(f"Outlier fish {outlier_fish} was dropped"
                      " in run number: {counter}")
                # save the fish in a list
                outlier_list.append(outlier_fish[0])

        counter += 1
        # return the Cluster_assignments, list of outliers and the linkage
    return Cluster_series, outlier_list, link


# for the seaborn clustermap, colors which indicate the treatment
# next to the heatmap
def create_row_colors(colors_to_zip, df):
    # Make a dictionary out of the ID and color
    zipped_IDs_and_colors = zip(np.unique(IDs), colors_to_zip)
    color_dictionary = dict(zipped_IDs_and_colors)
    # then for seaborn.clustermap make a mapped pandas series out of it
    row_colors = pd.Series(IDs, index=df.index).map(color_dictionary)

    return row_colors, color_dictionary


def plot_clustermap(df, linkage, row_colors, color_dictionary, fig_title):
    # plot the clustermap with the linkage precalculated
    # don't cluster the columns
    fig = sns.clustermap(df, row_linkage=linkage, col_cluster=False,
                         cmap='Spectral_r',
                         row_colors=row_colors, figsize=(20, 20),
                         yticklabels=False,
                         xticklabels=False, cbar_kws={"shrink": 1.5}
                         )
    # create the treatment legend
    legend_TN = [mpatches.Patch(color=color_dictionary[label], label=label)
                 for label in color_dictionary]
    # then defining where the legend box is anchored (bbox_to_anchor)
    # and specify the location with loc=center
    l2 = fig.ax_heatmap.legend(
                               loc='center', bbox_to_anchor=(0.5, 1.1),
                               handles=legend_TN, frameon=True, ncol=5,
                               prop={'size': 15})
    # define the title of the legend
    l2.set_title(title='Treatment', prop={'size': 30})
    # define the title of the whole plot
    fig.fig.suptitle(f'Clustermap {fig_title}', x=0.6, fontsize=60)
    # then save the figure
    fig.savefig(fr'clustermap_{fig_title}.png', dpi=300)


def plot_clusters(df, Cluster_series, fig_title):
    # iterate over every cluster
    for cluster in Cluster_series['Cluster'].unique():
        print(f"Starting with cluster {cluster}.")
        # select the fish time series assigned to the cluster
        index_df = Cluster_series.index[Cluster_series['Cluster'] == cluster]
        subset_cluster = df.loc[index_df]
        # define size of the figure
        fig = plt.figure(figsize=(10, 5))
        ax = fig.add_subplot(1, 1, 1)
        # iterate over every fishs timeseries
        for idx, (ID, fish) in enumerate(subset_cluster.iterrows()):
            fish.reset_index(inplace=True, drop=True)
            if idx == 0:
                ax.plot(pd.Series(fish), color='blue', label='Time series')
            # and plot it
            else:
                ax.plot(pd.Series(fish), color='blue')
        # plot the mean of the cluster
        mean = subset_cluster.mean(axis=0)
        mean.reset_index(inplace=True, drop=True)
        ax.plot(mean, color='red', label='Mean')
        # set axis properties:
        # xaxis and yaxis limits
        ax.set(xlim=(0), ylim=(0, df.max().max() + (1/10 * df.mean().mean())),
               # figure title with filter applied, amount of fish and cluster
               title=(f"{fig_title.capitalize()}, "
                      f"n={len(subset_cluster.index)}, "
                      f"Cluster: {cluster}"),
               ylabel=f"{fig_title} [mm]",
               xlabel="Trial time")
        ax.set_xticklabels([])
        ax.legend(loc='upper right')
        fig.savefig(f"Plot_{fig_title}_Cluster{cluster}.png",
                    dpi=300)


def autolabel_bar(rectangle, ax, cluster_series, stacked_height):
    for idx, rect in enumerate(rectangle):
        # write the sample size on top of each bar
        ax.text(rect.get_x() + rect.get_width()/2, stacked_height[idx],
                f"n={(cluster_series == idx).sum()[0]}",
                ha='center', va='bottom', fontsize=14)


def plot_stacked_barplot(Cluster_series, crosstab, fig_title):
    # initialize a series that keeps track
    # of the barplot height in each iteration to stack it
    stacked_height = pd.Series(np.zeros(len(crosstab.index)))
    # initialize figure
    fig, ax = plt.subplots(figsize=(7, 7))
    # join the Cluster numbers [0,1,2,...] with "Cluster" to "Cluster 0"
    # for the x axis labeling
    xticks = [' '.join(['Cluster', str(cluster)])
              for cluster in crosstab.index]

    ax.set_ylabel('Absolute amount []', fontsize=14)
    ax.set_xticks(np.arange(len(crosstab.index)))
    ax.set_xticklabels(xticks, rotation=45, fontsize=14)
    # define the bar positions (number of clusters)
    positions_x = np.arange(0, len(crosstab.index))
    # get the colors we want to use
    colormap = mlb.cm.get_cmap('Set3')
    barwidth = 0.85

    # iterate over every column (treatment count per cluster)
    for idx, series in enumerate(crosstab.iteritems()):
        # if the idx=0, no "bottom" keyword needs to be defined
        if idx == 0:
            rect = ax.bar(positions_x, series[1],
                          color=colormap(idx),
                          width=barwidth,
                          label=series[0])
        # all other time, "bottom" defines
        # where the next bar should be plotted on top
        # Herefore the height of the iteration before is used as height
        else:
            rect = ax.bar(positions_x, series[1],
                          color=colormap(idx),
                          bottom=stacked_height,
                          width=barwidth,
                          label=series[0])
            # if its the last iteration,
            # the clusters sample size is plotted on top
        # define the position of the legend in the plot
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        # add up the values, so the stacked barplots are always
        # plotted on the height of the iteration before
        stacked_height += series[1]
        if idx == (crosstab.shape[1]-1):
            autolabel_bar(rect, ax, Cluster_series, stacked_height)
    fig.suptitle(f'{fig_title}', fontsize=18, y=0.95)

    return fig


def create_barplot_fig_title(fig_title):

    fig_title_split = ' '.join(fig_title.split('_')).lower()
    barplot_fig_title = f"Clustercomposition {fig_title_split}"

    return barplot_fig_title


# Set the scripts location as working directory
script_location = os.path.dirname(os.path.abspath(__file__))
os.chdir(script_location)
path_list = get_file_paths(os.getcwd())
# depending on whether its the rolling mean/stddev the figure-title is defined
for path in path_list:
    if 'mean' in path.name.lower():
        fig_title = 'moving average'
    elif 'stddev' in path.name.lower():
        fig_title = 'moving standard deviation'
# read in the file
print(f"loading in file: {path.name}.")
df = pd.read_csv(path, sep=',', header=0, index_col=0)
# drop all na-values
df.dropna(inplace=True, axis=1)
# calculate the linkage and get the
# list of outliers, the linkage and the clusters
print(f"Calculating the linkage.")
Cluster_series, outlier_list, recursive_linkage = (
    calculate_hierarchy_linkage(df, get_amount_cluster())
    )
# drop the outliers from the original dataframe before continuing visualization
df.drop(outlier_list, axis=0, inplace=True)
# get treatment info from the IDs (e.g. 97_EtOH3_96hpf_2 becomes EtOH3-96hpf)
IDs = get_treatments_and_replace(df)
# colors which should be used for the row_colors
# need to match the amount of unique treatments
colors_to_zip = ['orange', 'yellow', 'black', 'springgreen',
                 'darkgreen', 'olive', 'deepskyblue', 'blue',
                 'rosybrown', 'red', 'darkviolet']

row_colors, color_dictionary = create_row_colors(colors_to_zip, df)
print(f"Plotting and saving the clustermap.")
plot_clustermap(df, recursive_linkage, row_colors, color_dictionary, fig_title)
# format the Cluster-assignment dataseries
Cluster_series.sort_values(inplace=True)
Cluster_series = Cluster_series.to_frame('Cluster')
print("Saving the Clustering results as csv file.")
Cluster_series.to_csv(f"{fig_title}_HClustering_results.csv",
                      sep=',',
                      header=True)
print("Plotting and saving every cluster.")
plot_clusters(df, Cluster_series, fig_title)
# Get the IDs with "neg control" again for the sorted Cluster_series dataframe
IDs = get_treatments_and_replace(Cluster_series)
# Get a copy with the updated ID
Cluster_series_new_ID = Cluster_series.set_index(IDs).copy()
# calculate the amount of each treatment present in each cluster
ctb = pd.crosstab(Cluster_series_new_ID['Cluster'],
                  Cluster_series_new_ID.index)
# plot a stacked barplot of the clusters,
# showing the composition of each
print(f"Plotting now the stacked barplot for the composition of each cluster")
fig_title = create_barplot_fig_title(fig_title)
barplot = plot_stacked_barplot(Cluster_series, ctb, fig_title)
barplot.savefig(f'{fig_title}_stacked_barplot.png',
                bbox_inches='tight',
                dpi=300)
