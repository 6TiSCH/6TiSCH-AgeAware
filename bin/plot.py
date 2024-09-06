"""
Plot a stat over another stat.

Example:
    python plot.py --inputfolder simData/numMotes_50/ -x chargeConsumed --y aveLatency
"""
from __future__ import print_function

# =========================== imports =========================================

# standard
from builtins import range
import os
import argparse
import json
import glob
from collections import OrderedDict
import numpy as np
import time

# third party
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

# ============================ defines ========================================

KPIS = [
    'latency_max_s',
    'latency_avg_s',
    'latencies',
    'lifetime_AA_years',
    'sync_time_s',
    'join_time_s',
    'upstream_num_lost',
    'aoi'
]

# ============================ main ===========================================

def main(options):

    # hard-coded config when debugging
    #options.inputfolder = '/home/sekiro/6tisch-arshi/6tisch-simulator-a015d94da186/bin/simData'
    #options.inputfolder = 'C:/Arshia/6TiCSH-AgeAware/bin/simData'

    bin_folder = os.getcwd()
    if os.path.exists(os.path.join(bin_folder, 'bin')) == False:
        bin_folder = os.path.join(bin_folder, '..')
        if os.path.exists(os.path.join(bin_folder, 'bin')) == False:
            bin_folder = os.path.join(bin_folder, '..')
            if os.path.exists(os.path.join(bin_folder, 'bin')) == False:
                bin_folder = os.path.join(bin_folder, '..')
                if os.path.exists(os.path.join(bin_folder, 'bin')) == False:
                    print('ERROR: could not find bin folder:', bin_folder)
                    return
    
    bin_folder = os.path.join(bin_folder, 'bin')
    options.inputfolder = os.path.join(bin_folder, 'simData')

    # init
    data = OrderedDict()

    # chose lastest results
    subfolders = list(
        [os.path.join(options.inputfolder, x) for x in os.listdir(options.inputfolder)]
    )
    subfolder = max(subfolders, key=os.path.getmtime)

    for key in options.kpis:
        # load data
        for file_path in sorted(glob.glob(os.path.join(subfolder, '*.kpi'))):
            curr_combination = os.path.basename(file_path)[:-8] # remove .dat.kpi
            with open(file_path, 'r') as f:

                # read kpi file
                kpis = json.load(f)

                # init data list
                data[curr_combination] = []

                # fill data list
                for run in kpis.values():
                    for mote in run.values():
                        if key in mote:
                            data[curr_combination].append(mote[key])

        # plot
        try:
            if key in ['lifetime_AA_years', 'latencies']:
                plot_cdf(data, key, subfolder)
            elif key == 'aoi':
                plot_aoi(data,subfolder)
            else:
                plot_box(data, key, subfolder)

        except TypeError as e:
            print("Cannot create a plot for {0}: {1}.".format(key, e))
    print("Plots are saved in the {0} folder.".format(subfolder))

# =========================== helpers =========================================

def plot_cdf(data, key, subfolder):
    for k, values in data.items():
        # convert list of list to list
        if type(values[0]) == list:
            values = sum(values, [])

        values = [None if value == 'N/A' else value for value in values]
        # compute CDF
        sorted_data = np.sort(values)
        yvals = np.arange(len(sorted_data)) / float(len(sorted_data) - 1)
        plt.plot(sorted_data, yvals, label=k)

    plt.xlabel(key)
    plt.ylabel("CDF")
    plt.legend()
    savefig(subfolder, key + ".cdf")
    plt.clf()

def plot_box(data, key, subfolder):
    plt.boxplot(list(data.values()))
    plt.xticks(list(range(1, len(data) + 1)), list(data.keys()))
    plt.ylabel(key)
    savefig(subfolder, key)
    plt.clf()

def plot_aoi(data,subFolder):
    
    for k, values in data.items():
        for i in range (0 ,len(values)):
            asn_values = [int(item["asn"]) for item in values[i]]
            aoi_values = [int(item["aoi"]) for item in values[i]]

            # Plotting the data
            plt.figure(figsize=(10, 6))
            plt.plot(asn_values, aoi_values, marker='o', linestyle='-', color='blue')

            # Adding labels and title
            plt.xlabel('Time')
            plt.ylabel('AOI')
            plt.title('Time vs AOI')
            plt.grid(True)
            savefig(subFolder, "aoi_"+str(i)+ ".cdf")
            plt.close()

            plot_moving_average(asn_values, aoi_values,subFolder=subFolder,index=i)
            
            plot_peaks(asn_values, aoi_values,subFolder=subFolder,index=i)

            time.sleep(0.2)  

def plot_moving_average(asn_values, aoi_values, window_size=5, subFolder="", index=0):
    if len(aoi_values) == 0:
        print("Warning: AOI values are empty, skipping moving average plot.")
        return
    
    moving_averages = np.convolve(aoi_values, np.ones(window_size)/window_size, mode='valid')
    
    if len(moving_averages) == 0:
        print("Warning: Moving averages result is empty, skipping plot.")
        return
    
    plt.figure(figsize=(10, 6))
    plt.plot(asn_values[:len(moving_averages)], moving_averages, marker='o', linestyle='-', color='green')
    plt.xlabel('Time')
    plt.ylabel('Moving Average of AOI')
    plt.title('Moving Average (window size={}) of AOI'.format(window_size))
    plt.grid(True)
    
    savefig(subFolder, "moving_average_aoi"+str(index)+ ".cdf")

    plt.close()

def plot_peaks(asn_values, aoi_values, window_size=5, subFolder="", index=0):
    if len(aoi_values) == 0:
        print("Warning: AOI values are empty, skipping peak plot.")
        return
    
    moving_averages = np.convolve(aoi_values, np.ones(window_size)/window_size, mode='valid')
    
    if len(moving_averages) == 0:
        print("Warning: Moving averages result is empty, skipping plot.")
        return
    
    peaks, _ = find_peaks(moving_averages)
    
    if len(peaks) == 0:
        print("Warning: No peaks found, skipping peak plot.")
        return
    
    plt.figure(figsize=(10, 6))
    plt.plot(asn_values[peaks], moving_averages[peaks], marker='o', linestyle='-', color='red', label='Peaks')
    
    plt.xlabel('Time')
    plt.ylabel('Peaks of Moving Average of AOI')
    plt.title('Peaks of Moving Average (window size={}) of AOI'.format(window_size))
    plt.grid(True)
    
    savefig(subFolder, "peaks_moving_average_aoi_" + str(index) + ".cdf")

    plt.close()


def savefig(output_folder, output_name, output_format="png"):
    # check if output folder exists and create it if not
    if not os.path.isdir(output_folder):
        os.makedirs(output_folder)

    # save the figure
    plt.savefig(
        os.path.join(output_folder, output_name + "." + output_format),
        bbox_inches     = 'tight',
        pad_inches      = 0,
        format          = output_format,
    )

def parse_args():
    # parse options
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--inputfolder',
        help       = 'The simulation result folder.',
        default    = 'simData',
    )
    parser.add_argument(
        '-k','--kpis',
        help       = 'The kpis to plot',
        type       = list,
        default    = KPIS
    )
    parser.add_argument(
        '--xlabel',
        help       = 'The x-axis label',
        type       = str,
        default    = None,
    )
    parser.add_argument(
        '--ylabel',
        help       = 'The y-axis label',
        type       = str,
        default    = None,
    )
    parser.add_argument(
        '--show',
        help       = 'Show the plots.',
        action     = 'store_true',
        default    = None,
    )
    return parser.parse_args()

if __name__ == '__main__':

    options = parse_args()

    main(options)
