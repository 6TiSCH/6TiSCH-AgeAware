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
    'aoi',
    'Asf_Averages',
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

    # print input folder read from data
    print("Input folder: {0}".format(subfolder))

    # read the config.json file to get settings values
    with open(os.path.join(subfolder, 'config.json'), 'r') as f:
        settings = json.load(f)

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

                if key == 'Asf_Averages':
                    for run in kpis.values():
                        data[curr_combination].append(run['global-stats']['aoi_feedbacks'][2]['value'])

        # plot
        try:
            if key in ['lifetime_AA_years', 'latencies']:
                plot_cdf(data, key, subfolder)
            elif key == 'aoi':
                plot_aoi(data, settings, subfolder)
            elif key == 'Asf_Averages':
                plot_aoi_feedback(data, settings, subfolder)
            else:
                plot_box(data, key, subfolder)

        except TypeError as e:
            print("Cannot create a plot for {0}: {1}.".format(key, e))
    
    
    #plot aoi in 100 run
    list_of_aoi = []
    for key, item in kpis.items():
        aoi_stats = item['global-stats']['aoi_stats'][0]['value']
        if aoi_stats:
            list_of_aoi.append((key, aoi_stats))
    list_of_aoi = sorted(list_of_aoi, key=lambda x: int(x[0]))
    print(list_of_aoi)
    plot_aoi_per_run(list_of_aoi, subFolder=subfolder)
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

def plot_aoi(data, settings, subFolder):
    slot_duration = settings['settings']['regular']['tsch_slotDuration']

    for k, values in data.items():
        asn_values = [int(item["asn"]) * slot_duration for item in values[0]]
        aoi_values = [int(item["aoi"]) * slot_duration for item in values[0]]

        # Plotting the data
        plt.figure(figsize=(10, 6))
        plt.plot(asn_values, aoi_values, marker='o', linestyle='-', color='blue')

        # Adding labels and title
        plt.xlabel('Time')
        plt.ylabel('AOI')
        plt.title('Time vs AOI')
        plt.grid(True)
        savefig(subFolder, "average_aoi")
        plt.close()

        plot_moving_average(asn_values, aoi_values, subFolder=subFolder,index=0)
        
        # Plotting the peak data only
        plot_peak_data(asn_values, aoi_values,subFolder=subFolder,index=0)
        plot_peaks_moving_average(asn_values, aoi_values,subFolder=subFolder,index=0)

        # plotting the variance of the data
        # I want to see how data is near the min and distribution
        plot_variance_near_min(asn_values, aoi_values,subFolder=subFolder,index=0)

        time.sleep(0.2)  

def plot_moving_average(asn_values, aoi_values, subFolder="", index=0):
    if len(aoi_values) == 0:
        print("Warning: AOI values are empty, skipping moving average plot.")
        return
    
    moving_averages = []

    current_aoi = 0
    aoi_sum = 0
    count = 0
    for index, aoi in enumerate(aoi_values):
        current_aoi = aoi
        aoi_sum += current_aoi
        count += 1

        moving_averages.append((aoi_sum / count))

        if index != len(aoi_values) - 1:
            for i in range(asn_values[index], asn_values[index+1]):
                current_aoi += 1
                aoi_sum += current_aoi
                count += 1
                moving_averages.append((aoi_sum / count))

    asn_values = np.arange(1, len(moving_averages) + 1)

    plt.figure(figsize=(10, 6))
    plt.plot(asn_values, moving_averages, marker='o', linestyle='-', color='green')
    plt.xlabel('Time')
    plt.ylabel('Cumulative Moving Average of AOI')
    plt.title('Cumulative Moving Average of AOI')
    plt.grid(True)
    
    savefig(subFolder, "cumulative_moving_average_aoi")

    plt.close()

def plot_peaks_moving_average(asn_values, aoi_values, subFolder="", index=0):
    if len(aoi_values) == 0:
        print("Warning: AOI values are empty, skipping peak plot.")
        return

    # Find peaks in the AOI values
    peaks, _ = find_peaks(aoi_values)
    peak_values = np.array(aoi_values)[peaks]

    if len(peak_values) == 0:
        print("Warning: No peaks found in AOI values, skipping peak moving average plot.")
        return

    # Calculate the moving average for the peaks
    moving_averages = []
    aoi_sum = 0
    count = 0
    for peak in peak_values:
        aoi_sum += peak
        count += 1
        moving_averages.append(aoi_sum / count)

    peak_asn_values = np.arange(1, len(moving_averages) + 1)

    plt.figure(figsize=(10, 6))
    plt.plot(peak_asn_values, moving_averages, label='Moving Average of Peaks', color='green')

    plt.xlabel('Peak Index')
    plt.ylabel('Moving Average of AOI Peaks')
    plt.title('Moving Average of AOI Peaks')
    plt.legend()
    plt.grid(True)

    savefig(subFolder, "peaks_moving_average_aoi")

    plt.close()

def plot_peak_data(asn_values, aoi_values, subFolder="", index=0):
    if len(aoi_values) == 0:
        print("Warning: AOI values are empty, skipping peak data plot.")
        return

    # Find peaks in the AOI values
    peaks, _ = find_peaks(aoi_values)

    plt.figure(figsize=(10, 6))
    plt.plot(asn_values, aoi_values, label='AOI', color='blue')
    plt.plot(np.array(asn_values)[peaks], np.array(aoi_values)[peaks], "x", label='Peaks', color='red')

    plt.xlabel('Time')
    plt.ylabel('AOI')
    plt.title('AOI Peaks')
    plt.legend()
    plt.grid(True)

    savefig(subFolder, "peak_data_aoi")

    plt.close()

def plot_variance_near_min(asn_values, aoi_values, subFolder="", index=0):
    if len(aoi_values) == 0:
        print("Warning: AOI values are empty, skipping variance plot.")
        return

    # Calculate the global minimum
    global_min = np.min(aoi_values)

    # Calculate the variance of AOI values compared to the global minimum
    variances = [(value - global_min) ** 2 for value in aoi_values]

    plt.figure(figsize=(10, 6))
    plt.plot(asn_values, variances, marker='o', linestyle='-', color='purple')
    plt.xlabel('Time')
    plt.ylabel('Variance from Global Minimum AOI')
    plt.title('Variance from Global Minimum AOI over Time')
    plt.grid(True)

    savefig(subFolder, "variance_from_global_min_aoi_" + str(index))

    plt.close()

def plot_aoi_per_run(data_tuples, subFolder=""):
    if len(data_tuples) == 0:
        print("Warning: Data tuples are empty, skipping moving average plot.")
        return

    # Extracting the index and AOI values from the tuples
    indices = [int(t[0]) for t in data_tuples]  # Extract index values (x-axis)
    aoi_values = [t[1] for t in data_tuples]    # Extract AOI values (y-axis)

    # Compute the moving average with a sliding window
    moving_averages = []
    current_avg_aoi = 0
    aoi_avg_sum = 0
    count = 0
    for index,aoi in enumerate(aoi_values):
        current_avg_aoi = aoi
        aoi_avg_sum += current_avg_aoi
        count += 1
        moving_averages.append((aoi_avg_sum / count))


    aoi_values=np.arange(1,len(moving_averages)+1)
   
    plt.figure(figsize=(10, 6))

    plt.plot(aoi_values, moving_averages, linestyle='-', color='blue')
    
    plt.xlabel('Runs')
    plt.ylabel('Moving Average (AOI values)')
    plt.title('Average of AOI per Run')
    plt.grid(True)    
    # Save the figure if necessary (optional)
    savefig(subFolder, "aoi_per_run")

    plt.close()

def plot_aoi_feedback(data, settings, subFolder):
    # use data in 'asn' , 'aoi_average' and plot it
    for k, values in data.items():
        asn_values = [int(item["asn"]) for item in values[0]]
        aoi_values = [int(item["aoi_average"]) for item in values[0]]

        # Plotting the data
        plt.figure(figsize=(10, 6))
        plt.plot(asn_values, aoi_values, marker='o', color='blue')

        # also add a line indicating the average of the AOI values
        average_aoi = sum(aoi_values) / len(aoi_values)
        plt.axhline(y=average_aoi, color='r', linestyle='--', label='Average AOI')
        # and show value of the average in the line
        plt.text(asn_values[-1], average_aoi, 'Average AOI: {0}'.format(average_aoi))

        # Adding labels and title
        plt.xlabel('Time')
        plt.ylabel('AOI')
        plt.title('AOI Average Calculated in Root')
        plt.grid(True)
        savefig(subFolder, "average_aoi_feedback_in_root")
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
