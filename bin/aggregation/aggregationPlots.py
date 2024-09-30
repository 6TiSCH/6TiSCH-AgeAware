import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import json
import matplotlib.dates as mdates
from matplotlib.dates import DateFormatter
from matplotlib.ticker import MaxNLocator
from collections import OrderedDict

def read_files(options,file_names):
    
    for filename in file_names:

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


        #option.inputfolder is the path to the aggregationKPI directory
        options.inputfolder = os.path.join(bin_folder, 'aggregation/aggregationKPI')
        #after findeing the path directory , it should read from filename
        file_path = os.path.join(options.inputfolder, filename)
        
        with open(os.path.join('E:\\6TiSCH-AgeAware\\bin\\', 'config.json'), 'r') as f:
            settings = json.load(f)
        # init
        data = OrderedDict()

        for key in options.kpis:
            # load data
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
                    break
        # plot_aoi(data,filename)
        plot_moving_average(data,filename,settings)

def plot_aoi(data,label):
    x_values = []
    y_values = []

    for item in data.values():
        x_values.extend([int(entry["asn"]) for entry in item[0]])
        y_values.extend([int(entry["aoi"]) for entry in item[0]])

    # Plot data from the file with a unique color and label
    plt.plot(x_values, y_values, label=label)


def plot_moving_average(data,label,settings, index=0):
    slot_duration = settings['settings']['regular']['tsch_slotDuration']
    aoi_values = []
    asn_values = []

    for item in data.values():
        asn_values.extend([int(entry["asn"]) for entry in item[0]])
        aoi_values.extend([int(entry["aoi"]) for entry in item[0]])
    moving_averages = []

    current_aoi = 0
    aoi_sum = 0
    count = 0
    for index, aoi in enumerate(aoi_values):
        current_aoi = aoi
        aoi_sum += current_aoi
        count += 1

        moving_averages.append((aoi_sum / count) * slot_duration)

        if index != len(aoi_values) - 1:
            for i in range(asn_values[index], asn_values[index+1]):
                current_aoi += 1
                aoi_sum += current_aoi
                count += 1
                moving_averages.append((aoi_sum / count) * slot_duration)

    asn_values = np.arange(1, len(moving_averages) + 1)

    plt.plot(asn_values, moving_averages, label=label)
    

def main():
    file_names=['Stack.dat.kpi','Queue.dat.kpi']

    options = type('', (), {})()
    options.kpis = ['aoi']

    plt.figure()

    read_files(options, file_names)

    # plt.xlabel('ASN')
    # plt.ylabel('AOI')
    # plt.title('AOI vs ASN (Feedback MODE)')
    

    plt.xlabel('Time')
    plt.ylabel('Cumulative Moving Average of AOI')
    plt.title('Cumulative Moving Average of AOI(Feedback MODE)')


    plt.legend()
    plt.grid()

    # plt.savefig('aoi_vs_asn_plot.png')
    plt.savefig('cumulative_moving_average_aoi_plot.png')
    plt.show()

if __name__ == '__main__':
    main()
