from __future__ import division
from __future__ import print_function

# =========================== adjust path =====================================

import os
import sys

import netaddr

if __name__ == '__main__':
    here = sys.path[0]
    sys.path.insert(0, os.path.join(here, '..'))

# ========================== imports ==========================================

import json
import glob
import numpy as np

from SimEngine import SimLog
import SimEngine.Mote.MoteDefines as d

# =========================== defines =========================================

DAGROOT_ID = 0  # we assume first mote is DAGRoot
DAGROOT_IP = 'fd00::1:0'
BATTERY_AA_CAPACITY_mAh = 2821.5

# =========================== decorators ======================================

def openfile(func):
    def inner(inputfile):
        with open(inputfile, 'r') as f:
            return func(f)
    return inner

# =========================== helpers =========================================

def mean(numbers):
    return float(sum(numbers)) / max(len(numbers), 1)

def init_mote():
    return {
        'upstream_num_tx': 0,
        'upstream_num_rx': 0,
        'upstream_num_lost': 0,
        'join_asn': None,
        'join_time_s': None,
        'sync_asn': None,
        'sync_time_s': None,
        'charge_asn': None,
        'upstream_pkts': {},
        'latencies': [],
        'hops': [],
        'charge': None,
        'lifetime_AA_years': None,
        'avg_current_uA': None,
    }

# =========================== KPIs ============================================

@openfile
def kpis_all(inputfile):

    allstats = {} # indexed by run_id, mote_id
    golabi = {}
    packet_received_time = {}

    feedback_addition = {}
    feedback_deletion = {}

    file_settings = json.loads(inputfile.readline())  # first line contains settings
    exec_numSlotframesPerRun = file_settings['exec_numSlotframesPerRun']
    tsch_slotframeLength = file_settings['tsch_slotframeLength']-1
    slot_duration = file_settings['tsch_slotDuration']

    # === gather raw stats

    for line in inputfile:
        logline = json.loads(line)

        # shorthands
        run_id = logline['_run_id']

        if '_asn' in logline: # TODO this should be enforced in each line
            asn = logline['_asn']
        if '_mote_id' in logline: # TODO this should be enforced in each line
            mote_id = logline['_mote_id']

        # populate
        if run_id not in allstats:
            allstats[run_id] = {}
            golabi[run_id] = {}
            packet_received_time[run_id] = {}
            feedback_addition[run_id] = 0
            feedback_deletion[run_id] = 0

        if (
                ('_mote_id' in logline)
                and
                (mote_id not in allstats[run_id])
                # and
                # (mote_id != DAGROOT_ID)
            ):
            allstats[run_id][mote_id] = init_mote()
            golabi[run_id][mote_id] = []
            packet_received_time[run_id][mote_id] = []
        # sync new mote (sync'ed)
        if   logline['_type'] == SimLog.LOG_TSCH_SYNCED['type']:

            # shorthands
            mote_id    = logline['_mote_id']

            # only log non-dagRoot sync times
            if mote_id == DAGROOT_ID:
                continue

            allstats[run_id][mote_id]['sync_asn']  = asn
            allstats[run_id][mote_id]['sync_time_s'] = asn*file_settings['tsch_slotDuration']

        # joined new mote
        elif logline['_type'] == SimLog.LOG_SECJOIN_JOINED['type']:

            # shorthands
            mote_id    = logline['_mote_id']

            # only log non-dagRoot join times
            if mote_id == DAGROOT_ID:
                continue

            # populate
            assert allstats[run_id][mote_id]['sync_asn'] is not None
            allstats[run_id][mote_id]['join_asn']  = asn
            allstats[run_id][mote_id]['join_time_s'] = asn*file_settings['tsch_slotDuration']

        # packet transmission
        elif logline['_type'] == SimLog.LOG_APP_TX['type']:

            # shorthands
            mote_id    = logline['_mote_id']
            dstIp      = logline['packet']['net']['dstIp']
            appcounter = logline['packet']['app']['appcounter']

            # only log upstream packets
            if dstIp != DAGROOT_IP:
                continue

            # populate
            assert allstats[run_id][mote_id]['join_asn'] is not None

            # add this pkt to transfering pkts of this mote
            if appcounter not in allstats[run_id][mote_id]['upstream_pkts']:
                allstats[run_id][mote_id]['upstream_pkts'][appcounter] = {
                    'hops': 0,
                }

            allstats[run_id][mote_id]['upstream_pkts'][appcounter]['tx_asn'] = asn
        # packet reception
        elif logline['_type'] == SimLog.LOG_APP_RX['type']:

            # shorthands
            # find mote_id of sender from srcIp
            mote_id    = netaddr.IPAddress(logline['packet']['net']['srcIp']).words[-1]
            mote_id_reciever =logline['_mote_id']
            dstIp      = logline['packet']['net']['dstIp']
            hop_limit  = logline['packet']['net']['hop_limit']
            appcounter = logline['packet']['app']['appcounter']
            generationTime =logline['packet']['app']['timestamp']
            # only log upstream packets
            if dstIp != DAGROOT_IP:
                # packet_received_time[run_id][0].append({'asn':asn,'tx_asn':generationTime})
                continue

            allstats[run_id][mote_id]['upstream_pkts'][appcounter]['hops']   = (
                d.IPV6_DEFAULT_HOP_LIMIT - hop_limit + 1
            )
            allstats[run_id][mote_id]['upstream_pkts'][appcounter]['rx_asn'] = asn

            golabi[run_id][mote_id].append({'asn':asn, 'tx_asn':generationTime})
            packet_received_time[run_id][mote_id_reciever].append(
            {
                'asn': asn,
                'tx_asn': generationTime,
                'generator_id': mote_id
            })

        # calculate charge usage
        elif logline['_type'] == SimLog.LOG_RADIO_STATS['type']:
            # shorthands
            mote_id    = logline['_mote_id']

            # only log non-dagRoot charge
            if mote_id == DAGROOT_ID:
                continue

            charge =  logline['idle_listen'] * d.CHARGE_IdleListen_uC
            charge += logline['tx_data_rx_ack'] * d.CHARGE_TxDataRxAck_uC
            charge += logline['rx_data_tx_ack'] * d.CHARGE_RxDataTxAck_uC
            charge += logline['tx_data'] * d.CHARGE_TxData_uC
            charge += logline['rx_data'] * d.CHARGE_RxData_uC
            charge += logline['sleep'] * d.CHARGE_Sleep_uC

            allstats[run_id][mote_id]['charge_asn'] = asn
            allstats[run_id][mote_id]['charge']     = charge

        elif logline['_type'] == SimLog.LOG_ASF_ACTION['type']:
            action = logline['action']
            if action == d.SIXP_FEEDBACK_ACTION_ADD:
                feedback_addition[run_id] += 1
            elif action == d.SIXP_FEEDBACK_ACTION_DELETE:
                feedback_deletion[run_id] += 1

    # === aoi stats
   
    aoi_vector = {} 
    for (run_id, per_mote_stats) in list(allstats.items()):
        aoi_vector[run_id] = {}
        mote_id = 0
        aoi_vector[run_id][mote_id] = []

        # hold the most fresh pkt received by each sender
        prev_packets = {}

        for received_packet in packet_received_time[run_id][mote_id]:
            sender_id = received_packet['generator_id']
            prev_packet = prev_packets.get(sender_id, None)

            if prev_packet == None:
                # first packet
                aoi_vector[run_id][mote_id].append({'asn':received_packet['asn'], 'aoi':(received_packet['asn'] - received_packet['tx_asn'])})
                # prev_packet = received_packet
                prev_packets[sender_id] = received_packet
                continue

            if received_packet['tx_asn'] > prev_packet['tx_asn']:
                # up the hill
                aoi_vector[run_id][mote_id].append({'asn':received_packet['asn'], 'aoi':(received_packet['asn'] - prev_packet['tx_asn'])})

                # down the hill
                aoi_vector[run_id][mote_id].append({'asn':received_packet['asn'], 'aoi':(received_packet['asn'] - received_packet['tx_asn'])})

                # prev_packet = received_packet
                prev_packets[sender_id] = received_packet

            # if prev_packet != None:
            #     # up the hill
            #     aoi_vector[run_id][mote_id].append({'asn':item['asn'], 'aoi':(item['asn'] - prev_packet['tx_asn'])})
            #     aoi_stats[run_id][mote_id].append({'asn':item['asn'],'aoi':(item['asn'] - prev_packet['tx_asn'])})
            
            # # down the hill
            # aoi_vector[run_id][mote_id].append({'asn':item['asn'], 'aoi':(item['asn'] - item['tx_asn'])})
            # aoi_stats[run_id][mote_id].append({'asn':item['asn'],'aoi':(item['asn'] - item['tx_asn'])})
            # prev_packet=item  
        allstats[run_id][mote_id]['aoi'] = aoi_vector[run_id][mote_id]

    # === compute advanced motestats

    for (run_id, per_mote_stats) in list(allstats.items()):
        for (mote_id, motestats) in list(per_mote_stats.items()):
            if mote_id != 0:

                # avg_current, lifetime_AA
                if (motestats['sync_asn'] is not None) and (motestats['charge_asn'] is not None):
                    if (
                            (motestats['charge'] <= 0)
                            or
                            (motestats['charge_asn'] <= motestats['sync_asn'])
                        ):
                        motestats['lifetime_AA_years'] = 'N/A'
                    else:
                        motestats['avg_current_uA'] = motestats['charge']/float((motestats['charge_asn']-motestats['sync_asn']) * file_settings['tsch_slotDuration'])
                        assert motestats['avg_current_uA'] > 0
                        motestats['lifetime_AA_years'] = (BATTERY_AA_CAPACITY_mAh*1000/float(motestats['avg_current_uA']))/(24.0*365)

                # latencies, upstream_num_tx, upstream_num_rx, upstream_num_lost
                if motestats['join_asn'] is not None:
                    for (appcounter, pktstats) in list(allstats[run_id][mote_id]['upstream_pkts'].items()):
                        motestats['upstream_num_tx']      += 1
                        # if rx_asn exsists it means pkt was received, otherwise it lost
                        if 'rx_asn' in pktstats:
                            motestats['upstream_num_rx']  += 1
                            thislatency = (pktstats['rx_asn']-pktstats['tx_asn'])*file_settings['tsch_slotDuration']
                            motestats['latencies']  += [thislatency]
                            motestats['hops']       += [pktstats['hops']]
                        else:
                            motestats['upstream_num_lost'] += 1
                    if (motestats['upstream_num_rx'] > 0) and (motestats['upstream_num_tx'] > 0):
                        motestats['latency_min_s'] = min(motestats['latencies'])
                        motestats['latency_avg_s'] = sum(motestats['latencies'])/float(len(motestats['latencies']))
                        motestats['latency_max_s'] = max(motestats['latencies'])
                        motestats['upstream_reliability'] = motestats['upstream_num_rx']/float(motestats['upstream_num_tx'])
                        motestats['avg_hops'] = sum(motestats['hops'])/float(len(motestats['hops']))

    # === network stats
    for (run_id, per_mote_stats) in list(allstats.items()):

        #-- define stats

        app_packets_sent = 0
        app_packets_received = 0
        app_packets_lost = 0
        joining_times = []
        us_latencies = []
        current_consumed = []
        lifetimes = []
        average_aoi_seconds = 0

        #-- compute stats

        for (mote_id, motestats) in list(per_mote_stats.items()):
            if mote_id == DAGROOT_ID:
                continue

            # counters

            app_packets_sent += motestats['upstream_num_tx']
            app_packets_received += motestats['upstream_num_rx']
            app_packets_lost += motestats['upstream_num_lost']

            # joining times

            if motestats['join_asn'] is not None:
                joining_times.append(motestats['join_asn'])

            # latency

            us_latencies += motestats['latencies']

            # current consumed

            current_consumed.append(motestats['charge'])
            if motestats['lifetime_AA_years'] is not None:
                lifetimes.append(motestats['lifetime_AA_years'])
            current_consumed = [
                value for value in current_consumed if value is not None
            ]

            # average values for aoi in every asns
            current_aoi = 0
            aoi_sum = 0
            asn_count = 0
            for index, event in enumerate(aoi_vector[run_id][0]):
                current_aoi = event['aoi']
                aoi_sum += current_aoi
                asn_count += 1

                if index != len(aoi_vector[run_id][0]) - 1:
                    for i in range(event['asn'], aoi_vector[run_id][0][index + 1]['asn']):
                        current_aoi += 1
                        aoi_sum += current_aoi
                        asn_count += 1

            average_aoi_seconds = (aoi_sum / asn_count) * slot_duration
            average_aoi_asn = aoi_sum / asn_count

            # variance of aoi
            variance_aoi = 0
            variance_count = 0
            min_aoi = min([event['aoi'] for event in aoi_vector[run_id][0]]) 
            for index, event in enumerate(aoi_vector[run_id][0]):
                variance_aoi += (event['aoi'] - min_aoi) ** 2
                variance_count += 1

                if index != len(aoi_vector[run_id][0]) - 1:
                    for i in range(event['asn'], aoi_vector[run_id][0][index+1]['asn']):
                        variance_aoi += (current_aoi - min_aoi) ** 2
                        variance_count += 1

            variance_aoi = (variance_aoi / variance_count) * slot_duration

        #-- save stats

        allstats[run_id]['global-stats'] = {
            'e2e-upstream-delivery': [
                {
                    'name': 'E2E Upstream Delivery Ratio',
                    'unit': '%',
                    'value': (
                        1 - app_packets_lost / app_packets_sent
                        if app_packets_sent > 0 else 'N/A'
                    )
                },
                {
                    'name': 'E2E Upstream Loss Rate',
                    'unit': '%',
                    'value': (
                        app_packets_lost / app_packets_sent
                        if app_packets_sent > 0 else 'N/A'
                    )
                }
            ],
            'e2e-upstream-latency': [
                {
                    'name': 'E2E Upstream Latency',
                    'unit': 's',
                    'mean': (
                        mean(us_latencies)
                        if us_latencies else 'N/A'
                    ),
                    'min': (
                        min(us_latencies)
                        if us_latencies else 'N/A'
                    ),
                    'max': (
                        max(us_latencies)
                        if us_latencies else 'N/A'
                    ),
                    '99%': (
                        np.percentile(us_latencies, 99)
                        if us_latencies else 'N/A'
                    )
                },
                {
                    'name': 'E2E Upstream Latency',
                    'unit': 'slots',
                    'mean': (
                        mean(us_latencies) / slot_duration
                        if us_latencies else 'N/A'
                    ),
                    'min': (
                        min(us_latencies) / slot_duration
                        if us_latencies else 'N/A'
                    ),
                    'max': (
                        max(us_latencies) / slot_duration
                        if us_latencies else 'N/A'
                    ),
                    '99%': (
                        np.percentile(us_latencies, 99) / slot_duration
                        if us_latencies else 'N/A'
                    )
                }
            ],
            'current-consumed': [
                {
                    'name': 'Current Consumed',
                    'unit': 'mA',
                    'mean': (
                        mean(current_consumed)
                        if current_consumed else 'N/A'
                    ),
                    '99%': (
                        np.percentile(current_consumed, 99)
                        if current_consumed else 'N/A'
                    )
                }
            ],
            'network_lifetime':[
                {
                    'name': 'Network Lifetime',
                    'unit': 'years',
                    'min': (
                        min(lifetimes)
                        if lifetimes else 'N/A'
                    ),
                    'total_capacity_mAh': BATTERY_AA_CAPACITY_mAh,
                }
            ],
            'joining-time': [
                {
                    'name': 'Joining Time',
                    'unit': 'slots',
                    'min': (
                        min(joining_times)
                        if joining_times else 'N/A'
                    ),
                    'max': (
                        max(joining_times)
                        if joining_times else 'N/A'
                    ),
                    'mean': (
                        mean(joining_times)
                        if joining_times else 'N/A'
                    ),
                    '99%': (
                        np.percentile(joining_times, 99)
                        if joining_times else 'N/A'
                    )
                }
            ],
            'app-packets-sent': [
                {
                    'name': 'Number of application packets sent',
                    'total': app_packets_sent
                }
            ],
            'app_packets_received': [
                {
                    'name': 'Number of application packets received',
                    'total': app_packets_received
                }
            ],
            'app_packets_lost': [
                {
                    'name': 'Number of application packets lost',
                    'total': app_packets_lost
                }
            ],
            'aoi_stats': [
                {
                    'name': 'Average Age of Information',
                    'unit': 's',
                    'value': average_aoi_seconds
                },
                {
                    'name': 'Average Age of Information',
                    'unit': 'asn',
                    'value': average_aoi_asn
                },
                {
                    'name': 'Variance of Age of Information',
                    'unit': 's',
                    'value': variance_aoi
                }
            ],
            'aoi_feedbacks': [
                {
                    'name': 'Number of Addition Requests',
                    'unit': 'Packets',
                    'value': feedback_addition[run_id]
                },
                {
                    'name': 'Number of Deletion Requests',
                    'unit': 'Packets',
                    'value': feedback_deletion[run_id]
                }
            ]
        }

    # === remove unnecessary stats

    for (run_id, per_mote_stats) in list(allstats.items()):
        for (mote_id, motestats) in list(per_mote_stats.items()):
            if 'sync_asn' in motestats:
                del motestats['sync_asn']
            if 'charge_asn' in motestats:
                del motestats['charge_asn']
                del motestats['charge']
            if 'join_asn' in motestats:
                del motestats['upstream_pkts']
                del motestats['hops']
                del motestats['join_asn']

    return allstats

# =========================== main ============================================

def main():
    #inputfolder = os.listdir('simData')

    # hard-coded config when debugging
    #inputfolder = 'C:/Arshia/6TiCSH-AgeAware/bin/simData'

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
    inputfolder = os.path.join(bin_folder, 'simData')

    # FIXME: This logic could be a helper method for other scripts
    # Identify simData having the latest results. That directory should have
    # the latest "mtime".
    subfolders = list(
        [os.path.join(inputfolder, x) for x in os.listdir(inputfolder)]
    )
    subfolder = max(subfolders, key=os.path.getmtime)

    for infile in glob.glob(os.path.join(subfolder, '*.dat')):
        print('generating KPIs for {0}'.format(infile))

        # gather the kpis
        kpis = kpis_all(infile)

        # print on the terminal
        print(json.dumps(kpis, indent=4))

        # add to the data folder
        outfile = '{0}.kpi'.format(infile)
        with open(outfile, 'w') as f:
            f.write(json.dumps(kpis, indent=4))
        print('KPIs saved in {0}'.format(outfile))

if __name__ == '__main__':
    main()
