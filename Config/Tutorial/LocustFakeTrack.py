#!/usr/bin/env python

'''
Script to run Locust fake track simulation many times. Each time it generates two waterfall root files: A fake track with and without noise
Run: python LocustFakeTrack.py n_sims [-h] [-w WORKING_DIR] [-l LOCUST_BIN] [-k KATYDID_BIN] [-c CONFIG]
                          

Author: L. Saldana
Date: 10/30/2018
'''

import argparse
import json
import random
import subprocess
import os
import sys

def RunSimulation(n_sims,working_dir,locust_binary_path,katydid_binary_path,katydid_config_path):

    # Generator lists
    generators = ['fake-track','lpf-fft','decimate-signal','digitizer']
    generators_wnoise = ['fake-track','lpf-fft','decimate-signal','gaussian-noise','digitizer']

    # Generator properties
    fake_track_gen = {'lo-frequency': 20e9,
                      'ntracks-mean': 2.0,
                      'signal-power': 1e-15,
                      'slope-mean': 0.6,
                      'slope-std': 0.025,
                      'start-frequency-max': 20.053e9,
                      'start-frequency-min': 20.049e9,
                      'start-time-max': 0.003,
                      'start-time-min': 0.001,
                      'start-vphase': 0.0,
                      'track-length-mean': 0.001}
    gaussian_noise_gen = {'domain': 'time',
                          'noise-floor': 4e-21}
    simulation_gen = {'n-channels': 1,
                      'n-records': 1,
                      'record-size': 4194304}
    simulation_wnoise_gen = simulation_gen.copy()
    digitizer_gen = {'v-offset': -0.5e-05,
                     'v-range': 1e-05}

    # Generate Locust JSON config files
    faketrack_config = {'generators': generators,
                        'digitizer': digitizer_gen}
    faketrack_wnoise_config = {'generators': generators_wnoise,
                               'gaussian-noise': gaussian_noise_gen,
                               'digitizer': digitizer_gen}

    # Run Simulations
    for ii_sim in range(n_sims):

        rand_seed = random.randint(540559518,1325542009) # choose random seed

        # Setting up JSON configs
        fake_track_gen['random-seed'] = rand_seed
        locust_egg_file = working_dir + 'locust_faketrack' + '_{}'.format(ii_sim) + '.egg'
        locust_egg_wnoise_file = working_dir + 'locust_faketrack_wnoise' + '_{}'.format(ii_sim) + '.egg'
        waterfall_file = working_dir + 'locust_faketrack' + '_{}'.format(ii_sim) + '.root'
        waterfall_wnoise_file = working_dir + 'locust_faketrack_wnoise' + '_{}'.format(ii_sim) + '.root'
        simulation_gen['egg-filename'] = locust_egg_file
        simulation_wnoise_gen['egg-filename'] = locust_egg_wnoise_file
        faketrack_config['fake-track'] = fake_track_gen
        faketrack_wnoise_config['fake-track'] = fake_track_gen
        faketrack_config['simulation'] = simulation_gen
        faketrack_wnoise_config['simulation'] = simulation_wnoise_gen
        faketrack_config_json = json.dumps(faketrack_config,indent=4,sort_keys=True,separators=(',',':'))
        faketrack_wnoise_config_json = json.dumps(faketrack_wnoise_config,indent=4,sort_keys=True,separators=(',',':'))
        locust_config_path = '{}LocustFakeTrack_{}.json'.format(working_dir,ii_sim)
        locust_wnoise_config_path = '{}LocustFakeTrack_wnoise_{}.json'.format(working_dir,ii_sim)
        
        # Creating JSON file, run simulation, process with Katydid; no noise
        locust_config = open(locust_config_path,'a')
        locust_config.write(faketrack_config_json)
        locust_config.close()
        print('Created: {}'.format(locust_config_path))
        print('\tRunning simultion without noise')
        try: # Locust
            output = subprocess.check_output("{} config={}".format(locust_binary_path,locust_config_path), shell=True, stderr=subprocess.STDOUT)
            os.remove(locust_config_path)  
            print('\t\tSucessful! Removed JSON file') 
        except subprocess.CalledProcessError as e:
            print("Error: {}".format(e.output))
            break
        try: # Katydid
            print('\tRunning Katydid...')
            output = subprocess.check_output("{} -c {} -e {} --waterfall-writer.output-file={}".format(katydid_binary_path,katydid_config_path,locust_egg_file,waterfall_file),shell=True,stderr=subprocess.STDOUT)
            print('\t\tSucessfully created: {}'.format(waterfall_file))
        except subprocess.CalledProcessError as e:
            print("Error: {}".format(e.output))
            break

        # Creating JSON file, run simulation, process with Katydid; with noise
        locust_wnoise_config = open(locust_wnoise_config_path,'a')
        locust_wnoise_config.write(faketrack_wnoise_config_json)
        locust_wnoise_config.close()
        print('Created: {}'.format(locust_wnoise_config_path)) 
        print('\tRunning simultion with noise')
        try: # Locust
            output_wnoise = subprocess.check_output("{} config={}".format(locust_binary_path,locust_wnoise_config_path), shell=True, stderr=subprocess.STDOUT)
            os.remove(locust_wnoise_config_path)  
            print('\t\tSucessful! Removed JSON file') 
        except subprocess.CalledProcessError as e:
            print("Error: {}".format(e.output))
            break
        try: # Katydid
            print('\tRunning Katydid...')
            output = subprocess.check_output("{} -c {} -e {} --waterfall-writer.output-file={}".format(katydid_binary_path,katydid_config_path,locust_egg_wnoise_file,waterfall_wnoise_file),shell=True,stderr=subprocess.STDOUT)
            print('\t\tSucessfully created: {}'.format(waterfall_wnoise_file))
        except subprocess.CalledProcessError as e:
            print("Error: {}".format(e.output)) 
            break

        return

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="generate fake tracks in Locust and obtain Katydid waterfall spectrograms")

    parser.add_argument('n_sims',
                        help='Number of simulations',
                        type=int) 
    parser.add_argument('-w','--working_dir',
                        help="Path to working directory to save egg and root files",
                        type=str,
                        default='/home/les67/temp/')
    parser.add_argument('-l','--locust_bin',
                        help="Path to Locust binary",
                        type=str,
                        default='/home/les67/locust_mc/build/bin/LocustSim')
    parser.add_argument('-k','--katydid_bin',
                        help="Path to Katydid binary",
                        type=str,
                        default='/home/les67/katydid/build/bin/Katydid')
    parser.add_argument('-c','--config',
                        help="Path to Katydid config file",
                        type=str,
                        default='/home/les67/locust_mc/Config/Tutorial/katydid_faketrack.json')     

    args = parser.parse_args()
    if not args.working_dir.endswith(os.path.sep): # make sure working dir path is correctly formatted
        args.working_dir += os.path.sep

    print('--------RUNNING {} FAKE TRACK SIMULATIONS!--------\n'.format(args.n_sims))
    RunSimulation(args.n_sims,args.working_dir,args.locust_bin,args.katydid_bin,args.config)    
    print('--------SIMULATION DONE--------')

    sys.exit(0)
