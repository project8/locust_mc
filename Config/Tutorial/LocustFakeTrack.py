import json

# Generator list
generators = ['fake-track','lpf-fft','decimate-signal','gaussian-noise','digitizer']

# Generator properties
fake_track_gen = {'lo-frequency': 20e9,
                  'ntracks-mean': 2.0,
                  'random-seed': 0,
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
simulation_gen = {'egg-filename': '/home/les67/locust_faketrack.egg',
                  'n-channels': 1,
                  'n-records': 1,
                  'record-size': 4194304}
digitizer_gen = {'v-offset': -5e-06,
                 'v-range': 1e-05}

# Generate Locust JSON config file
faketrack_config = {'generators': generators,
                    'fake_track': fake_track_gen,
                    'simulation': simulation_gen,
                    'gaussian-noise': gaussian_noise_gen,
                    'digitizer': digitizer_gen}
faketrack_config_json = json.dumps(faketrack_config,indent=4,sort_keys=True,separators=(',',':'))
config_file = open('/home/les67/locust_mc/Config/Tutorial/LocustFakeTrack_temp.json','a')
config_file.write(faketrack_config_json)
config_file.close()
