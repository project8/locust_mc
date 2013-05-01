import json

def generate_configs(file_prefix=False, num_events=100, duration=1e-4, dfdt=1e8, power=1e-15):
    '''
        generate a config json file for locust input

        <num_events> is number of chirps
        <duration> is time length of each chirp
        <dfdt> is the chirp rate
        <power> is the signal power
    '''
    if not file_prefix:
        file_prefix = (str(num_events) + 'events_' + str(duration) + 'dur_' +
            str(dfdt) + 'dfdt_' + str(power) + 'power')
    outfile = open(file_prefix + '_config.json', 'w')
    out_dict = {
        "amp1_noise_temperature": 30,
        "amp2_noise_temperature": 30,
        "receiver1_noise_temperature": 30,
        "receiver2_noise_temperature": 30,
        "BField": 0.9,
        "waveguide_length": 20,
        "phase_delay_length": 5,
        "hf_mixing_frequency": 24.2e9,
        "lf_mixing_frequency": 500e6,
        "datafile_duration": 0.02,
        "egg_outfile_name": file_prefix + ".egg",
        "mcinfo_outfile_name": file_prefix + ".mcinfo",
        "number_of_events": num_events,
        "events": []
    }
    for start in range(1, num_events+1):
        out_dict['events'].append({
            "start_time": start*1e-4,
            "duration": duration,
            "start_frequency": 24.750e9,
            "dfdt": dfdt,
            "power": power
            })
    print(json.dumps(out_dict, indent=4))
    json.dump(out_dict, outfile, indent=4)
