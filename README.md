locust_mc
=========

Monte Carlo including correlation between channels

Required libraries:
yajl
monarch
fftw

Required config parameters

waveguide_setup
if DOUBLEAMP
  receiver1_noise_temperature
  receiver2_noise_temperature
  amp1_noise_temperature
  amp2_noise_temperature
  waveguide_length
  phase_delay_length
if SINGELAMP
  waveguide_length
  amp1_noise_temperature
  receiver1_noise_temperature
  distance_to_short
