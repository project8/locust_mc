Directory Structure
-------------------
*  Contains headers/ source code for the generators for the receiver chain. The "SignalGenerators" actually generate the actual signals. 

*  LMCFakeTrackSignalGenerator generates a "fake" track with a specified slope, start frequency, start time, bypassing the kassiopeia interface
*  LMCKassSignalGenerator is the Phase I/ Phase II signal generator (configurable), using kassiopeia as a particle tracker.
*  LMCPatchSignalGenerator and LMCFreeFieldGenerator are the Phase III signal generators, using particles simulated in kassiopeia. 
*  LMCButterworthLPFGenerator and LMCLowPassFilterGenerator are the low pass filters. 
*  LMCDecimateSignalGenerator keeps every 10th entry in time (see DSP decimation), limiting the bandwidth post-filtering.
*  LMCGaussianNoiseGenerator generates white noise in either the frequency or time domains. 
*  LMCTemplateGenerator is not functional, gives a template from which to create a new generator.
