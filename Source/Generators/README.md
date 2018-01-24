Directory Structure
-------------------
*   Contains headers/ source code for the generators for the receiver chain. As indicated, LMCKassSignalGenerator and LMCFreeFieldGenerator are the signal generators, using the particle simulated in kassiopeia. These generators are appropriate for Phases II and III, respectively.
*  LMCButterworthLPFGenerator and LMCLowPassFilterGenerator are the low pass filters. Presently the butterworth filter is not optimal, but still functioning.
*  LMCGaussianNoiseGenerator generates white noise in either the frequency or time domains. 
