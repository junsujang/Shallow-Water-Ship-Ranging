# Statistical WI-based Ranging 
The code in this repository implements the statistical range (and waveguide invariant) estimation methd proposed in Jang *et al.* (2025) along with other reference methods. The code is available to those who are interested in implementing the proposed method. Since the data used is not publically available, the codes in this repository will not run properly. One could generate a spectrogram using a simulator, e.g., KRAKEN, to process it.

## Description of each directories
### Acoustic Data
Would contain the acoustic data (both real and simulated) to be processed. The data are not provide in this repository but see [Contacts](#contacts).

### Libraries
Functions that are used to process and perform range estimation based on the acoustic data.

### RealDataProcessing
Scripts to perform range estimation using the acoustic data are here. 

### RealDataProcessing
Scripts to perform WI estimation and analysis (pairwise) using the acoustic data exist here. 

## References
The proposed method is presented in 
- J.Jang, W. H. Hodkiss and F. Meyer, "Ranging of a Moving Ship Using a Single Acoustic Receiver in Shallow Water," submitted to J. Acoust. Soc. Am. (2025)
which can be found in [arxiv link]().

Other references used for implementation are as follows
- K. L. Cockrell and H. Schmidt, “Robust passive range estimation using the waveguide invariant,” J. Acoust. Soc. Am. 127(5), 2780–2789 (2010)
- A. H. Young, H. A. Harms, G. W. Hickman, J. S. Rogers, and J. L. Krolik, “Waveguide-invariant-based ranging and receiver localization using tonal sources of opportunity,” IEEE J. Ocean. Eng. 45(2), 631–644 (2020).
- Y. Yao, C. Sun, X. Liu, and G. Jiang, “Robust range estimation based on the slope of the two-dimensional Fourier transform ridge using radon transform,” IEEE Access 9, 24093–24104 (2021).
- J. Jang and F. Meyer, “Navigation in shallow water using passive acoustic ranging,” in Proc. Int. Conf. Inf. Fusion, Charleston, SC, USA (2023), pp. 1–8.
- J. Jang and F. Meyer, “A new statistical model for waveguide invariant-based range estimation in shallow water,” in 2024 58th Asilomar Conference on Signals, Systems, and Computers (2024), pp. 1–5.
- Y. Le Gall and J. Bonnel, “Passive estimation of the waveguide invariant per pair of modes,” J. Acoust. Soc. Am. 134(2), EL230–EL236 (2013).

### Contacts
To inquire about the real and simulated data used in SBCEX17 or questions about the code, please reach out to Junsu Jang (junsu.jang94@gmail.com).