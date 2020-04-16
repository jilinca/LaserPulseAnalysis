Takes in raw measured intensity data. The functions are looking for occurences 
of pulses in the one dimensionsal data. Functions can find the waveform of the specific 
pulse and recontruct it via different fitting functions. Finds the peak value
of the pulse from the fitted pulses to provide a more accurate measurement result
for the your measuring device. This code is primarily meant for assecing RAW lidar
waveform data.


The Levenberg Marquardt fitting procedure depends on 
https://github.com/pngts/Nonlinear-parameter-estimation-in-thermodynamic-models/blob/master/LevenbergMarquardt.m
https://www.mathworks.com/matlabcentral/fileexchange/53440-jacobian-toolbox
