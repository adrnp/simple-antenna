function antenna = convertToMatlabElement(ant)

% get the full pattern from the simpant antenna
patterndB = ant.PatterndB;

% convert the thetoetical pattern to the frame for matlab
[patchPatternAzEl, az, el] = phitheta2azelpat(patterndB', unique(ant.PhiMesh(:)), unique(ant.ThetaMesh(:)));

% create custom antenna element to work with the toolbox framework
freqVector = ant.Frequency*1e9 * [0.9 1.1];  % Frequency range for element pattern
antenna = phased.CustomAntennaElement('FrequencyVector', freqVector, ...
                              'AzimuthAngles', az, ...
                              'ElevationAngles', el, ...
                              'MagnitudePattern', patchPatternAzEl, ...
                              'PhasePattern', zeros(size(patchPatternAzEl)));