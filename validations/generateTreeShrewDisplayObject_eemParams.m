function presentationDisplay = generateTreeShrewDisplayObject_eemParams()
    % Specify files containing the SPDs of the R,G,B guns
    spdDataFile = 'demoSPD.mat';

    % Specify file containing the ambient SPD (when R=G=B=0) 
    ambientSPDDataFile = 'demoAmbientSPD.mat';

    % Generate fake gamma table - or import the actual gamma table
    gammaTableLength = 256;
    gammaExponent = 2.0;
    gammaTable = repmat((linspace(0,1,gammaTableLength)').^gammaExponent, [1 3]);

    % Specify display resolution - UPDATED
    % Screen resolution: 1920x1200
    % Screen dimensions: 152x94.5mm
    dotsPerInch = 322;

    % Set viewing distance between display and animal - UPDATED
    viewingDistanceMeters = 0.08;

    % Visualize the SPDs and the gamma table
    plotDisplayCharacteristics = false;

    % Load the RGB gun spectral power distributions
    fprintf('Loading SPDs from %s\n', spdDataFile);
    load(spdDataFile, 'spd');

    % Check data consistency
    assert(size(spd, 2) == 4, 'The rgbSPD matrix must be an N x 4 matrix, with the first column being the spectral support');
    wavelengthSupport = spd(:,1);
    rgbSPDs = spd(:,2:4);

    % Load the ambient SPD
    fprintf('Loading ambient SPD from %s\n', ambientSPDDataFile);
    load(ambientSPDDataFile, 'spd');
    assert(size(spd, 2) == 2, 'The ambientSPD matrix must be an N x 2 matrix, with the first column being the spectral support');
    ambientWavelengthSupport = spd(:,1);
    ambientSPD = spd(:,2);

    size(rgbSPDs);
    size(ambientSPD);
    % Check data consistency
    assert(size(rgbSPDs,1) == size(ambientSPD,1), 'The ambient SPD must have the same wavelength entries as the display SPD');
    assert(all(ambientWavelengthSupport == wavelengthSupport), 'The ambient wavelength support must match the wavelength support of the display SPD');
     
    presentationDisplay = generateCustomDisplay(...
           'dotsPerInch', dotsPerInch, ...
           'wavelengthSupportNanoMeters', wavelengthSupport, ...
           'spectralPowerDistributionWattsPerSteradianM2NanoMeter', rgbSPDs, ...
           'ambientSPDWattsPerSteradianM2NanoMeter', ambientSPD, ...
           'gammaTable', gammaTable, ...
           'plotCharacteristics', plotDisplayCharacteristics);
    
    % Set the viewing distance
    presentationDisplay = displaySet(presentationDisplay, ...
        'viewing distance', viewingDistanceMeters );
    
end