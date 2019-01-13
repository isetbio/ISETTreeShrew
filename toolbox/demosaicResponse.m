function dStruct = demosaicResponse(theMosaicResponse, theConeMosaic, demosaicLatencySeconds, imageSizePixels, interpolationMethod)
    
    dStruct = struct(...
        'theDemosaicedLconeIsomerizationMap', [], ...
        'theDemosaicedMconeIsomerizationMap', [], ...
        'theDemosaicedSconeIsomerizationMap', [], ...
        'demosaicedMapSupportDegs', [], ...
        'theLConeIsomerizations', [], ...
        'theLConeXlocsDegs', [], ...
        'theLConeYlocsDegs', [], ...
        'theMConeIsomerizations', [], ...
        'theMConeXlocsDegs', [], ...
        'theMConeYlocsDegs', [], ...
        'theSConeIsomerizations', [], ...
        'theSConeXlocsDegs', [], ...
        'theSConeYlocsDegs', []);
    
        timeAxis = theConeMosaic.timeAxis;
        timeAxis = timeAxis - mean(timeAxis);
        [~, centerTimeBin] = min(abs(timeAxis-demosaicLatencySeconds));
        %fprintf('Demosaicing mean (across trials) response at t = %2.1f ms\n', timeAxis(centerTimeBin)*1000);
        
        %[nTrials, nConesNum, nTimeBins] = size(theMosaicResponse)
        
        % Demosaic the mean response across all trials
        meanResponse = theConeMosaic.reshapeHex1DmapToHex2Dmap(squeeze(mean(theMosaicResponse(:,:,centerTimeBin),1)));
        
        demosaicingSampleSpacingMicrons = 1;
        
        % L-cone interpolation
        [dStruct.theDemosaicedLconeIsomerizationMap, ...
         dStruct.demosaicedMapSupportDegs, ...
         dStruct.theLConeIsomerizations, ...
         dStruct.theLConeXlocsDegs, ...
         dStruct.theLConeYlocsDegs] = ...
         theConeMosaic.demosaicConeTypeActivationFromFullActivation(...
            'L-cones', meanResponse, demosaicingSampleSpacingMicrons, ...
            'interpolationMethod', interpolationMethod); 
        
        
        % S-cone interpolation
        [dStruct.theDemosaicedSconeIsomerizationMap, ...
         dStruct.demosaicedMapSupportDegs, ...
         dStruct.theSConeIsomerizations, ...
         dStruct.theSConeXlocsDegs, ...
         dStruct.theSConeYlocsDegs] = ...
         theConeMosaic.demosaicConeTypeActivationFromFullActivation(...
            'S-cones', meanResponse, demosaicingSampleSpacingMicrons, ...
            'interpolationMethod', interpolationMethod); 
        
        % Resample to imageSizePixels (the input image size)
        [X,Y] = meshgrid(dStruct.demosaicedMapSupportDegs,dStruct.demosaicedMapSupportDegs);
        
        dStruct.demosaicedMapSupportDegs = ...
            linspace(dStruct.demosaicedMapSupportDegs(1),dStruct.demosaicedMapSupportDegs(end), imageSizePixels(1));
        
        [Xq,Yq] = meshgrid(dStruct.demosaicedMapSupportDegs,dStruct.demosaicedMapSupportDegs);
        dStruct.theDemosaicedLconeIsomerizationMap = interp2(X,Y,dStruct.theDemosaicedLconeIsomerizationMap,Xq,Yq, 'linear');

        [Xq,Yq] = meshgrid(dStruct.demosaicedMapSupportDegs,dStruct.demosaicedMapSupportDegs);
        dStruct.theDemosaicedSconeIsomerizationMap = interp2(X,Y,dStruct.theDemosaicedSconeIsomerizationMap,Xq,Yq, 'linear');
end

