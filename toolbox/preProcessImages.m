function [inputImageArray, sceneArray] = preProcessImages(imDir, imFileNames, ...
    presentationDisplay, resampledImageSize, ...
    resampledImageSizeVisualAngle, meanLuminanceCdPerM2)
    
    
    mRows = resampledImageSize(1);
    nCols = resampledImageSize(2);
    
    resampledImage = zeros(mRows, nCols,3);
    xDestRange = (0:(nCols-1))/nCols;
    yDestRange = (0:(mRows-1))/mRows;
    [Xdest, Ydest] = meshgrid(xDestRange, yDestRange);
    
    % Initialize outputs
    inputImageArray = cell(1, numel(imFileNames));
    sceneArray = cell(1, numel(imFileNames));
    
    for imIndex = 1:numel(imFileNames)
        imFileName = sprintf('%s.jpg',imFileNames{imIndex});
        rawSourceImage = double(imread(fullfile(imDir,imFileName)));
        nCols = size(rawSourceImage,2);
        mRows = size(rawSourceImage,1);
        if (ndims(rawSourceImage) == 2)
           rawSourceImage = repmat(rawSourceImage, [1 1 3]);
        end
        minRGB = min(rawSourceImage(:));
        maxRGB = max(rawSourceImage(:));
        rawSourceImage = (rawSourceImage-minRGB)/(maxRGB-minRGB);

        % save input RGB image
        inputImageArray{imIndex} = rawSourceImage;
        
        
        xSourceRange = (0:(nCols-1))/nCols;
        ySourceRange = (0:(mRows-1))/mRows;
        [X, Y] = meshgrid(xSourceRange, ySourceRange);
        for channelIndex = 1:size(rawSourceImage,3)
            resampledImage(:,:,channelIndex) = interp2(X,Y,...
                squeeze(rawSourceImage(:,:,channelIndex)), ...
                Xdest, Ydest, 'linear');
        end
        
        % Generate ISETBio scene from resampled input RGB image
        realizedScene = sceneFromFile(resampledImage, 'rgb', ...
            meanLuminanceCdPerM2, presentationDisplay);
        realizedScene = sceneSet(realizedScene, 'fov', resampledImageSizeVisualAngle);
        realizedScene = sceneSet(realizedScene, 'name', imFileName);
        
        % Save output scene
        sceneArray{imIndex} = realizedScene;
        
        if (1==2)
        figure(1); clf;
        subplot(1,3,1);
        size(inputImageArray{imIndex})
        image(1:nCols, 1:mRows, inputImageArray{imIndex}); axis 'image';
        title('input RGB');

        subplot(1,3,2);
        image(1:resampledImageSize(1), 1:resampledImageSize(2), resampledImage); axis 'image';
        title(sprintf('input RGB - resampled [%d x %d]', resampledImageSize(1), resampledImageSize(2)));
        
       
        
        realizedSceneRGB = sceneGet(sceneArray{imIndex}, 'rgb');
        spatialSupportMilliMeters = sceneGet(sceneArray{imIndex}, 'spatial support', 'mm');
        spatialSupportX = squeeze(spatialSupportMilliMeters(1,:,1));
        spatialSupportY = squeeze(spatialSupportMilliMeters(:,1,2));
        
        subplot(1,3,3);
        image(spatialSupportX, spatialSupportY, realizedSceneRGB); axis 'image';
        xlabel('space (mm on presentation display)');
        title('ISETBio scene');
        pause
        end
    end
    
end

