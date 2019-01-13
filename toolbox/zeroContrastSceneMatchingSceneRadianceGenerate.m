function zeroContrastScene = zeroContrastSceneMatchingSceneRadianceGenerate(sourceScene, visualizeScenes)
% Generate a uniform scene whose radiance matches that of the source scene
%
% Syntax:
%   zeroContrastScene = ...
%    zeroContrastSceneMatchingSceneRadianceGenerate(sourceScene, visualizeScenes);
%
% Description:
%    Generate a spatially-uniform scene whose radiance (photon rate)
%    matches that of the sourceScene. The spatial extent and all other
%    params of the generated scene match the sourceScene  
%
% Inputs:
%    sourceScene                    ISETBio scene (source)
%    visualizeScenes                Boolean, whether to visualize the
%                                   source and generated scene
%
% Outputs:
%    zeroContrastScene              The generated zero contrast scene
%

% History:
%    1/24/2019  NPC   Wrote it

    photons = sceneGet(sourceScene, 'photons');
    meanPhotons = bsxfun(@plus, photons*0, mean(mean(photons,1),2));
    zeroContrastScene = sourceScene;
    zeroContrastScene = sceneSet(zeroContrastScene, 'photons', meanPhotons);
    
    if (visualizeScenes)
        figure();
        subplot(1,2,1)
        imagesc(sceneGet(sourceScene, 'rgb image'));
        axis image;
        subplot(1,2,2)
        imagesc(sceneGet(zeroContrastScene, 'rgb image'));
        axis image;
    end
end

