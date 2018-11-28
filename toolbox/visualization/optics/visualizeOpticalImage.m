function visualizeOpticalImage(opticalImage)
  
% retrieve the spatial support of the scene(in millimeters)
spatialSupportMM = oiGet(opticalImage, 'spatial support', 'mm');
% Convert spatial support in degrees
optics = oiGet(opticalImage, 'optics');
focalLength = opticsGet(optics, 'focal length');
mmPerDegree = focalLength*tand(1)*1e3;
spatialSupportDegs = spatialSupportMM/mmPerDegree;

% retrieve the sRGB components of the optical image (just for visualization)
rgbImage = oiGet(opticalImage, 'rgb image');
% retrieve the XYZ tristimulus components of the optical image
XYZmap = oiGet(opticalImage, 'xyz');
% compute the xy-chroma and the luminance maps
luminanceMap = squeeze(XYZmap(:,:,2));
xMap = squeeze(XYZmap(:,:,1))./sum(XYZmap,3);
yMap = squeeze(XYZmap(:,:,2))./sum(XYZmap,3);
% Compute mean chromaticity
meanChromaticity = [mean(xMap(:)) mean(yMap(:))];
visualizeSceneRGB(spatialSupportDegs, 'degs', rgbImage, ...
    [], meanChromaticity, 'Gabor optical image');

 
% retrieve the retinal irradiance (photon rate in photons/pixel/sec/nm)
retinalImagePhotonRate = oiGet(opticalImage, 'photons');
% retrieve wavelength support
wavelengthSupport = oiGet(opticalImage, 'wave');
wavelengthBandsToVisualize = 400:30:700;
%
visualizeSceneRadiance(spatialSupportDegs, 'degs', ...
    retinalImagePhotonRate, wavelengthSupport, wavelengthBandsToVisualize);
% % visualize the optical image as RGB

end
