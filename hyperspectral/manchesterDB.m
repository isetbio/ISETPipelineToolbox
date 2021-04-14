function [allData, nScene] = manchesterDB()

hyperspectralRootDir = '/resources/scenes/hyperspectral';
rdt = RdtClient('isetbio');
rdt.crp(hyperspectralRootDir);

% Select a database
dataBaseName = 'manchester_database';
nScene = 8;

allData = struct();
allData.image = {}; allData.wave = {};

count = 1;
for family = {'2002', '2004'}
    rdt.crp(fullfile(hyperspectralRootDir, dataBaseName, family{1}));
    for idx = 1 : nScene
        sceneName = strcat('scene', num2str(idx), '_', family{1}, '_iset');
        data = rdt.readArtifact(sceneName);
        
        scene = sceneFromBasis(data);
        
        allData.image{count} = scene.data.photons;
        allData.wave{count}  = scene.spectrum.wave;
        count = count + 1;
    end
end

% Use the same wavelength band for all images
refWave = 420 : 20 : 700;
for idx = 1 : (count - 1)
    image = allData.image{idx};
    wave  = allData.wave{idx};
    
    [h, w, d] = size(image);
    image = reshape(image, h * w, d);
    image = interp1(wave, image', refWave);
    
    assert(sum(isnan(image(:))) == 0);
    
    allData.image{idx} = reshape(image', [h, w, length(refWave)]);
end

allData.wave = refWave;
nScene = count - 1;

end


