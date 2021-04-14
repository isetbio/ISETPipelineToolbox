function [allData, nScene] = harvardDB()

hyperspectralRootDir = '/resources/scenes/hyperspectral';
rdt = RdtClient('isetbio');
rdt.crp(hyperspectralRootDir);

% Select a database
dataBaseName = 'harvard_database';

prefix = {'', 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h'};
number = {[1, 6], [1, 8], [1, 9], [1, 9], [0, 9], ...
    [0, 7], [1, 8], [0, 9], [0, 7]};

% Select a scene family
rdt.crp(fullfile(hyperspectralRootDir, dataBaseName));

allData = struct();
allData.image = {}; allData.wave = {};

count = 1;
for idx = 1 : length(prefix)
    bound = number{idx};
    fprintf('prefix %s \n', prefix{idx});
    for suffix = bound(1) : bound(2)
        sceneName = strcat('img', prefix{idx}, num2str(suffix));
        sceneData = rdt.readArtifact(sceneName);
        
        allData.image{count} = sceneData.scene.data.photons;
        allData.wave{count}  = sceneData.scene.spectrum.wave;
        count = count + 1;
        
        fprintf('image count %d \n', count - 1);
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

