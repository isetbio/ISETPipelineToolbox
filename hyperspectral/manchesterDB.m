function allData = manchesterDB(showPlot)

hyperspectralRootDir = '/resources/scenes/hyperspectral';
rdt = RdtClient('isetbio');
rdt.crp(hyperspectralRootDir);

% Select a database
dataBaseName = 'manchester_database';
nScene = 8;

if showPlot
    figure();
end

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
 
        if showPlot
            image(sceneGet(scene, 'rgbimage'));
        end
    end
end

end

