%% Read scene from database
[allData, nScene] = manchesterDB();

%% Show scenes
scene = sceneCreate('whitenoise');
scene.spectrum.wave = allData.wave;
for idx = 1 : nScene
    scene.data.photons = allData.image{idx};
    imshow(sceneGet(scene, 'rgbimage')); pause(0.5);
end

%% Read scene from database
[allData, nScene] = harvardDB();

%% Show scenes
scene = sceneCreate('whitenoise');
scene.spectrum.wave = allData.wave;
for idx = 1 : nScene
    scene.data.photons = allData.image{idx};
    size(scene.data.photons)
    imshow(sceneGet(scene, 'rgbimage')); pause(1.0);
end