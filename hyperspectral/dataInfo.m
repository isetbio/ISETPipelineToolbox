%% Dataset Browser
rdt = RdtClient('isetbio');
rdt.crp('/resources/scenes/hyperspectral');
rdt.openBrowser

%% Manchster DB
hyperspectralRootDir = '/resources/scenes/hyperspectral';
rdt = RdtClient('isetbio');
rdt.crp(hyperspectralRootDir);

% Select a database
dataBaseName = 'manchester_database';

% Select a scene family - Note: some data bases do not have this
dataBaseSceneFamily = '2002';
rdt.crp(fullfile(hyperspectralRootDir, dataBaseName, dataBaseSceneFamily));

% Select a scene name whithin the family
% sceneName = 'FreshFruit_Cx';
sceneName = 'scene1_2002_iset';
% Fetch the hyperspectral data
hyperspectralData = rdt.readArtifact(sceneName);

% Generate an ISETbio from the hyperspectralData
theScene = sceneFromBasis(hyperspectralData);
% Plot it
figure();
image(sceneGet(theScene, 'rgbimage'));

%% Stanford DB
hyperspectralRootDir = '/resources/scenes/hyperspectral';
rdt = RdtClient('isetbio');
rdt.crp(hyperspectralRootDir);

% Select a database
dataBaseName = 'stanford_database';

% Select a scene family - Note: some data bases do not have this
dataBaseSceneFamily = 'landscape';
rdt.crp(fullfile(hyperspectralRootDir, dataBaseName, dataBaseSceneFamily));

% Select a scene name whithin the family
sceneName = 'StanfordMemorial';
% Fetch the hyperspectral data
hyperspectralData = rdt.readArtifact(sceneName);

% Generate an ISETbio from the hyperspectralData
theScene = sceneFromBasis(hyperspectralData);
% Plot it
figure();
image(sceneGet(theScene, 'rgbimage'));

%% Penn DB 
hyperspectralRootDir = '/resources/scenes/hyperspectral';
rdt = RdtClient('isetbio');
rdt.crp(hyperspectralRootDir);

% Select a database
dataBaseName = 'penn_database';

% Select a scene family - Note: some data bases do not have this
rdt.crp(fullfile(hyperspectralRootDir, dataBaseName));

% Select a scene name whithin the family
sceneName = 'BearFruitGrayY';
% Fetch the hyperspectral data
scene = rdt.readArtifact(sceneName);

% Plot it
figure();
image(sceneGet(scene.scene, 'rgbimage'));

%% Harvard DB
hyperspectralRootDir = '/resources/scenes/hyperspectral';
rdt = RdtClient('isetbio');
rdt.crp(hyperspectralRootDir);

% Select a database
dataBaseName = 'harvard_database';

% Select a scene family - Note: some data bases do not have this
rdt.crp(fullfile(hyperspectralRootDir, dataBaseName));

% Select a scene name whithin the family
sceneName = 'imgh0';

% Fetch the hyperspectral data
scene = rdt.readArtifact(sceneName);

% Plot it
figure();
image(sceneGet(scene.scene, 'rgbimage'));
