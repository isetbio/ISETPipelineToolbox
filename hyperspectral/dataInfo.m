%% Dataset Browser
rdt = RdtClient('isetbio');
rdt.crp('/resources/scenes/hyperspectral');
rdt.openBrowser

%% Read and load images
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

%% Read and load images another example
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
