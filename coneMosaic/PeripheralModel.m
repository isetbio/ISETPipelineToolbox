classdef PeripheralModel
    
    methods(Static)
        
        function [theConeMosaic, theOI] = eyeModelEcc(eccXrange, eccYrange, fovDegs, pupilDiam)
            % Select data for a particular subject
            % Use [] to indicate average subject
            subID = 6;
            desiredPupilDiamMM = pupilDiam;
            applyCentralCorrection = true;
            % Get a struct with the Polans data
            d = PeripheralModel.Polans2015Data(applyCentralCorrection);
            
            wavefrontSpatialSamples = 201;
            wavelengthsListToCompute = 450:100:750;
            
            for eccYindex = 1:numel(eccYrange)
                for eccXindex = 1:numel(eccXrange)
                    
                    % The eccentricity in degrees
                    eccXY = [eccXrange(eccXindex) eccYrange(eccYindex)];
                    
                    % Get cone spacing at this eccentricity
                    eccRadiusDegs = sqrt(sum(eccXY.^2,2));
                    eccAngleDegs = atan2d(eccXY(2), eccXY(1));
                    [coneSpacingInMeters, coneApertureInMeters] = coneSizeReadData('eccentricity', eccRadiusDegs, ...
                        'angle', eccAngleDegs, ...
                        'eccentricityUnits', 'deg', ...
                        'angleUnits', 'deg', ...
                        'whichEye', d.eye, ...
                        'useParfor', false);
                    
                    % Generate a 1.0 x 1.0 deg regular hex cone mosaic with ecc-adjusted cone separation and aperture
                    theConeMosaic = coneMosaicHex(13, ...
                        'fovDegs', fovDegs, ...
                        'customLambda', coneSpacingInMeters*1e6, ...
                        'customInnerSegmentDiameter', coneApertureInMeters*1e6);
                    
                    
                    % Get zCoeffs for this eccentricity
                    [zCoeffs, nearestEccXY] = PeripheralModel.zCoeffsForSubjectAndEccentricity(d, subID, eccXY);
                    
                    % Generate oi at this eccentricity
                    theOI = PeripheralModel.makeCustomOI(zCoeffs, d.measurementPupilDiameMM, d.measurementWavelength, ...
                        desiredPupilDiamMM, wavelengthsListToCompute, wavefrontSpatialSamples, nearestEccXY, d.eye);
                    
                end
            end
            
        end
        
        function [mosaic, psfObj, psfData] = eyeModelCmosaic(eccX, eccY, fovDegs, pupilDiam, subjectID)
            mosaicEcc = [eccX, eccY];
            
            % Generate mosaic centered at target eccentricity
            mosaic = cMosaic(...
                'sizeDegs', [1 1] * fovDegs, ...
                'eccentricityDegs', mosaicEcc, ...
                'noiseFlag', 'none', ...
                'integrationTime', 0.1, ...
                'useParfor', false);
            
            % Generate optics appropriate for the mosaic's eccentricity
            [oiEnsemble, psfEnsemble] = mosaic.oiEnsembleGenerate(mosaicEcc, ...
                'zernikeDataBase', 'Polans2015', ...
                'subjectID', subjectID, ...
                'pupilDiameterMM', pupilDiam, ...
                'subtractCentralRefraction', true);
            
            psfObj = oiEnsemble{1};
            psfData = psfEnsemble{1};
        end
        
        function theOI = makeCustomOI(zCoeffs, measPupilDiameterMM, measWavelength, ...
                desiredPupilDiamMM, wavelengthsListToCompute, wavefrontSpatialSamples, nearestEccXY, whichEye)
            
            showTranslation = false;
            
            [~, theOTF, xSfCyclesDeg, ySfCyclesDeg] = ...
                computePSFandOTF(zCoeffs, ...
                wavelengthsListToCompute, wavefrontSpatialSamples, ...
                measPupilDiameterMM, desiredPupilDiamMM, ...
                measWavelength, showTranslation);
            
            for waveIndex = 1:numel(wavelengthsListToCompute)
                theWaveOTF = squeeze(theOTF(:,:,waveIndex));
                theOTF(:,:,waveIndex) = ifftshift(theWaveOTF);
            end
            
            umPerDegree = 300;
            theOI = oiCreate('wvf human', desiredPupilDiamMM,[],[], umPerDegree);
            optics = oiGet(theOI,'optics');
            optics = opticsSet(optics, 'otfwave', wavelengthsListToCompute);
            
            % Update optics with new OTF data
            xSfCyclesPerMM = 1000*xSfCyclesDeg / umPerDegree;
            ySfCyclesPerMM = 1000*ySfCyclesDeg / umPerDegree;
            customOptics = opticsSet(optics,'otf data',theOTF);
            customOptics = opticsSet(customOptics, 'otffx',xSfCyclesPerMM);
            customOptics = opticsSet(customOptics,'otffy',ySfCyclesPerMM);
            
            % Update theOI with custom optics
            theOI = oiSet(theOI,'optics', customOptics);
            if (strcmp(whichEye, 'right'))
                if (nearestEccXY(1)<0)
                    xEccLabel = sprintf('x:%2.0f^o(N)', nearestEccXY(1));
                elseif(nearestEccXY(1)>0)
                    xEccLabel = sprintf('x:%2.0f^o(T)', nearestEccXY(1));
                elseif(nearestEccXY(1)==0)
                    xEccLabel = sprintf('x:%2.0f^o', nearestEccXY(1));
                end
            elseif (strcmp(whichEye, 'left'))
                if (nearestEccXY(1)<0)
                    xEccLabel = sprintf('x:%2.0f^o(T)', nearestEccXY(1));
                elseif(nearestEccXY(1)>0)
                    xEccLabel = sprintf('x:%2.0f^o(N)', nearestEccXY(1));
                elseif(nearestEccXY(1)==0)
                    xEccLabel = sprintf('x:%2.0f^o', nearestEccXY(1));
                end
            end
            yEccLabel = sprintf('y:%2.0f^o', nearestEccXY(2));
            theOI = oiSet(theOI,'name', sprintf('%s  %s', xEccLabel, yEccLabel));
        end
        
        function [theZcoeffs, nearestEccXY] = zCoeffsForSubjectAndEccentricity(d, subjectIndex, eccXY)
            dX = d.eccXgrid - eccXY(1);
            dY = d.eccYgrid - eccXY(2);
            [~,indexOfNearestEcc] = min(dX.^2+dY.^2);
            nearestEccXY = [d.eccXgrid(indexOfNearestEcc) d.eccYgrid(indexOfNearestEcc)];
            
            if (strcmp(d.source, 'Polans_2015'))
                if (isempty(subjectIndex))
                    % mean over all subjects
                    theMeasuredZcoeffs = mean(squeeze(d.zCoeffs(:,indexOfNearestEcc,:)),1);
                else
                    theMeasuredZcoeffs = squeeze(d.zCoeffs(subjectIndex,indexOfNearestEcc,:));
                end
            else
                eyeIndex = 1;
                if (isempty(subjectIndex))
                    theMeasuredZcoeffs = squeeze(d.zCoeffsMean(eyeIndex,indexOfNearestEcc,:));
                else
                    theMeasuredZcoeffs = squeeze(d.zCoeffs(eyeIndex,subjectIndex,indexOfNearestEcc,:));
                end
                
            end
            
            % Place zCoeffs in right bin
            theZcoeffs = zeros(1,21);
            theZcoeffs(d.zCoeffOSAIndices+1) = theMeasuredZcoeffs;
        end
        
        function d = Polans2015Data(applyCentralCorrection)
            % Load the raw data from Polans et al (2015), "Wide-field optical model
            %     of the human eye with asymmetrically tilted and decentered lens
            %     that reproduces measured ocular aberrations", Optica, 2(2), 2015,
            %     pp.124-134
            %
            
            % Load the matfile
            allData = rawDataReadData('zCoefsPolans2015', ...
                'datatype', 'isetbiomatfileonpath');
            allData = allData.data;
            
            d.source = 'Polans_2015';
            d.subjectsNum = size(allData,1);
            % Retrieve the Y-eccentricity grid
            subjectIndex = 1; entryIndex = 1;
            d.eccYgrid = allData(subjectIndex,:,entryIndex);
            
            
            % Retrieve the X-eccentricity grid
            entryIndex = 2;
            d.eccXgrid = allData(subjectIndex,:,entryIndex);
            
            % Retrieve the z-coeffs
            entryIndices = 3:20;
            d.zCoeffs = allData(:,:,entryIndices);
            if (any(isnan(d.zCoeffs(:))))
                fprintf(2,'Found z-coeff with NaN values');
            end
            
            
            % Add zcoeff OSA indices and names
            d.zCoeffOSAIndices = entryIndices;
            d.zCoeffNames = {...
                'Z3',  'oblique astigmatism'; ...
                'Z4',  'defocus'; ...
                'Z5',  'vertical astigmatism'; ...
                'Z6' , 'vertical trefoil'; ...
                'Z7',  'vertical coma'; ...
                'Z8',  'horizontal coma'; ...
                'Z9',  'oblique trefoil'; ...
                'Z10', 'oblique quadrafoil'; ...
                'Z11', 'oblique secondary astigmatism'; ...
                'Z12', 'primary spherical'; ...
                'Z13', 'vertical secondary astigmatism'; ...
                'Z14', 'vertical quadrafoil'; ...
                'Z15', '??'; ...
                'Z16', '??'; ...
                'Z17', '??'; ...
                'Z18', '??'; ...
                'Z19', '??'; ...
                'Z20', '??' ...
                };
            
            % Measurement params
            % The Polans 2015 measurements were all done in the right eye
            d.eye = 'right';
            
            % Using a 4 mm pupil
            d.measurementPupilDiameMM = 4;
            
            % And 550 nm light
            d.measurementWavelength = 550;
            
            if (applyCentralCorrection)
                d.zCoeffs = removeCentralRefractionError(d);
            end
            function zCoeffs = removeCentralRefractionError(d)
                zCoeffs = d.zCoeffs;
                % Find the zero ecc index
                zeroEccIndex = find(d.eccXgrid == 0 & d.eccYgrid == 0);
                % Find the z-coeff index corresponding to the defocus term
                theDefocusCoeffIndex = find(strcmp(d.zCoeffNames(:,2), 'defocus'));
                % Central corrections for all subjects
                allSubjectsCentralRefractionCorrection = zCoeffs(:, zeroEccIndex, theDefocusCoeffIndex);
                % Subtract central correction for all subjects and all eccentricities
                zCoeffs(:,:,theDefocusCoeffIndex) = bsxfun(@minus, ...
                    zCoeffs(:,:,theDefocusCoeffIndex), allSubjectsCentralRefractionCorrection);
            end
        end
        
    end
    
end