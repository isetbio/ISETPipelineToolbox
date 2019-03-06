classdef CrossValidation < handle
%CROSSVALIDATION Generic tool for evaluation of an estimator.
%
% Syntax: cv = CrossValidation(testInput, testOutput, nTest)
%
% Description:
%   This class implements some generic function for evaluation and cross
%   validation of estimator.
%
% CrossValidation Methods:
%   eval          - Compute the reconstruction and compute MSE for a given
%                   estimator, and a pair of input - groundtruth.
%   sampleTest    - Run obj.eval on a random sample from test set.
%   evalTest      - Run obj.eval on the entire test set, compute total MSE
%                   and the distribution (i.e. a list) of MSE.
%   crossValidate - Run obj.evalTest with different hyperparameters
%                   registered in Estimator.regParaList
%
% Inputs:
%   testConeVec   - Input test set.
%   testImage     - Output test set.
%   nTest         - Number of test samples we are interested in running
%
% Outputs:
%   CrossValidation object.
    
    properties
        testImage;
        testConeVec;        
        nTest;

        imageDim = [32, 32, 3];        
    end
    
    methods
        
        function obj = CrossValidation(testConeVec, testImage, nTest)
            obj.testConeVec = testConeVec;
            obj.testImage   = testImage;            
            obj.nTest = nTest;
        end
        
        function [recon, mse] = eval(obj, estimator, input, image, visualization)
            recon = estimator.estimate(input);
            mse = sum((recon - image) .^ 2) / prod(obj.imageDim);
            
            if (visualization)
            subplot(1, 2, 1);
            imshow(reshape(image, obj.imageDim), 'InitialMagnification', 300);
            
            subplot(1, 2, 2);
            imshow(reshape(recon, obj.imageDim), 'InitialMagnification', 300);
            end            
            
        end
        
        function [recon, testImg, mse] = sampleTest(obj, estimator, visualization)
            testIdx = randi(obj.nTest, 1, 1);
            coneVec = obj.testConeVec(testIdx, :);
            testImg = obj.testImage(testIdx, :);
            
            [recon, mse] = obj.eval(estimator, coneVec, testImg, visualization);
        end
        
        function [totalMSE, listMSE] = evalTest(obj, estimator)
            mses = zeros(1, obj.nTest);
            parfor idx = 1:obj.nTest
                try
                    [~, mses(idx)] = obj.eval(estimator, obj.testConeVec(idx, :), obj.testImage(idx, :), false);                
                catch E
                    fprintf('Error in Reconstruction: %s for input %d \n', E.identifier, idx);
                    mses(idx) = NaN;
                end                
            end
            listMSE  = mses;
            totalMSE = sum(mses);          
        end
        
        function [paraList, mse] = crossValidate(obj, estimator)
            paraList = estimator.cvIterList();
            nPara = length(paraList);
            
            mse = zeros(1, nPara);
            for idx = 1:nPara
                fprintf('Cross Validate, Iteration %d / %d \n', idx, nPara);
                estimator.setRegPara(paraList(idx));
                mse(idx) = obj.evalTest(estimator);
            end
            
            figure;
            plot(paraList, mse, '-o', 'LineWidth', 2); grid on;
            title('Cross Validation');
            xlabel('Hyperparameter Value');
            ylabel('Total Mean Squared Error')            
        end
        
    end
    
end

