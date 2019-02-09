classdef CrossValidation < handle
    %CROSSVALIDATION Some useful tool for evaluation of an estimator.
    
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
        
        function [recon, mse] = sampleTest(obj, estimator)
            testIdx = randi(obj.nTest, 1, 1);
            coneVec = obj.testConeVec(testIdx, :);
            testImg = obj.testImage(testIdx, :);
            
            [recon, mse] = obj.eval(estimator, coneVec, testImg, true);
        end
        
        function [totalMSE, listMSE] = evalTest(obj, estimator)
            mses = zeros(1, obj.nTest);
            for idx = 1:obj.nTest
                [~, mses(idx)] = obj.eval(estimator, obj.testConeVec(idx, :), obj.testImage(idx, :), false);                
            end
            listMSE  = mses;
            totalMSE = sum(mses);            
        end
        
        function [paraList, mse] = crossValidate(obj, estimator)
            paraList = estimator.cvIterList();
            nPara = length(paraList);
            
            mse = zeros(1, nPara);
            for idx = 1:nPara
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

