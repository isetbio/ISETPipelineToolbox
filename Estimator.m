classdef (Abstract) Estimator < handle
    %ESTIMATOR Estimate image from cone reseponse.
    
    properties        
        regParaList;
        lastRecon;        
    end
    
    methods
                
        function paraList = cvIterList(obj)
            paraList = obj.regParaList;
        end
        
        function visualizeRecon(obj)
            figure;
            reconImage = reshape(obj.lastRecon, [32, 32, 3]);
            imshow(reconImage, 'InitialMagnification', 300);
        end
    end
    
    methods (Abstract)
        
        setRegPara(obj, para)
        
        reconImage = estimate(obj, input)                    
            
    end
    
end

