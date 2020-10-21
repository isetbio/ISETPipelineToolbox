classdef (Abstract) Estimator < handle
    %ESTIMATOR Estimate image from cone reseponse.
    
    properties        
        regParaList;            
    end
    
    methods
                
        function paraList = cvIterList(obj)
            paraList = obj.regParaList;
        end
        
        function setParaList(obj, paraList)
            obj.regParaList = paraList;
        end
                
    end
    
    methods (Abstract)
        
        setRegPara(obj, para)
        
        reconImage = estimate(obj, input)                    
            
    end
    
end

