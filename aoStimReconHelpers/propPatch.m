function mosaicConeInfo = propPatch(mosaicConeInfo) 
% Temporary function to correct for previous error in how mosaic cone
% proportions was being calculated

if (~isfield(mosaicConeInfo,'version') | mosaicConeInfo.version < 1)
    for i = 1:length(mosaicConeInfo.targetPropsL)
        propL(i) = mosaicConeInfo.targetPropsL(i);
        propS(i) = mosaicConeInfo.targetPropsS(i);

        mosaicConeInfo.targetPropsL(i) = propL(i);
        mosaicConeInfo.targetPropsM(i) = 1 - propL(i);
        mosaicConeInfo.targetPropsS(i) = propS(i);
        
        mosaicConeInfo.targetPropsLOfWholeMosaic(i) = propL(i)*(1-propS(i));
        mosaicConeInfo.targetPropsMOfWholeMosaic(i) = (1 - propL(i))*(1-propS(i));
        mosaicConeInfo.targetPropsSOfWholeMosaic(i) = propS(i);
        
        mosaicConeInfo.achievedPropsL(i) = mosaicConeInfo.numL(i) / (mosaicConeInfo.numL(i) + mosaicConeInfo.numM(i));
        mosaicConeInfo.achievedPropsM(i) = mosaicConeInfo.numM(i) / (mosaicConeInfo.numL(i) + mosaicConeInfo.numM(i));
        mosaicConeInfo.achievedPropsS(i) = mosaicConeInfo.numS(i) / mosaicConeInfo.numTotal(i);
        
        mosaicConeInfo.achievedPropsLOfWholeMosaic(i) = mosaicConeInfo.numL(i) / mosaicConeInfo.numTotal(i);
        mosaicConeInfo.achievedPropsMOfWholeMosaic(i) = mosaicConeInfo.numM(i) / mosaicConeInfo.numTotal(i);
        mosaicConeInfo.achievedPropsSOfWholeMosaic(i) = mosaicConeInfo.numS(i) / mosaicConeInfo.numTotal(i);
    end

    clear propL propS
end
end