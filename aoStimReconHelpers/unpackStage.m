function st = unpackStage(pr, cnv, stage)
% Description:
%    Unpacks the pr and cnv struct parameters specific to either the forward
%    or recon conditions. Bypasses the need to send in a bunch of variables
%    at once, with the drawback that don't immediately know what gets used.
%    Only applied to the builders, in other files manually distinguish
%    forward from recon.
%
% See also: buildRenderStruct, buildMosaicMontage, aoStimReconRunMany
%
% History:
%   03/15/24  chr  branched off into callable function
%% Establish the struct and apply the stage
st = struct;

switch stage
    case "forward"
        st.pupilDiamMM = cnv.forwardPupilDiamMM;
        st.aoRender = pr.forwardAORender;
        st.noLCA = pr.forwardNoLCA;
        st.defocusDiopters = pr.forwardDefocusDiopters;
        st.randSeed = pr.forwardRandSeed;
        st.eccVars = pr.forwardEccVars;
        st.subjectID = pr.forwardSubjectID;
        st.zernikeDataBase = pr.forwardZernikeDataBase;
        st.renderDirFull = cnv.forwardRenderDirFull;
%         st.montageDirFull = cnv.forwardMontageDirFull;
    case "recon"
        st.pupilDiamMM = cnv.reconPupilDiamMM;
        st.aoRender = pr.reconAORender;
        st.noLCA = pr.reconNoLCA;
        st.defocusDiopters = pr.reconDefocusDiopters;
        st.randSeed = pr.reconRandSeed;
        st.eccVars = pr.reconEccVars;
        st.subjectID = pr.reconSubjectID;
        st.zernikeDataBase = pr.reconZernikeDataBase;
        st.renderDirFull = cnv.reconRenderDirFull;
%         st.montageDirFull = cnv.reconMontageDirFull;
    otherwise
        error('Unrecognized stage')
end
end