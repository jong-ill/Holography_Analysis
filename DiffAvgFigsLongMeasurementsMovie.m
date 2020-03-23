function [ PreStimAvg, PostStimAvg, AvgMovie,DffAvgMovie] = DiffAvgFigsLongMeasurementsMovie(Stack ,frameidxtemp ,prefigs ,postfigs,Omitpost )
%UNTITLED3 Summary of this function goes here
%   This function takes the fluorescense signal and average it over all
%   trials before and after the stimulation onset



PreStim = zeros(size(Stack,1),size(Stack,2),prefigs,size(frameidxtemp,2));
PostStim = zeros(size(Stack,1),size(Stack,2),postfigs,size(frameidxtemp,2));
for stim=1:1:size(frameidxtemp,2) % Changed to 2, in case not enough frames before first stim
    PreStim(:,:,:,stim) = Stack(:,:,frameidxtemp(stim)-prefigs:frameidxtemp(stim)-1); % omit 1 frames before stim to reject artefacts
    PostStim(:,:,:,stim) = Stack(:,:,frameidxtemp(stim)+Omitpost:frameidxtemp(stim)+(Omitpost+postfigs-1)); % omit frames after stim to reject artefacts
end

PreStimAvg = mean(PreStim,4);
PostStimAvg = mean(PostStim,4);
AvgMovie = zeros(size(Stack,1),size(Stack,2),size(PreStimAvg,3)+size(PostStimAvg,3));
AvgMovie(:,:,1:size(PreStimAvg,3)) = PreStimAvg;
AvgMovie(:,:,size(PreStimAvg,3)+1:size(PreStimAvg,3)+size(PostStimAvg,3)) = PostStimAvg;

DffAvgMovie = (AvgMovie-mean(PreStimAvg,3))./mean(PreStimAvg,3);
% AvgMovie = [PreStimAvg PostStimAvg];
% DffStimAvg = (
% Diff=PostStim-PreStim;
% normDiff=Diff./PreStim;
% normDiff(isinf(normDiff))=mean(mean(mean(normDiff(~isinf(normDiff)))));

end
% dFvec=dFvec';
    