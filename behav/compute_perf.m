%
% compute_perf
% Part of the JoNmusic2020 code (github.com/octaveEtard/JoNmusic2020)
% Author: Octave Etard, 2020
%
% Compute behavioural performance on the vibrato detection task (True/False
% Positive Rate), and save results.
%
% For each keyboard key press, the closest preceding vibrato is determined.
% Only inputs within at most maxRT seconds after a vibrato onset are
% considered to be valid.
%
% For each keyboard key press, if the closest preceding vibrato:
%   - was in the target track, it is a true positive;
%   - was in the ignored track, it is a false positive;
%   - no vibrato within the preceding maxRT seconds, it is counted
% as 'unprompted'.
%
%% Parameters
maxRT = 2; % in s
nVib = 10; % we expect 10 vibratos / track

% change this to the folder containing the behavioural data
baseFolder = JoNmusic2020.getPath('behav');
% folder containing time of the vibratos
vibTimeFolder = fullfile(baseFolder,'vibTime');

% will save the results here:
saveFolder = baseFolder;
saveName = sprintf('clickPerformance_RT_%.1f.mat', maxRT);

% subject IDs to use
SID = 2:17; % #1 missing ; % 'EBIP02' to 'EBIP17'

instrumentOrder = {'guitar','piano'};
% data from each condition divided in parts corresponding to different
% inventions: data parts indices to use here
% (1 was training block, not used)
parts = 2:7;

%%
% index to subject ID
SID = arrayfun(@(idx) sprintf('EBIP%02i',idx),SID,'UniformOutput',false);
nParts = numel(parts);
nSub = numel(SID);

% number of key presses... (nSub x target instrument x nParts)
nkPress = zeros(nSub,2,nParts);
% reaction time for ... (at most nVib true/false positive)
rtTP = nan(nVib,nSub,2,nParts); % ... true positive
rtFP = nan(nVib,nSub,2,nParts); % ... false positive
rtUmprompted = cell(nSub,2,nParts); % ... umprompted (time to closest preceding vibrato in this case)

for iiPart = 1:nParts
    
    % load true time of vibratos for the guitar & piano tracks
    iPart = parts(iiPart);
    vibTime = load(fullfile(vibTimeFolder,sprintf('vibTime_%i.mat',iPart)));
    
    % ensure vibrato times are in the order specified by instrumentOrder
    iOrder = 1:2;
    if ~all(strcmp(instrumentOrder,vibTime.instrumentOrder))
        iOrder = 2:1;
    end
    instrumentOrder_ = vibTime.instrumentOrder(iOrder);
    vibTime = vibTime.vibTiming(:,iOrder); % vibrato onset times
    assert(all(strcmp(instrumentOrder,instrumentOrder_)),'Invalid intrument order');
    
    for iSub = 1:nSub
        currentSID = SID{iSub};
        currentFolder = fullfile(baseFolder,currentSID);
        
        for iInstru = 1:2 % attended/target instrument
            condition = sprintf('f%sc',upper(instrumentOrder{iInstru}(1)));
            
            % vibTime(:,iInstru) now correspond to time of vibrato to
            % detect ; vibTime(:,3-iInstru) are the ones to ignore
            
            % key press times
            kPressTime = load(fullfile(currentFolder,sprintf('%s_keyboardInputs_%s_%i.mat',currentSID,condition,iPart)));
            kPressTime = kPressTime.timePressed;
            nInputs = numel(kPressTime);
            nkPress(iSub,iInstru,iiPart) = nInputs;
            
            if nInputs < 1
                % no key press to analyse
                continue;
            end
            % make sure these are sorted
            kPressTime = sort(kPressTime,'ascend');
            detectedVibratos = false(nVib,2);
            nUmp_ = 0;
            rtUmp_ = nan(nInputs,1);
            
            for ik = 1:nInputs
                % time between key press & each vibrato
                rt = kPressTime(ik) - vibTime;
                rt( rt < 0 ) = nan; % vibrato after key press
                % time to closest vibrato preceding key press in any track
                [rt,iVib] = min(rt,[],1:2,'omitnan','linear');
                % index of the corresponding vibrato, and track
                [iVib,iTrack] = ind2sub([nVib,2],iVib);
                
                % key press later than maxRT or before any vibrato
                % or if vibrato already detected --> unprompted
                if isnan(rt) || maxRT < rt || detectedVibratos(iVib,iTrack)
                    nUmp_ = nUmp_ + 1;
                    rtUmp_(nUmp_) = rt;
                else
                    detectedVibratos(iVib,iTrack) = true;
                    if iTrack == iInstru % true positive vibrato
                        rtTP(iVib,iSub,iInstru,iiPart) = rt;
                    else
                        rtFP(iVib,iSub,iInstru,iiPart) = rt;
                    end
                end
            end
            % Note that the rt are stored at indices indicating their
            % corresponding vibrato
            rtUmprompted{iSub,iInstru,iiPart} = rtUmp_(1:nUmp_);
        end
    end
end


%%
% number of TP/FP for each subject & part & each attended instrument
nTP = squeeze(sum(~isnan(rtTP),1));
nFP = squeeze(sum(~isnan(rtFP),1));
nUmprompted = cellfun(@numel,rtUmprompted);

% sanity check
assert(all(nTP + nFP + nUmprompted == nkPress,'all'));

% overall number of TP/FP across all parts for each subjects & instrument
nTP = sum(nTP,3);
nFP = sum(nFP,3);
% average RT
meanRT_TP = squeeze(sum(rtTP,[1,4],'omitnan')) ./ nTP;
meanRT_FP = squeeze(sum(rtFP,[1,4],'omitnan')) ./ nFP;
% TP/FP Rate
nVibTot = nParts * nVib;
TPR = nTP / nVibTot;
FPR = nFP / nVibTot;


%% Saving
d = struct();
d.instrumentOrder = instrumentOrder;
d.parts = parts;
d.SID = SID;

d.maxRT = maxRT;
d.nVib = maxRT;

d.nTP = nTP;
d.nFP = nFP;
d.nUmprompted = nUmprompted;

d.TPR = TPR;
d.FPR = FPR;

d.rtTP = rtTP;
d.rtFP = rtFP;
d.rtUmprompted = rtUmprompted;

d.meanRT_TP = meanRT_TP;
d.meanRT_FP = meanRT_FP;


LM.save(d,saveName,saveFolder);
%
%