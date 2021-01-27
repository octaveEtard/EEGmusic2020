%
% figure_2
% Part of the EEGmusic2020 code (github.com/octaveEtard/EEGmusic2020)
% Author: Octave Etard, 2020
%
% Plot Figure 2.
%
% 1 - Load MIDI files and compute pitch distribution.
% 2 - Load EEG files and compute auto-correlation.
% 3 - Plot.
%
% This code requires the miditoolbox and EEGLAB toolboxes.
%
%% get data to plot track pitch distribution
% pathToMIDIfiles = EEGmusic2020.getPath('MIDI');
pathToMIDIfiles = '../../data/stimuli/inventionsMIDI';

tracks = [2:4,7:9]; % invention ID that were played
nTracks = numel(tracks);

pitch_parts = cell(nTracks,2);

bpm = 60;
% 3rd column in matrix indicates the track which corresponds to instrument
iChan = 3;
% 4th column is MIDI note
iPitch = 4;

baseF_piano = 440;
baseF_guitar = 220;

% this the order in the MIDI files, do not change (see below)
instrumentOrder = {'Guitar','Piano'};
nInstru = numel(instrumentOrder);

for iiTrack = 1:nTracks
    iTrack = tracks(iiTrack);
    
    [nmat, mstr] = readmidi( fullfile(pathToMIDIfiles,...
        sprintf('invent%i',iTrack),...
        sprintf('invent%i_%ibpm.mid',iTrack,bpm)) );
    
    % this column 
    nmatGuitar = nmat(~nmat(:,iChan),:);
    nmatPiano = nmat(logical(nmat(:,iChan)),:);
    
    % 1 is guitar, 2 is piano
    pitch_parts{iiTrack,1} = EEGmusic2020.midiNoteToPitch(nmatGuitar(:,iPitch),baseF_guitar);
    pitch_parts{iiTrack,2} = EEGmusic2020.midiNoteToPitch(nmatPiano(:,iPitch),baseF_piano);
end

% compute pitch distribution
pitch = [];
pitch{1} = vertcat(pitch_parts{:,1});
pitch{2} = vertcat(pitch_parts{:,2});

%
pitch_pdf = cell(nInstru,1); % actually a PMF
pitch_cdf = cell(nInstru,1);
pitch_uval = cell(nInstru,1);

for iInstru = 1:nInstru
    
    % pitch takes on a set of discrete values
    pitch_uval{iInstru} = unique( pitch{iInstru} );
    be = [pitch_uval{iInstru};pitch_uval{iInstru}(end)+1]; % bin Edges
    
    pitch_pdf{iInstru} = histcounts(pitch{iInstru},be,'Normalization','probability'); % PMF
    pitch_cdf{iInstru} = histcounts(pitch{iInstru},be,'Normalization','cdf');
end


%% Compute auto-correlations
Fs = 5000; % sampling rate
% xcorr parameters
maxLagT = 1;
maxLag = ceil(maxLagT * Fs);
% recording indices
allParts = 2:7;
nParts = numel(allParts);


%% Compute feature (music) auto-correlation
featureOpt = [];
featureOpt.Fs = Fs;
% where are features stored
featureOpt.baseFolder = EEGmusic2020.getPath('features');
% processing applied to the features
featureOpt.proc = 'LP-2000';
featureOpt.typeName = 'waveform';
% field 'attended' will correspond to the 1st instrument and field
% 'ignored' to the other one
featureOpt.condition = ['f',instrumentOrder{1}(1),'c'];
fields = {'attended','ignored'};

p = EEGmusic2020.makePathFeatureFolder(featureOpt);


axc_music = nan(2*maxLag+1,nParts,nInstru);

for iPart = 1:nParts
    featureOpt.part = allParts(iPart);
    f = EEGmusic2020.makeNameFeatureFile(featureOpt);

    d = load(fullfile(p,f));

    for iInstru = 1:nInstru
        axc_music(:,iPart,iInstru) = xcorr(d.(fields{iInstru}),maxLag,'normalized');
    end
end

axc_music = squeeze(mean(axc_music,3));


%% Compute EEG auto-correlation
% subject IDs to use
allSID = 1:17;
% 'EBIP01' to 'EBIP17'
allSID = arrayfun(@(idx) sprintf('EBIP%02i',idx),allSID,'UniformOutput',false);

allConditions = {'fGs','fPs'};

EEGopt = [];
% processing applied to the EEG
EEGopt.proc = 'HP-130';
% where is the EEG data stored
EEGopt.baseFolder = EEGmusic2020.getPath('EEG','processed');
% file type to load
EEGopt.ext = '.set';
EEGopt.Fs = Fs;

nSub = length(allSID);
nCond = length(allConditions);

axc_EEG = nan(2*maxLag+1,2,nParts,nSub,nCond);

for iCond = 1:nCond
    EEGopt.condition = allConditions{iCond};
    
    for iSub = 1:nSub
        EEGopt.SID = allSID{iSub};
        
        % EEG data folder
        p = EEGmusic2020.makePathEEGFolder(EEGopt);
        
        for iPart = 1:nParts
            EEGopt.part = allParts(iPart);
            % EEG file
            f = EEGmusic2020.makeNameEEGDataFile(EEGopt);
            
            [EEG,iB,iE] = EEGmusic2020.loadEEG({{p,f}});
            % data during stiumlus
            EEG = EEG(iB:iE,:);
            
            for iChan = 1:2
                axc_EEG(:,iChan,iPart,iSub,iCond) = xcorr(EEG(:,iChan),maxLag,'normalized');
            end
        end
    end
end

axc_EEG = squeeze(mean(axc_EEG,2:5));

%% Plotting
[col,col_null] = EEGmusic2020.plotStyleArgs(); % colors

fts = 11;
lwd = 1;
mks = 3;

% 3 row plot
xSpan = [1,1,2];
ySpan = [5,5,4];
% 
% 3 row plot
nTiles_x = xSpan(end);
nTiles_y = ySpan(1)+ySpan(3);

fig = figure;
tl = tiledlayout(nTiles_y,nTiles_x);

axx = gobjects(1,3);

axx(1) = nexttile(tl,[ySpan(1),xSpan(1)]); hold on;

for iInstru = 1:2
    plot(pitch_uval{iInstru},pitch_pdf{iInstru},'-o',...
        'Color',col(iInstru,:),'LineWidth',lwd,...
        'MarkerFaceColor','w','MarkerSize',mks);
end


axx(2) = nexttile(tl,[ySpan(2),xSpan(2)]); hold on;

for iInstru = 1:2
    plot(pitch_uval{iInstru},pitch_cdf{iInstru},'-o',...
        'Color',col(iInstru,:),'LineWidth',lwd,...
        'MarkerFaceColor','w','MarkerSize',mks);
end

% auto-correlation plots
tms = 1e3*(-maxLag:maxLag)/Fs;

% EEG
axx(end) = nexttile(tl,[ySpan(end),xSpan(end)]); hold on;

for iInstru = 1:2
    plot(tms,axc_EEG,'Color',col_null(1,:),'LineWidth',lwd);
end

% formatting
axx(1).YAxis.Label.String = 'Probability';

m = min(vertcat(pitch_uval{:}));
M = max(vertcat(pitch_uval{:}));

for ax = axx(1:2)
    ax.XAxis.Label.String = 'Pitch (Hz)';
    ax.XAxis.Limits = [m,M];
    ax.XAxis.TickValues = [100,300,500];
end

lg = legend(axx(2), instrumentOrder,'Box','off','Location','southeast',...
    'FontSize',fts);

for ax = axx(3:end)
    ax.XAxis.Limits = 50*[-1,1];
    ax.YAxis.Limits = [-0.6,1];
    ax.XAxis.TickValues = -45:15:45;
    ax.YAxis.Label.String = 'Correlation';
end
axx(end).XAxis.Label.String = 'Time (ms)';

pltools.formatAxisLabels(axx,fts,lwd);
tl.Padding = 'none';
tl.TileSpacing = 'none';

pltools.topLeftLabel(axx(1),'A',fts);
pltools.topLeftLabel(axx(2),'B',fts);
pltools.topLeftLabel(axx(3),'C',fts);

width = 9;
height = 7;
fileName = 'figure_2';
% pltools.printFigure(fig,'',fileName,600,width,height,1,1,1,'pdf',0);
%
%