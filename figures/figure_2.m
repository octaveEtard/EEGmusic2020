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
pathToMIDIfiles = EEGmusic2020.getPath('MIDI');
tracks = [2:4,7:9]; % invention ID that were played
nTracks = numel(tracks);

pitch_parts = cell(nTracks,2);

bpm = 60;
% 3rd column in matrix is channel indicator
iChan = 3;
% 4th column is MIDI notes
iPitch = 4;

baseF_piano = 440;
baseF_guitar = 220;

for iiTrack = 1:nTracks
    iTrack = tracks(iiTrack);
    
    [nmat, mstr] = readmidi( fullfile(pathToMIDIfiles,...
        sprintf('invent%i',iTrack),...
        sprintf('invent%i_%ibpm.mid',iTrack,bpm)) );
    
    % this column indicate the track which corresponds to instrument
    iChan = 3;
    nmatGuitar = nmat(~nmat(:,iChan),:);
    nmatPiano = nmat(logical(nmat(:,iChan)),:);
    
    pitch_parts{iiTrack,1} = EEGmusic2020.midiNoteToPitch(nmatGuitar(:,iPitch),baseF_guitar);
    pitch_parts{iiTrack,2} = EEGmusic2020.midiNoteToPitch(nmatPiano(:,iPitch),baseF_piano);
    
end

% compute pitch distribution
pitch{1} = vertcat(pitch_parts{:,1});
pitch{2} = vertcat(pitch_parts{:,2});

%
pitch_pdf = cell(2,1); % actually a PMF
pitch_cdf = cell(2,1);
pitch_uval = cell(2,1);

for iInstru = 1:2
    
    % pitch takes on a set of discrete value
    pitch_uval{iInstru} = unique( pitch{iInstru} );
    be = [pitch_uval{iInstru};pitch_uval{iInstru}(end)+1]; % bin Edges
    
    pitch_pdf{iInstru} = histcounts(pitch{iInstru},be,'Normalization','probability'); % PMF
    pitch_cdf{iInstru} = histcounts(pitch{iInstru},be,'Normalization','cdf');
end


%% Compute EEG auto-correlation
% where is the EEG data stored
baseFolder = EEGmusic2020.getPath('EEG','processed');
% subject IDs to use
allSID = 1:17;
% 'EBIP01' to 'EBIP17'
allSID = arrayfun(@(idx) sprintf('EBIP%02i',idx),allSID,'UniformOutput',false);

allConditions = {'fGs','fPs'};
allParts = 2:7; % EEG parts

proc = {'HP-175','HP-115'};
Fs = 5000;

nSub = length(allSID);
nCond = length(allConditions);
nParts = length(allParts);

nEEGParts = nCond * nParts;

maxLagT = 1;
maxLag = ceil(maxLagT * Fs);

axc = nan(2*maxLag+1,2,nParts,nSub,nCond);

for iCond = 1:nCond
    for iSub = 1:nSub
        SID = allSID{iSub};
        condition = allConditions{iCond};
        
        % EEG data folder
        p = EEGmusic2020.makePathEEGFolder(baseFolder,proc{iCond},Fs,SID);

        for iPart = 1:nParts
            part = allParts(iPart);
            f = EEGmusic2020.makeNameEEGDataFile(proc{iCond},Fs,SID,condition,part,'.set');
            
            [EEG,iB,iE] = EEGmusic2020.loadEEG({{p,f}});
            % data during stiumlus
            EEG = EEG(iB:iE,:);
            
            for iChan = 1:2
                axc(:,iChan,iPart,iSub,iCond) = xcorr(EEG(:,iChan),maxLag,'normalized');
            end
        end
    end
end


%% Plotting
col = EEGmusic2020.plotStyleArgs(); % colors

fts = 11;
lwd = 1;
mks = 3;

% 2 row plot
xSpan = [1,1,2];
ySpan = [5,5,4];

% 2 row plot
nTiles_x = xSpan(end);
nTiles_y = ySpan(1)+ySpan(3);

fig = figure;
tl = tiledlayout(nTiles_y,nTiles_x);

ax1 = nexttile(tl,[ySpan(1),xSpan(1)]); hold on;

for iInstru = 1:2
    plot(pitch_uval{iInstru},pitch_pdf{iInstru},'-o',...
        'Color',col(iInstru,:),'LineWidth',lwd,...
        'MarkerFaceColor','w','MarkerSize',mks);
end


ax2 = nexttile(tl,[ySpan(2),xSpan(2)]); hold on;

for iInstru = 1:2
    plot(pitch_uval{iInstru},pitch_cdf{iInstru},'-o',...
        'Color',col(iInstru,:),'LineWidth',lwd,...
        'MarkerFaceColor','w','MarkerSize',mks);
end


ax3 = nexttile(tl,[ySpan(3),xSpan(3)]); hold on;

maxc = squeeze(mean(axc,2:4));
tms = 1e3*(-maxLag:maxLag)/Fs;

for iInstru = 1:2
    plot(ax3,tms,maxc(:,iInstru),'Color',col(iInstru,:),'LineWidth',lwd);
end

for iInstru = 1:2
    plot(tms,maxc(:,iInstru),'Color',col(iInstru,:),'LineWidth',lwd);
end

% formatting
ax1.YAxis.Label.String = 'Probability';

m = min(vertcat(pitch_uval{:}));
M = max(vertcat(pitch_uval{:}));

for ax = [ax1,ax2]
    ax.XAxis.Label.String = 'Pitch (Hz)';
    ax.XAxis.Limits = [m,M];
    ax.XAxis.TickValues = [100,300,500];
end

lg = legend(ax2, {'Guitar','Piano'},'Box','off','Location','southeast',...
    'FontSize',fts);

ax3.XAxis.Limits = 50*[-1,1];
ax3.YAxis.Limits = [-0.6,1];
ax3.XAxis.TickValues = -45:15:45;
ax3.YAxis.Label.String = 'Correlation';
ax3.XAxis.Label.String = 'Time (ms)';

pltools.formatAxisLabels([ax1,ax2,ax3],fts,lwd);
tl.Padding = 'none';
tl.TileSpacing = 'none';

pltools.topLeftLabel(ax1,'A',fts);
pltools.topLeftLabel(ax2,'B',fts);
pltools.topLeftLabel(ax3,'C',fts);
%
%