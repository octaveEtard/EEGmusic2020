function iSound = findIndexSoundChannel(EEG)
iSound = find(strcmp({EEG.chanlocs(:).labels},'Sound'),1);
end