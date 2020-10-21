function iE = findLatencyEvent(EEG, targetEvent)
% findLatencyStory Return latency (in samples) corresponding to the (first)
% event of type 'targetEvent'
%
eventType = {EEG.event(:).type};

iE = find(strcmp(eventType, targetEvent), 1);
if ~isempty(iE)
    iE = EEG.event(iE).latency;
end

end
