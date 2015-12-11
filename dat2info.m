function allExpInfos = dat2info(filename)

if nargin<1
    [filename, folder] = uigetfile('G:\');
    filename = fullfile(folder, filename);
end

fprintf('Loading data from %s ...\n', filename)
data = load(filename); % loading the dat structure
data = data.dat; % otherwise structure dat is conflicting with the +dat package

% figureing out which experiments are included in these data

animal = data.ops.mouse_name;
expDate = data.ops.date;
experiments = data.ops.expts;
nExps = length(experiments);

ExpRefs = cell(nExps, 1);
expExists = nan(nExps, 1);
for iExp = 1:nExps
    ExpRefs{iExp} = dat.constructExpRef(animal, expDate, experiments(iExp));
    expExists(iExp) = dat.expExists(ExpRefs{iExp});
end

fprintf('In this file we have data from the following experiments:\n');
disp(ExpRefs);

if sum(~expExists)
    
    fprintf('Oh, wow, it looks like not all of the experiments above actually exist...\n');
    fprintf('Specifically, the missing experiments are:\n');
    disp(ExpRefs(~expExists));
    fprintf('I will now try to figure out if there was a confusion with the animal name...\n');
    suspects = {};
    list = dat.listSubjects;
    for iSubject = 1:length(list)
        if ~isempty(strfind(list{iSubject}, animal)) || ~isempty(strfind(animal, list{iSubject}))
            if ~strcmp(list{iSubject}, animal)
                suspects{end+1} = list{iSubject};
            end
        end
    end
    
    missingExps = experiments(find(~expExists));
    ref = cell(length(suspects), length(missingExps));
    isFound = nan(size(ref));
    for iExp = 1:length(missingExps)
        for iSuspect = 1:length(suspects)
            ref{iSuspect, iExp} = dat.constructExpRef(suspects{iSuspect}, expDate, missingExps(iExp));
            isFound(iSuspect, iExp) = dat.expExists(ref{iSuspect, iExp});
%         iExp
        end
    end
    
    iSuspect = find(all(isFound, 2));
    if isempty(iSuspect)
        fprintf('Hmmm, I wasn''t able to find these missing experiments, aborting ...\n');
        return;
    elseif length(iSuspect)>1
        fprintf('Hmmm, there is too much ambiguity with experiment names,\n you will need to look closer into it, aborting ...\n');
        return;
    end
    altName = suspects{iSuspect};
    for iExp = 1:length(missingExps)
            ExpRefs{find(experiments == missingExps(iExp))} = dat.constructExpRef(altName, expDate, missingExps(iExp));
    end
    
    fprintf('OK, my best guess is that these experiments are actually:\n');
    disp(ExpRefs(~expExists));
    fprintf('If that is wrong abort the script and check the data\n');
    fprintf('Will extract data from the following experiments now:\n');
    disp(ExpRefs);
    
end

%% Now, after figuring out the experiments, let's to the job of building info structures

fprintf('Extracting data..');
allExpInfos = cell(nExps, 1);
startFrame = cumsum([1; data.ops.Nframes(:)]);
endFrame = cumsum(data.ops.Nframes(:));

for iExp = 1:nExps
    roiIdx = find(data.cl.iscell)';
    frameIdx = [startFrame(iExp):endFrame(iExp)]';
    evalc('info = ppbox.infoPopulate(ExpRefs{iExp});');
    iCh = 1;
    info.chData(iCh).baseName = sprintf('%s_plane%03.0f', info.basename2p, data.ops.iplane);
    info.chData(iCh).tiffFrames = info.nChannels * (data.ops.iplane - 1) + iCh + ...
        ([1:data.ops.Nframes(iExp)]'-1)*info.nPlanes*info.nChannels;
    info.chData(iCh).meanIntensity = nan(data.ops.Nframes(iExp), 1);
    info.chData(iCh).targetFrame = data.ops.mimg1;
    info.chData(iCh).registeredAVG = data.ops.mimg1(data.ops.yrange(:), data.ops.xrange(:));
    info.iPlane = data.ops.iplane;
    tiffName = fullfile(info.folder2p, [info.basename2p, '_001_001.tif']);
    [~, info.planeHeaders] = img.loadFrames(tiffName, 1, 1, 1);
    info.meanIntensity = nan(data.ops.Nframes(iExp), 1);
    info.basenameRaw = sprintf('%s_plane%03.0f_raw', info.basename2p, data.ops.iplane);
    info.basenamePlane = sprintf('%s_plane%03.0f', info.basename2p, data.ops.iplane);
    info.planeFrames = data.ops.iplane + ([1:data.ops.Nframes(iExp)]'-1) * info.nPlanes;
    info.registrationChannel = data.ops.gchannel;
    info.targetFrame = data.ops.mimg1;
    [nY, nX] = size(info.targetFrame);
    info.targetFrameCrop = {info.targetFrame};
    info.cropPosSource.ymin = 1;
    info.cropPosSource.ymax = nY;
    info.cropPosSource.xmin = 1;
    info.cropPosSource.xmax = nX;
    info.cropPosDest.ymin = 1;
    info.cropPosDest.ymax = nY;
    info.cropPosDest.xmin = 1;
    info.cropPosDest.xmax = nX;
    info.nCrops = 1;
    info.collage = info.targetFrame;
    info.validX = data.ops.xrange(:);
    info.validY = data.ops.yrange(:);
    info.dx = data.ops.DS(frameIdx, 2);
    info.dy = data.ops.DS(frameIdx, 1);
    info.headers = info.planeHeaders;
    nROIs = length(roiIdx);
    info.ROI.CellMaps = {data.stat(roiIdx).ipix};
    info.ROI.CellClasses = repmat({'s'}, 1, nROIs);
    rectSource = [data.ops.yrange(1), data.ops.xrange(1), data.ops.yrange(end), data.ops.xrange(end)];
    rectDest = [1, 1, size(data.ops.mimg1)]; % full size of the targetFrame
    info.targetFrameROI = ppbox.remapROI(info.ROI.CellMaps, rectSource, rectDest);
    info = ppbox.getROIExtras(info);
    info.F = double(data.F.Fcell{iExp}(roiIdx, :)');
    info.basenameCrop = sprintf('%s_plane%03.0f', info.basename2p, data.ops.iplane);
    info.registeredMIP = data.ops.mimg1(data.ops.yrange(:), data.ops.xrange(:));
    info.collageMIP = data.ops.mimg1;
    info.mergeInfo.iCrop = ones(1, nROIs);
    info.mergeInfo.iCell = 1:nROIs;
    allExpInfos{iExp} = info;
end
fprintf('.done\n');
