% dat2info wrapper script (including saving files)

[filename, folder] = uigetfile('G:\', '', '', 'multiselect', 'on');
if ~iscell(filename)
    filename = {filename};
end

nFiles = length(filename);

for iFile = 1:nFiles
    allInfos = dat2info(fullfile(folder, filename{iFile}));
    nInfos = length(allInfos);
    for iInfo = 1:nInfos
        meta = allInfos{iInfo};
        targetFolder = meta.folderProcessed;
        targetFile = fullfile(meta.folderProcessed, [meta.basenamePlane, '_ROI.mat']);
        if ~exist(targetFolder, 'dir')
            mkdir(targetFolder);
        end
        save(targetFile, 'meta');
    end
end