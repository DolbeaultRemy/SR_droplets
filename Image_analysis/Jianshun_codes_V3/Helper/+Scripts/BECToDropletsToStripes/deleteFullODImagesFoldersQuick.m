function deleteFullODImagesFoldersQuick(targetDir)
    % 快速删除（无确认提示）
    % 谨慎使用！
    
    if nargin < 1
        targetDir = pwd;
    end
    
    if ~exist(targetDir, 'dir')
        error('目标文件夹不存在');
    end
    
    % 查找文件夹
    dirInfo = dir(fullfile(targetDir, 'FullODImages*'));
    folderNames = {dirInfo([dirInfo.isdir]).name};
    folderNames = folderNames(~ismember(folderNames, {'.', '..'}));
    
    if isempty(folderNames)
        fprintf('未找到符合条件的文件夹。\n');
        return;
    end
    
    % 直接删除
    for i = 1:length(folderNames)
        folderPath = fullfile(targetDir, folderNames{i});
        try
            rmdir(folderPath, 's');
            fprintf('已删除: %s\n', folderPath);
        catch
            fprintf('删除失败: %s\n', folderPath);
        end
    end
end