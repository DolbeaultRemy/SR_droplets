


function [od_imgs,Date] = run_data_remy(folderPath, run, cam, groupList, span, center, fraction ,shotwindow, removeFringes)
    

    %% input example:
    % folderPath   = "\\DyLabNAS\Data\TwoDGas\2025\12\02\\";
    % run          = "0000";
    % groupList    = ["/images/MOT_3D_Camera/in_situ_absorption", "/images/ODT_1_Axis_Camera/in_situ_absorption", "/images/ODT_2_Axis_Camera/in_situ_absorption", "/images/Horizontal_Axis_Camera/in_situ_absorption", "/images/Vertical_Axis_Camera/in_situ_absorption"];
    % cam          = 5; %% (vertical)
    % span         = [150, 150]; %% ROI stuff (I larger but if zou are mainly doing in-situ... )
    % center       = [1420,2050]; %% ROI stuff (For recent times should be around here)
    % fraction     = [0.1, 0.1]; %% Chooses fraction in corner
    % shotwindow   = 1 %% amountof shots in each run
    % removeFringes = false;

    folderPath     = strcat(folderPath, run);    
    filePattern = fullfile(folderPath, '*.h5');
    files       = dir(filePattern);
    k = 1;
    textprogressbar('Reading Data: ');
    Date        = files(1).date;
    for j = 1:shotwindow
        baseFileName = files(j).name;
        fullFileName = fullfile(files(j).folder, baseFileName);

        % fprintf(1, 'Now reading %s\n', fullFileName);
        
        atm_img  = double(imrotate(h5read(fullFileName, append(groupList(cam), "/atoms")), 0)); % im2double rescales values to between [0, 1], use double instead
        bkg_img  = double(imrotate(h5read(fullFileName, append(groupList(cam), "/background")), 0));
        dark_img = double(imrotate(h5read(fullFileName, append(groupList(cam), "/dark")), 0));

        refimages(:,:,k)  = subtractBackgroundOffset(cropODImage(bkg_img, center, span), fraction)';
        absimages(:,:,k)  = subtractBackgroundOffset(cropODImage(calculateODImage(atm_img, bkg_img, dark_img), center, span), fraction)';
        
        textprogressbar(k/shotwindow*100)
        k = k + 1;
    end
    textprogressbar(' Done');
    if removeFringes
        fprintf(1, '--Removing Fringes--\n');
        optrefimages               = removefringesInImage(absimages, refimages);
        absimages_fringe_removed   = absimages(:, :, :) - optrefimages(:, :, :);

        nimgs                      = size(absimages_fringe_removed,3);
        od_imgs                    = cell(1, nimgs);
        for i = 1:nimgs
            od_imgs{i}             = absimages_fringe_removed(:, :, i);
        end
    else
        nimgs                      = shotwindow;
        od_imgs                    = cell(1, nimgs);
        for i = 1:nimgs
            od_imgs{i}             = smoothdata(absimages(:, :, i));
        end
    end
end