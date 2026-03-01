folderPath = "Y:\StructuralPhaseTransition\2026\02\11\";
ROIcenter = [1500 1950];
measurementName = "Stripes_90°";
run = "0001";
cam = 5;
groupList    = ["/images/MOT_3D_Camera/in_situ_absorption", "/images/ODT_1_Axis_Camera/in_situ_absorption", "/images/ODT_2_Axis_Camera/in_situ_absorption", "/images/Horizontal_Axis_Camera/in_situ_absorption", "/images/Vertical_Axis_Camera/in_situ_absorption"];
span         = [400, 400];
center       = [1500, 1950];
fraction     = [0.1, 0.1];
shotwindow   = 380;
removeFringes = false;
detectionParams = setDetectionParameters();

[od_imgs,Date] = run_data_remy(folderPath, run, cam, groupList, span, center, fraction ,shotwindow, removeFringes);

for i = 1:length(od_imgs)
    [patchProps, patchCentroidsGlobal, imgCropped, xStart, yStart, binaryMask, CC] = ...
                detectStructure(od_imgs{i}, detectionParams, "Saved_masks/", "1");
    save("Saved_masks/"+string(i)+'.mat', "binaryMask", "xStart", "yStart");
end


