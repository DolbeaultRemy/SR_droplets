% combine_plot_data.m
function combine_plot_data(data_folders, output_folder)
    % Merge and plot data from multiple datasets
    % Input:
    %   data_folders - cell array containing paths to data folders
    %   output_folder - output directory path for plots
    
    % Ensure output folder exists
    if ~exist(output_folder, 'dir')
        mkdir(output_folder);
    end
    
    % Initialize arrays to store all datasets
    all_datasets = struct();  % Changed to structure, stored by numerical value
    all_dataset_values = [];  % Store numerical labels of all datasets
    
    % Load each dataset
    for i = 1:length(data_folders)
        data_file = fullfile(data_folders{i}, 'plotting_data.mat');
        
        if exist(data_file, 'file')
            % Load data
            load(data_file);
            
            % Check if data label exists
            if ~isfield(plot_data, 'dataset_label')
                warning('Dataset %s does not have dataset_label field, using default label', data_folders{i});
                dataset_label_str = sprintf('Dataset_%d', i);
            else
                dataset_label_str = plot_data.dataset_label;
            end
            
            % Try to convert label to floating point number
            try
                dataset_value = str2double(dataset_label_str);
                
                % Check if conversion is successful
                if isnan(dataset_value)
                    error('Cannot convert to floating point number');
                end
            catch
                warning('Dataset label "%s" cannot be converted to floating point, using index %d', dataset_label_str, i);
                dataset_value = i;  % Use index as numerical value
            end
            
            % Generate valid field name (avoid decimal points)
            field_name = create_valid_fieldname(dataset_value);
            
            % Check if dataset with same numerical value already exists
            if ~isfield(all_datasets, field_name)
                % Create new dataset group
                all_datasets.(field_name) = struct();
                all_datasets.(field_name).value = dataset_value;
                all_datasets.(field_name).label = dataset_label_str;
                all_datasets.(field_name).indices = i;  % Store original index
                all_datasets.(field_name).folders = {data_folders{i}};
                
                % Initialize data storage
                all_datasets.(field_name).number_data = [];
                all_datasets.(field_name).area_data = [];
                all_datasets.(field_name).aspect_data = [];
                all_datasets.(field_name).spacing_data = [];
                
                % Add value to list
                all_dataset_values = [all_dataset_values; dataset_value];
            else
                % Add to existing dataset group
                all_datasets.(field_name).indices = [all_datasets.(field_name).indices, i];
                all_datasets.(field_name).folders = [all_datasets.(field_name).folders, data_folders{i}];
            end
            
            % Get current dataset group
            current_dataset = all_datasets.(field_name);
            
            % Extract data
            number_data = struct(...
                'name', sprintf('%.4f', dataset_value), ...
                'x', plot_data.detected_number.x, ...
                'mean', plot_data.detected_number.mean, ...
                'std', plot_data.detected_number.std);
            
            area_data = struct(...
                'name', sprintf('%.4f', dataset_value), ...
                'x', plot_data.detected_area.x, ...
                'mean', plot_data.detected_area.mean, ...
                'std', plot_data.detected_area.std);
            
            aspect_data = struct(...
                'name', sprintf('%.4f', dataset_value), ...
                'x', plot_data.aspect_ratio.x, ...
                'mean', plot_data.aspect_ratio.mean, ...
                'std', plot_data.aspect_ratio.std);
            
            spacing_data = struct(...
                'name', sprintf('%.4f', dataset_value), ...
                'x', plot_data.average_spacing.x, ...
                'mean', plot_data.average_spacing.mean, ...
                'std', plot_data.average_spacing.std);
            
            % Store data to current dataset group
            all_datasets.(field_name).number_data = [current_dataset.number_data; number_data];
            all_datasets.(field_name).area_data = [current_dataset.area_data; area_data];
            all_datasets.(field_name).aspect_data = [current_dataset.aspect_data; aspect_data];
            all_datasets.(field_name).spacing_data = [current_dataset.spacing_data; spacing_data];
            
            fprintf('Loaded dataset: %s (value: %.4f, field name: %s)\n', dataset_label_str, dataset_value, field_name);
        else
            warning('Data file not found: %s', data_file);
        end
    end
    
    % Sort by numerical value
    all_dataset_values = sort(unique(all_dataset_values));
    
    % Create merged data structure
    all_number_data = {};
    all_area_data = {};
    all_aspect_data = {};
    all_spacing_data = {};
    
    % For each unique value, create merged data
    for i = 1:length(all_dataset_values)
        dataset_value = all_dataset_values(i);
        field_name = create_valid_fieldname(dataset_value);
        
        if isfield(all_datasets, field_name)
            dataset_group = all_datasets.(field_name);
            
            % If multiple datasets with same value, calculate average
            if length(dataset_group.indices) > 1
                fprintf('Merging %d datasets with same value (%.4f)\n', length(dataset_group.indices), dataset_value);
                
                % Merge detected area number data
                merged_number_data = merge_dataset_group(dataset_group.number_data, dataset_value);
                all_number_data{end+1} = merged_number_data;
                
                % Merge detected area size data
                merged_area_data = merge_dataset_group(dataset_group.area_data, dataset_value);
                all_area_data{end+1} = merged_area_data;
                
                % Merge aspect ratio data
                merged_aspect_data = merge_dataset_group(dataset_group.aspect_data, dataset_value);
                all_aspect_data{end+1} = merged_aspect_data;
                
                % Merge average spacing data
                merged_spacing_data = merge_dataset_group(dataset_group.spacing_data, dataset_value);
                all_spacing_data{end+1} = merged_spacing_data;
            else
                % Only one dataset, use directly
                all_number_data{end+1} = dataset_group.number_data;
                all_area_data{end+1} = dataset_group.area_data;
                all_aspect_data{end+1} = dataset_group.aspect_data;
                all_spacing_data{end+1} = dataset_group.spacing_data;
            end
        end
    end
    
    % Create combined plots
    create_combined_plots(all_number_data, all_area_data, ...
                          all_aspect_data, all_spacing_data, ...
                          output_folder, all_dataset_values);
    
    % Save combined data
    save_combined_data(all_number_data, all_area_data, ...
                       all_aspect_data, all_spacing_data, ...
                       output_folder);
end

function field_name = create_valid_fieldname(value)
    % Create valid MATLAB field name
    % Convert numerical value to valid field name string, avoiding decimal points
    
    % Convert value to string with fixed format
    str = sprintf('val_%.6f', value);
    
    % Replace decimal points and other invalid characters
    str = strrep(str, '.', 'p');  % Replace decimal point with 'p'
    str = strrep(str, '-', 'n');  % Replace negative sign with 'n'
    str = strrep(str, '+', 'p');  % Replace positive sign with 'p'
    
    % Ensure field name starts with letter (MATLAB requirement)
    if ~isletter(str(1))
        str = ['v', str];  % If not starting with letter, add 'v' prefix
    end
    
    field_name = str;
end

function merged_data = merge_dataset_group(data_group, dataset_value)
    % Merge dataset groups with same numerical value
    
    if isempty(data_group)
        merged_data = [];
        return;
    end
    
    % Check if all data have same x values
    x_values = data_group(1).x;
    all_same_x = true;
    
    for i = 2:length(data_group)
        if ~isequal(data_group(i).x, x_values)
            all_same_x = false;
            break;
        end
    end
    
    if all_same_x && length(data_group) > 1
        % All datasets have same x values, calculate mean and standard deviation
        all_means = zeros(length(data_group), length(x_values));
        all_stds = zeros(length(data_group), length(x_values));
        
        for i = 1:length(data_group)
            all_means(i, :) = data_group(i).mean(:)';
            all_stds(i, :) = data_group(i).std(:)';
        end
        
        % Calculate average mean and standard deviation
        merged_mean = mean(all_means, 1);
        
        % Merge standard deviations: using error propagation formula
        merged_std = sqrt(sum(all_stds.^2, 1)) / length(data_group);
        
        % Create merged data structure
        merged_data = struct(...
            'name', sprintf('%.4f', dataset_value), ...
            'x', x_values, ...
            'mean', merged_mean(:), ...
            'std', merged_std(:), ...
            'num_datasets', length(data_group));  % Add dataset count information
    else
        % Different x values or only one dataset, use first dataset
        merged_data = data_group(1);
        merged_data.num_datasets = length(data_group);
    end
end

function create_combined_plots(all_number_data, all_area_data, ...
                               all_aspect_data, all_spacing_data, ...
                               output_folder, dataset_values)
    % Create combined plots
    
    % Use jet colormap
    colors = jet(length(dataset_values));
    
    % Define line styles and markers (for distinguishing multiple datasets with same value)
    line_styles = {'-', '--', ':', '-.'};
    markers = {'o', 's', 'd', '^', 'v', '>', '<', 'p', 'h'};
    
    % Create figure window, adjust width to 1.4x original to accommodate legends
    figure('Position', [100, 100, 1800, 900]);
    
    % Subplot 1: Number of detected areas
    subplot(2, 2, 1);
    [h1, legend_labels1] = plot_combined_data(all_number_data, colors, line_styles, markers, dataset_values);
    xlabel('α (deg)');
    ylabel('Detected Area Number');
    title('Number of Detected Areas');
    grid on;
    
    % Subplot 2: Size of detected areas
    subplot(2, 2, 2);
    [h2, legend_labels2] = plot_combined_data(all_area_data, colors, line_styles, markers, dataset_values);
    xlabel('α (deg)');
    ylabel('Detected Area Size');
    title('Size of Detected Areas');
    grid on;
    
    % Subplot 3: Aspect ratio
    subplot(2, 2, 3);
    [h3, legend_labels3] = plot_combined_data(all_aspect_data, colors, line_styles, markers, dataset_values);
    xlabel('α (deg)');
    ylabel('Average Aspect Ratio');
    title('Aspect Ratio');
    grid on;
    
    % Subplot 4: Average spacing
    subplot(2, 2, 4);
    [h4, legend_labels4] = plot_combined_data(all_spacing_data, colors, line_styles, markers, dataset_values);
    xlabel('α (deg)');
    ylabel('Average Spacing');
    title('Average Spacing');
    grid on;
    
    % Add overall title
    sgtitle('Comparative Analysis of Multiple Datasets (Sorted by Numerical Label)', 'FontSize', 16, 'FontWeight', 'bold');
    
    % Adjust layout to make space for legends
    theme("light");
    set(findall(gcf, '-property', 'FontSize'), 'FontSize', 12);
    
    % Add legend to each subplot (placed on the right side of subplot)
    add_legend_to_subplot(1, h1, legend_labels1, 'right');
    add_legend_to_subplot(2, h2, legend_labels2, 'right');
    add_legend_to_subplot(3, h3, legend_labels3, 'right');
    add_legend_to_subplot(4, h4, legend_labels4, 'right');
    
    % Save figure
    save_path = fullfile(output_folder, 'combined_plots.jpg');
    print(gcf, save_path, '-djpeg', '-r300');
    
    % Create individual comparison plots (one plot per metric)
    create_individual_combined_plots(all_number_data, all_area_data, ...
                                     all_aspect_data, all_spacing_data, ...
                                     output_folder, dataset_values);
end

function add_legend_to_subplot(subplot_idx, h_plots, legend_labels, position)
    % Add legend to subplot (placed outside subplot)
    
    % Select current subplot
    subplot(2, 2, subplot_idx);
    
    % Get subplot position
    subplot_pos = get(gca, 'Position');
    
    % Calculate legend position
    switch position
        case 'right'
            % Place legend on the right side of subplot, next to it
            legend_pos = [subplot_pos(1)+subplot_pos(3)+0.01, subplot_pos(2), ...
                          0.1, subplot_pos(4)];
        case 'left'
            % Place legend on the left side of subplot
            legend_pos = [subplot_pos(1)-0.11, subplot_pos(2), ...
                          0.1, subplot_pos(4)];
        otherwise
            % Default to right side
            legend_pos = [subplot_pos(1)+subplot_pos(3)+0.01, subplot_pos(2), ...
                          0.1, subplot_pos(4)];
    end
    
    % Create legend
    legend(h_plots, legend_labels, ...
           'Location', 'none', ...
           'Position', legend_pos, ...
           'Box', 'off', ...
           'FontSize', 10);
end

function [h_plots, legend_labels] = plot_combined_data(all_data, colors, line_styles, markers, dataset_values)
    % Plot combined data and return plot handles and legend labels
    
    hold on;
    
    % Get unique dataset values
    unique_values = dataset_values;
    
    % Preallocate handles and labels arrays
    h_plots = zeros(length(all_data), 1);
    legend_labels = cell(length(all_data), 1);
    
    for i = 1:length(all_data)
        data = all_data{i};
        
        % Find color index corresponding to data
        % Assume dataset name is numerical string
        data_value = str2double(data.name);
        color_idx = find(abs(unique_values - data_value) < 1e-6);
        
        if isempty(color_idx)
            color_idx = mod(i-1, size(colors, 1)) + 1;
        else
            color_idx = color_idx(1);  % Take first match
        end
        
        % Select style
        line_style_idx = mod(i-1, length(line_styles)) + 1;
        marker_idx = mod(i-1, length(markers)) + 1;
        
        % Create legend label
        if isfield(data, 'num_datasets') && data.num_datasets > 1
            legend_label = sprintf('%.4f (%d)', data_value, data.num_datasets);
        else
            legend_label = sprintf('%.4f', data_value);
        end
        legend_labels{i} = legend_label;
        
        % Plot error bar graph
        h = errorbar(data.x, data.mean, data.std, ...
                     'Color', colors(color_idx, :), ...
                     'LineStyle', line_styles{line_style_idx}, ...
                     'Marker', markers{marker_idx}, ...
                     'MarkerSize', 8, ...
                     'MarkerFaceColor', colors(color_idx, :), ...
                     'LineWidth', 2, ...
                     'DisplayName', legend_label);
        
        % Make markers visible
        set(h, 'MarkerFaceColor', colors(color_idx, :));
        
        % Store handle
        h_plots(i) = h;
    end
    
    hold off;
    
    % Set grid
    grid on;
end

function create_individual_combined_plots(all_number_data, all_area_data, ...
                                          all_aspect_data, all_spacing_data, ...
                                          output_folder, dataset_values)
    % Create individual comparison plots (one plot per metric)
    
    colors = jet(length(dataset_values));
    line_styles = {'-', '--', ':', '-.'};
    markers = {'o', 's', 'd', '^', 'v', '>', '<', 'p', 'h'};
    
    % 1. Number of detected areas comparison plot
    figure('Position', [100, 100, 1000, 600]);  % Increase width to accommodate legend
    [h1, legend_labels1] = plot_combined_data(all_number_data, colors, line_styles, markers, dataset_values);
    xlabel('α (deg)');
    ylabel('Detected Area Number');
    title('Number of Detected Areas - Multiple Datasets');
    grid on;
    theme("light");
    set(gca, 'FontSize', 14);
    
    % Adjust current axis position to make space for legend
    ax_pos = get(gca, 'Position');
    set(gca, 'Position', [ax_pos(1), ax_pos(2), ax_pos(3)*0.75, ax_pos(4)]);
    
    % Add legend (placed on the right side)
    legend(h1, legend_labels1, ...
           'Location', 'eastoutside', ...
           'Box', 'off', ...
           'FontSize', 12, ...
           'NumColumns', 1);  % Set to single column
    
    save_path = fullfile(output_folder, 'combined_number.jpg');
    print(gcf, save_path, '-djpeg', '-r300');
    close(gcf);
    
    % 2. Size of detected areas comparison plot
    figure('Position', [100, 100, 1000, 600]);
    [h2, legend_labels2] = plot_combined_data(all_area_data, colors, line_styles, markers, dataset_values);
    xlabel('α (deg)');
    ylabel('Detected Area Size');
    title('Size of Detected Areas - Multiple Datasets');
    grid on;
    theme("light");
    set(gca, 'FontSize', 14);
    
    % Adjust current axis position to make space for legend
    ax_pos = get(gca, 'Position');
    set(gca, 'Position', [ax_pos(1), ax_pos(2), ax_pos(3)*0.75, ax_pos(4)]);
    
    % Add legend (placed on the right side)
    legend(h2, legend_labels2, ...
           'Location', 'eastoutside', ...
           'Box', 'off', ...
           'FontSize', 12, ...
           'NumColumns', 1);
    
    save_path = fullfile(output_folder, 'combined_area.jpg');
    print(gcf, save_path, '-djpeg', '-r300');
    close(gcf);
    
    % 3. Aspect ratio comparison plot
    figure('Position', [100, 100, 1000, 600]);
    [h3, legend_labels3] = plot_combined_data(all_aspect_data, colors, line_styles, markers, dataset_values);
    xlabel('α (deg)');
    ylabel('Average Aspect Ratio');
    title('Aspect Ratio - Multiple Datasets');
    grid on;
    theme("light");
    set(gca, 'FontSize', 14);
    
    % Adjust current axis position to make space for legend
    ax_pos = get(gca, 'Position');
    set(gca, 'Position', [ax_pos(1), ax_pos(2), ax_pos(3)*0.75, ax_pos(4)]);
    
    % Add legend (placed on the right side)
    legend(h3, legend_labels3, ...
           'Location', 'eastoutside', ...
           'Box', 'off', ...
           'FontSize', 12, ...
           'NumColumns', 1);
    
    save_path = fullfile(output_folder, 'combined_aspect.jpg');
    print(gcf, save_path, '-djpeg', '-r300');
    close(gcf);
    
    % 4. Average spacing comparison plot
    figure('Position', [100, 100, 1000, 600]);
    [h4, legend_labels4] = plot_combined_data(all_spacing_data, colors, line_styles, markers, dataset_values);
    xlabel('α (deg)');
    ylabel('Average Spacing');
    title('Average Spacing - Multiple Datasets');
    grid on;
    theme("light");
    set(gca, 'FontSize', 14);
    
    % Adjust current axis position to make space for legend
    ax_pos = get(gca, 'Position');
    set(gca, 'Position', [ax_pos(1), ax_pos(2), ax_pos(3)*0.75, ax_pos(4)]);
    
    % Add legend (placed on the right side)
    legend(h4, legend_labels4, ...
           'Location', 'eastoutside', ...
           'Box', 'off', ...
           'FontSize', 12, ...
           'NumColumns', 1);
    
    save_path = fullfile(output_folder, 'combined_spacing.jpg');
    print(gcf, save_path, '-djpeg', '-r300');
    close(gcf);
    
    % 5. Create colorbar legend (showing value-color correspondence)
    create_colorbar_legend(dataset_values, colors, output_folder);
end

function create_colorbar_legend(dataset_values, colors, output_folder)
    % Create colorbar legend showing correspondence between values and colors
    
    figure('Position', [100, 100, 400, 500]);
    
    % Normalize values to [0,1] range
    if length(dataset_values) > 1
        normalized_values = (dataset_values - min(dataset_values)) / (max(dataset_values) - min(dataset_values));
    else
        normalized_values = 0.5;  % Only one value, place in middle
    end
    
    % Create colorbar
    colormap(jet);
    c = colorbar('eastoutside');
    c.Label.String = 'Dataset Value';
    c.Ticks = normalized_values;
    
    % Set tick labels to actual values
    tick_labels = cell(length(dataset_values), 1);
    for i = 1:length(dataset_values)
        tick_labels{i} = sprintf('%.4f', dataset_values(i));
    end
    c.TickLabels = tick_labels;
    
    % Add title
    title('Dataset Values - Color Mapping');
    axis off;
    
    % Save colorbar legend
    save_path = fullfile(output_folder, 'colorbar_legend.jpg');
    print(gcf, save_path, '-djpeg', '-r300');
    close(gcf);
end

function save_combined_data(all_number_data, all_area_data, ...
                            all_aspect_data, all_spacing_data, ...
                            output_folder)
    % Save combined data
    
    % Create data structure
    combined_data = struct();
    
    % Number of detected areas data
    combined_data.detected_number = all_number_data;
    
    % Size of detected areas data
    combined_data.detected_area = all_area_data;
    
    % Aspect ratio data
    combined_data.aspect_ratio = all_aspect_data;
    
    % Average spacing data
    combined_data.average_spacing = all_spacing_data;
    
    % Save to file
    save_filename = fullfile(output_folder, 'combined_plot_data.mat');
    save(save_filename, 'combined_data');
    
    fprintf('Combined data saved to: %s\n', save_filename);
end