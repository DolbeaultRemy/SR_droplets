function drawODOverlays(x1, y1, x2, y2)
%% drawODOverlays
% Author:       Karthik
% Date:         2025-09-12
% Version:      1.0
%
% Description:
%   Brief description of the script functionality.
%
% Notes:
%   Optional notes, references.

    % Parameters
    n_ticks      = 5;        % Fixed number of ticks along the diagonal
    tick_length  = 1;        % µm tick mark length
    line_color   = [0.5 0.5 0.5];
    tick_color   = [0.5 0.5 0.5];
    font_size    = 10;

    % Vector from start to end
    dx           = x2 - x1;
    dy           = y2 - y1;
    L            = sqrt(dx^2 + dy^2);

    % Unit direction vector along diagonal
    ux           = dx / L;
    uy           = dy / L;

    % Perpendicular unit vector for ticks
    perp_ux      = -uy;
    perp_uy      = ux;

    % Midpoint (center)
    xc           = (x1 + x2) / 2;
    yc           = (y1 + y2) / 2;

    % Dynamic tick spacing
    tick_spacing = L / (2 * n_ticks);

    % Draw main diagonal line
    plot([x1 x2], [y1 y2], '--', 'Color', line_color, 'LineWidth', 1.2);

    for i = -n_ticks+1:n_ticks-1   % exclude endpoints
        d = i * tick_spacing;
        xt = xc + d * ux;
        yt = yc + d * uy;

        % Tick line endpoints
        xt1 = xt - 0.5 * tick_length * perp_ux;
        yt1 = yt - 0.5 * tick_length * perp_uy;
        xt2 = xt + 0.5 * tick_length * perp_ux;
        yt2 = yt + 0.5 * tick_length * perp_uy;

        % Draw tick
        plot([xt1 xt2], [yt1 yt2], '--', 'Color', tick_color, 'LineWidth', 1);

        % Label: centered at tick, offset slightly along diagonal
        if i > -n_ticks && i < n_ticks && i ~= 0
            text(xt, yt, sprintf('%+d', round(d)), ...
                'Color', tick_color, ...
                'FontSize', font_size, ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'bottom', ...
                'Rotation', atan2d(dy, dx));
        end
    end
end
