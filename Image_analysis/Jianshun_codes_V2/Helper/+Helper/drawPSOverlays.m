function drawPSOverlays(kx, ky, k_min, k_max)
%% drawPSOverlays
% Author:       Karthik
% Date:         2025-09-12
% Version:      1.0
%
% Description:
%   Draw overlays on existing FFT plot.
%   - Radial lines every 30°
%   - Annular highlight with white (upper half) and gray (lower half) circles at k_min and k_max
%   - Horizontal white bands at ky=0 between k_min and k_max
%   - Scale ticks and labels every 1 μm⁻¹ along each radial line
%
% Inputs:
%   kx, ky   - reciprocal space vectors (μm⁻¹)
%   k_min    - inner annulus radius (μm⁻¹)
%   k_max    - outer annulus radius (μm⁻¹)
%
% Notes:
%   Optional notes, references.

    hold on

    % === Overlay Radial Lines + Scales ===
    max_kx = max(abs(kx));
    max_ky = max(abs(ky));

    for angle = 0 : pi/6 : pi
        x_line = [0, max_kx] * cos(angle);
        y_line = [0, max_ky] * sin(angle);

        % Plot radial lines
        plot(x_line,  y_line, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2);
        plot(x_line, -y_line, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2);

        % Draw scale ticks along both lines
        drawTicksAlongLine(0,0, x_line(2),  y_line(2));
        drawTicksAlongLine(0,0, x_line(2), -y_line(2));
    end

    % === Overlay Annular Highlight ===
    theta_full = linspace(0, 2*pi, 500);

    % Upper half: white dashed circles
    plot(k_min * cos(theta_full(theta_full <= pi)), ...
         k_min * sin(theta_full(theta_full <= pi)), 'k--', 'LineWidth', 1.2);
    plot(k_max * cos(theta_full(theta_full <= pi)), ...
         k_max * sin(theta_full(theta_full <= pi)), 'k--', 'LineWidth', 1.2);

    % Lower half: gray dashed circles
    plot(k_min * cos(theta_full(theta_full > pi)), ...
         k_min * sin(theta_full(theta_full > pi)), '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.0);
    plot(k_max * cos(theta_full(theta_full > pi)), ...
         k_max * sin(theta_full(theta_full > pi)), '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.0);

    % === Highlight horizontal band across k_y = 0 ===
    x_vals = kx;
    xW1 = x_vals((x_vals >= -k_max) & (x_vals < -k_min));
    xW2 = x_vals((x_vals >  k_min) & (x_vals <= k_max));

    plot(xW1, zeros(size(xW1)), 'k--', 'LineWidth', 1.2);
    plot(xW2, zeros(size(xW2)), 'k--', 'LineWidth', 1.2);

    hold off


    % --- Nested helper function to draw ticks along a radial line ---
    function drawTicksAlongLine(x_start, y_start, x_end, y_end)
        % Tick parameters
        tick_spacing = 1;      % spacing between ticks in μm⁻¹
        tick_length  = 0.05 * sqrt((x_end - x_start)^2 + (y_end - y_start)^2);
        tick_color   = [0.5 0.5 0.5];
        font_size    = 8;

        % Vector along the line
        dx = x_end - x_start;
        dy = y_end - y_start;
        L  = sqrt(dx^2 + dy^2);
        ux = dx / L;
        uy = dy / L;

        % Perpendicular vector for ticks
        perp_ux = -uy;
        perp_uy =  ux;

        % Number of ticks
        n_ticks = floor(L / tick_spacing);

        for i = 1:n_ticks
            xt = x_start + i * tick_spacing * ux;
            yt = y_start + i * tick_spacing * uy;

            % Tick endpoints
            xt1 = xt - 0.5 * tick_length * perp_ux;
            yt1 = yt - 0.5 * tick_length * perp_uy;
            xt2 = xt + 0.5 * tick_length * perp_ux;
            yt2 = yt + 0.5 * tick_length * perp_uy;

            % Draw tick
            plot([xt1 xt2], [yt1 yt2], '-', 'Color', tick_color, 'LineWidth', 1);

            % Label
            text(xt, yt, sprintf('%d', i), ...
                'Color', tick_color, ...
                'FontSize', font_size, ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'bottom', ...
                'Rotation', atan2d(dy, dx));
        end
    end
end