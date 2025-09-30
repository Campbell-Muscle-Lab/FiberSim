function fv_curves

% Variables
data_folders = ["../sim_data/release", ...
    "../sim_data/afterload"];

t_ticks = [0 0.25];
fit_t = [0.175 0.21 ; 0.175 0.21];

panels_wide = 2;
panels_high = 7;

% Make a figure
sp = layout_subplots( ...
        panels_wide = panels_wide, ...
        panels_high = panels_high, ...
        figure_width = 7, ...
        x_to_y_ratio = 4);

% Get the color_map
cm = return_matplotlib_default_colors();

% Loop through folders
for fold_i = 1 : numel(data_folders)

    sim_files = findfiles('txt', ...
        fullfile(data_folders(fold_i), "sim_output"), 1)'


    for file_i = 1 : numel(sim_files)

        d = readtable(sim_files{file_i});
        dn = d.Properties.VariableNames'

        ti = find((d.time > t_ticks(1)) & ...
                (d.time < t_ticks(end)));

        subplot(sp(fold_i));
        hold on;
        plot(d.time(ti), 10.^(-d.hs_1_pCa(ti)), '-');

        subplot(sp(fold_i + panels_wide));
        hold on;
        plot(d.time(ti), d.hs_1_force(ti), '-');

        subplot(sp(fold_i + 2*panels_wide));
        hold on;
        plot(d.time(ti), d.hs_1_length(ti), '-');

        % Fit
        fi = find((d.time > fit_t(fold_i, 1)) & ...
                (d.time < fit_t(fold_i, 2)));

        f(fold_i, file_i) = mean(d.hs_1_force(fi), 1);

        fit_line = fit_straight_line(d.time(fi), d.hs_1_length(fi));
        v(fold_i, file_i) = -fit_line.slope;

        plot(fit_line.x_fit, fit_line.y_fit, 'k-');

        % Thin
        subplot(sp(fold_i + (3*panels_wide)));
        for pop_i = 1 : 2
            pop_string = sprintf('hs_1_a_pop_%i', pop_i);
            plot(d.time(ti), d.(pop_string)(ti), '-', ...
                Color = cm(pop_i,:));
        end

        % Thick
        subplot(sp(fold_i + (4*panels_wide)));
        for pop_i = 1 : 4
            pop_string = sprintf('hs_1_m_pop_%i', pop_i);
            plot(d.time(ti), d.(pop_string)(ti), '-', ...
                Color = cm(pop_i,:));
        end

        % Mybpc
        subplot(sp(fold_i + (5*panels_wide)));
        for pop_i = 1 : 3
            pop_string = sprintf('hs_1_c_pop_%i', pop_i);
            plot(d.time(ti), d.(pop_string)(ti), '-', ...
                Color = cm(pop_i,:));
        end

    end
end

for i = 1 : 2

    vi = find(v(i,:) > 0);

    syms = {'o-', 's-'};

    subplot(sp((panels_high-1)*panels_wide + 1));
    hold on;
    plot(f(i, vi), v(i,vi), syms{i}, ...
        Color = cm(i,:));

    subplot(sp((panels_high-1)*panels_wide + 2));
    hold on;
    plot(f(i, vi), f(i, vi).*v(i, vi), syms{i}, ...
        Color = cm(i,:));
end