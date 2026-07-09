function collate_fv_curves

top_data_dir = "../sim_data/test_pCa50";
top_data_dir = "C:/Users/Ken/University of Kentucky/K. Campbell Grants - Submissions - Campbell/temp/test_45"
excel_file_string = "isotonic/sim_output/fv_analysis.xlsx";
field_strings = ["m_force", "m_f_to_f_max", ...
    "m_velocity_l0_per_s", "m_power", "m_rel_power"];
field2_strings = [...
    "fv_x_0", "fv_a", "fv_b", "fv_v_max", "fv_v_max_l0_per_s", ...
    "rel_fv_x_0", "rel_fv_a", "rel_fv_b", ...
    "pow_x_0", "pow_a", "pow_b", "x_at_max_power", ...
    "rel_pow_x_0", "rel_pow_a", "rel_pow_b", "x_at_max_rel_power"]
ofs = "../output/test_pCa50/summary.xlsx";

nn = 10;


data_dirs = return_sub_folders(top_data_dir)

output = [];
output2 = [];

figure(1);
counter  = 1;

waitbar(0);

for i = 1 : numel(data_dirs)

    % if (i > nn)
    %     break
    % end

    waitbar(i / numel(data_dirs))

    try
        dfs = fullfile(data_dirs(i), excel_file_string);
        d = readtable(dfs);
        dn = d.Properties.VariableNames';

        d2 = readtable(dfs, Sheet = 'curve_1');
    catch
        continue
    end

    for j = 1 : numel(field_strings)
        output.(field_strings{j}){counter} = d.(field_strings{j});
    end

    for j = 1 : numel(field2_strings) 
        output2.(field2_strings(j))(counter) = d2.(field2_strings(j))(1);
    end

    parts = split(dfs, filesep);
    output2.sample(counter) = parts(end-3);

    counter = counter + 1;
    
end

output2 = columnize_structure(output2);
output2 = struct2table(output2);
try
    delete(ofs);
end
writetable(output2, ofs);

% output = columnize_structure(output);
% output = struct2table(output)
% 
% return

subplots = layout_subplots( ...
    figure_width = 7, ...
    panels_wide = 2, ...
    panels_high = 2, ...
    x_to_y_ratio = 2)

color_map = parula(5);

n = numel(output.(field_strings{1}))
for i = 1 : n

    % if (i > nn)
    %     break
    % end

    f = output.m_f_to_f_max{i};
    v = output.m_velocity_l0_per_s{i};
    p = output.m_rel_power{i};

    cm = color_map(1 + mod(i, size(color_map, 1)), :);

    subplot(subplots(1));
    hold on
    [~,~,~,~, x_fit, y_fit] = fit_hyperbola('x_data', f, 'y_data', v);
    plot(f, v, 'o', Color = cm);
    plot(x_fit, y_fit, '-', Color = cm);

    subplot(subplots(2));
    hold on
    [~,~,~,~, x_fit, y_fit] = fit_power_curve(f, p);
    plot(f, p, 'o', Color = cm);
    plot(x_fit, y_fit, '-', Color = cm);


end

for i = 1 : 2
    subplot(subplots(i));
    xlim([0 1])
end

x = output2.x_at_max_rel_power'



