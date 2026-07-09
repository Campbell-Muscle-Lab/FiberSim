function peruse_fv_curves

top_data_dir = "../sim_data/test_pCa50";
excel_file_string = "isotonic/sim_output/fv_analysis.xlsx";
field_strings = ["m_force", "m_f_to_f_max", ...
    "m_velocity_l0_per_s", "m_power", "m_rel_power"];
field2_strings = [...
    "fv_x_0", "fv_a", "fv_b", "fv_v_max", "fv_v_max_l0_per_s", ...
    "rel_fv_x_0", "rel_fv_a", "rel_fv_b", ...
    "pow_x_0", "pow_a", "pow_b", "x_at_max_power", ...
    "rel_pow_x_0", "rel_pow_a", "rel_pow_b", "x_at_max_rel_power"]
ifs = "../output/test_pCa50/summary.xlsx";
pfs = "../generated/setup/parameter_values.xlsx"

% Code
d = readtable(ifs);
d.sample_id = str2double(extract(d.sample, digitsPattern));

p = readtable(pfs);
p.id = (1:size(p,1))';

d = innerjoin(d, p, LeftKey = "sample_id", RightKey = "id");

% Filter
fs = ["fv_v_max_l0_per_s", "x_at_max_rel_power"];
for i = 1 : numel(fs)
    oi = find(isoutlier(d.(fs(i))));
    d(oi, :) = []
end

% Look at v_max and x_max_power simultaneously
z_v_max = zscore(d.fv_v_max_l0_per_s);
z_x_max_p = zscore(d.x_at_max_rel_power);

d.r = hypot(z_v_max, z_x_max_p);
vi = find( (z_v_max > 0) & (z_x_max_p > 0) )
% [~, vim] = max(d.r(vi));
% vi = vi(vim)

mfa = 0.25

figure(4);
clf;
hold on;
scatter(z_v_max, z_x_max_p, 'bo');
scatter(z_v_max(vi), z_x_max_p(vi), ...
    'r', 'filled', ...
    MarkerFaceAlpha = mfa);
xlabel('z score for x at max power')
ylabel('z score for V_{max}')



dn = d.Properties.VariableNames'
d2 = removevars(d, dn([1:15 17 end]));
d2n = d2.Properties.VariableNames'



cn = d2n(3:end);
r = ceil(sqrt(numel(cn)))

sp3 = layout_subplots( ...
    figure_handle = 3, ...
    figure_width = 10, ...
    panels_wide = r, ...
    panels_high = r);

for i = 1 : numel(cn)

    x_field = cn{i};
    y_field = "x_at_max_rel_power";

    x = d2.(x_field);
    y = d2.(y_field);

    subplot(sp3(i));
    hold on;
    scatter(log10(x), y, 'b', 'filled', ...
        MarkerFaceAlpha = mfa)
    scatter(log10(x(vi)), y(vi), 'r', 'filled', ...
        MarkerFaceAlpha = 2*mfa)
    xlabel(x_field)
    ylabel(y_field)

end




 

return



sdfsf


figure(2);
clf
subplots = layout_subplots( ...
    figure_handle = 2, ...
    panels_wide = 2, ...
    panels_high = 2);


x = linspace(0, 1, 100);

for i = 1 : 10

    for j = 1 : 2
        if (j==1)
            vi = i;
        else
            vi = size(d,1) - (i-1);
        end
    
        x_0 = d.rel_fv_x_0(vi);
        a = d.rel_fv_a(vi);
        b = d.rel_fv_b(vi);
        v(vi,:) = return_hyperbola(x, x_0, a, b);
    
        plot_index = 2*(j-1) + 1;
        subplot(subplots(plot_index));
        hold on;
        plot(x, v(vi,:), 'b-');
    
        x_0 = d.rel_pow_x_0(vi);
        a = d.rel_pow_a(vi);
        b = d.rel_pow_b(vi);
        po(vi,:) = return_power_curve(x, x_0, a, b);
    
        plot_index = 2*(j-1) + 2;
        subplot(subplots(plot_index));
        hold on;
        plot(x, po(vi,:), 'b-');
    end

end

d2 = removevars(d, dn([1:15 17]))
d2n = d2.Properties.VariableNames'

cn = d2n(3:end);
r = ceil(sqrt(numel(cn)))

sp3 = layout_subplots( ...
        figure_handle = 3, ...
        figure_width = 7, ...
        panels_wide = r, ...
        panels_high = r);

for i = 1 : numel(cn)
    
    x = d2.(cn{i});
    y = d2.x_at_max_rel_power;
    
    subplot(sp3(i));
    plot(x, y, 'bo')
end
    


return



% 
% (((x0+a).*b)./(x_fit(i)+a))-b;
% return



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

x = output2.x_at_max_rel_power'



end

function y = return_hyperbola(x, x0, a, b)

    for i = 1 : numel(x)
        y(i) = (((x0+a)*b)/(x(i)+a))-b;
    end
end

function y=return_power_curve(x,x0,a,b)
    for i = 1 : numel(x)
        y(i) = x(i)*b*(((x0+a)/(x(i)+a))-1);
    end
end