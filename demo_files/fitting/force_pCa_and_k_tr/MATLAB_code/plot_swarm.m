function plot_swarm

dfs = '../progress/progress.xlsx'
no_of_particles = 6;
no_of_parameters = 3;


% Read in the progress table
d = readtable(dfs);
dn = d.Properties.VariableNames

% Deduce number of iterations
no_of_iterations = floor(size(d,1) / no_of_particles);

% Setup
particle_best_value = inf * ones(no_of_particles, 1);
particle_best_x = NaN * ones(no_of_particles, no_of_parameters);

global_best_value = inf;
global_best_x = NaN * ones(1, no_of_parameters);

% Make a figure
figure(1)
subplots = initialise_publication_quality_figure( ...
                'no_of_panels_wide', no_of_parameters-1, ...
                'no_of_panels_high', no_of_parameters-1, ...
                'right_margin', 1, ...
                'axes_padding_left', 0.25, ...
                'axes_padding_right', 0.25, ...
                'x_to_y_axes_ratio', 2);

% Loop through iterations
for iter = 1 : no_of_iterations

    for row = 1 : (no_of_parameters - 1)
        for col = 1 : row

            subplot(subplots(col + (row-1)*(no_of_parameters-1)));
            cla
            hold on;

            x_string = sprintf('p_%i', col);
            y_string = sprintf('p_%i', row+1);

            for part = 1 : no_of_particles
                idx = (iter-1)*no_of_particles + part;

                plot(d.(x_string)(idx), d.(y_string)(idx), 'bo');

                if ((row == 1) && (col==1))
                    if (d.error_total(idx) < particle_best_value(part))
                        particle_best_value(part) = d.error_total(idx);
                        for i = 1 : no_of_parameters
                            particle_best_x(part, i) = d.(sprintf('p_%i', i))(idx);
                        end
                    end
    
                    if (d.error_total(idx) < global_best_value)
                        global_best_value = d.error_total(idx);
                        for i = 1 : no_of_parameters
                            global_best_x(i) = d.(sprintf('p_%i', i))(idx);
                        end
                    end
                end

                plot(particle_best_x(part, col), particle_best_x(part, row+1), 'go');
            end
            plot(global_best_x(col), global_best_x(row+1), 'rp')

            xlim([0, 1]);
            ylim([0, 1]);
        end
    end

    drawnow;
    pause(0.2);
    
end
