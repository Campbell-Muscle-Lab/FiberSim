function plot_swarm

dfs = '../progress/progress.xlsx'

d = readtable(dfs)

xi = 1

dn = d.Properties.VariableNames';


c = 0
ec = []
for i = 1 : numel(dn)
    if (startsWith(dn{i}, 'p_'))
        c = c + 1;
        fn = dn{i}
        p = d.(fn)(xi);
        ec = [ec (p-(c*0.1))^2];
    end
end
ec
e_total = sum(ec)

r = d(xi,:)



% 
% n_particles = 10;
% 
% max_rows = size(d,1);
% 
% keep_going = true;
% 
% figure(1);
% clf;
% 
% while (keep_going)
% 
%     for i = 1 : n_particles
%         
%         r = 
% 
