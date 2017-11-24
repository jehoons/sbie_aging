clear; close all; clc;

% Load model file & Convert data type
load('Relation_emt4.mat');
load('Basal_emt4.mat');
node_names=  {'DNA_damage','ATM','EGF','EGFR','TGFb','TGFbR','TNFa','TNFR','WNT','NOTCH1','smad2/3/4','Grb2/Sos','Ras','Raf','MEK','ERK','MYC','PI3K','PTEN','AKT','GSK3-b','b-catenin_nuclear','b-catenin','EML4_ALK','STAT3','RHOA','ROCK','NFKB','IKB','miR200','miR34','SLUG','SNAIL','TWIST1','ZEB1/2','FOXA2','RB1','Mdm2','TP53','P21','cyclinD','cyclinE','E-cadherin','VIM','Nanog'};
node_map = containers.Map( node_names, 1:length(node_names) );

% Generate initial state ( row vector )
initial_state = zeros( size( node_names ) );

X_index =5;
% Set options for Helikar simulation
sim_options.simulation_time = 1000; % total simulation time
sim_options.transient_time = 700; % transient time
constraints = struct( ...
    'index', {5}, ... % node
    'range', {[0,100]} ... % ranged constraint of each node
);
sim_options.constraints = constraints;
sim_options.constraints_num = 100; % number of different constraints
%% 
    % Set simulation option
    simulation_time = sim_options.simulation_time;
    transient_time = sim_options.transient_time;
    indices = [sim_options.constraints.index];
    constraints_num = sim_options.constraints_num;
    [rand_constraints, rand_constraint_seqs] = processOption( sim_options );

    result = nan( constraints_num, length( initial_state ) );
% Run simulation with the options
for i = 1 : constraints_num
    input_seqs = squeeze( rand_constraint_seqs( :, i, : ) );

    trajectory = nan( simulation_time, length( initial_state ) );

    current_state = initial_state;
    for j = 1 : simulation_time
        if size(input_seqs,2)==1
            current_state( indices ) = input_seqs( j );
        else
            current_state( indices ) = input_seqs( :, j );
        end;
        weight_sum = transpose(Relation * transpose(current_state)+Basal); % Change Cell Type
        updated_state = current_state;
        updated_state(weight_sum>0)=1;
        updated_state(weight_sum<0)=0;
        if size(input_seqs,2)==1
            updated_state( indices ) = input_seqs( j );
        else
            updated_state( indices ) = input_seqs( :, j );
        end;
        trajectory( j, : ) = updated_state;
        current_state = updated_state;
    end;

    result( i, : ) = sum( trajectory( transient_time+1 : end, : ) ) / (simulation_time-transient_time);
end;

%X_index = 5;% Change X_index
layout = [6, 8];
num_plot = layout(1)*layout(2);
num_fig = ceil( length( node_names ) / num_plot );
node_counter = 1;

for i = 1 : num_fig
    figure( i );
    plot_counter = 1;
    for j = 1 : num_plot
        subplot( 6, 8, plot_counter );
        plot( result(:,X_index), result(:,node_counter), 'o' );
        title( strrep( node_names( node_counter ), '_', '\_' ) );
        xlim( [0, 1] ); ylim( [0, 1] );

        plot_counter = plot_counter + 1;

        node_counter = node_counter + 1;
        if( node_counter > length( node_names ) )
            break;
        end;

        if( plot_counter > num_plot )
            break;
        end;
    end;
end;