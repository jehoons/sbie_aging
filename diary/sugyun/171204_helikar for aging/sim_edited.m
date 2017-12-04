clear; close all; clc;

% Load model file & Convert data type
load('Relation.mat');
load('basal.mat');
node_names=  {'ABL1', 'AKT1', 'AMPK', 'ATG13', 'ATG5', 'ATM/ATR', 'ATM/ATR_2', 'BAX', 'BCL2', 'CASP3', 'CASP9', 'CDK4', 'CDKN1A', 'CDKN2A', 'CHK1/CHK2', 'DNAdamage', 'DNAdamage_high', 'E2F1', 'EIF4EBP1', 'FOXO3', 'FOXO3_ac', 'HDAC1', 'IGF1', 'IGF1R', 'IKBKB', 'IL1B', 'IL6', 'IRS1', 'MAP2K3/MAP2K6', 'MAPK14', 'MDM2', 'MTOR', 'NAD+', 'NFKB1', 'NFKBIE', 'PDK1', 'PIK3CA', 'PPARGC1A', 'PTEN', 'RB1', 'RHEB', 'ROS', 'S6K1', 'SGK1', 'SIRT1', 'SOD2', 'TNF', 'TP53', 'TP53_s15', 'TSC2', 'ULK1', 'lowNutrition'};
node_map = containers.Map( node_names, 1:length(node_names) );

% Generate initial state ( row vector )
initial_state = zeros( size( node_names ) );
%16,17,23,52
X_index = 52;
% Set options for Helikar simulation
sim_options.simulation_time = 100; % total simulation time, 한 점의 시뮬레이션 수
sim_options.transient_time = 5; % transient time
constraints = struct( ...
    'index', {X_index}, ... % node
    'range', {[0,100]} ... % ranged constraint of each node
);
sim_options.constraints = constraints;
sim_options.constraints_num = 100; % number of different constraints, width를 몇등분할 것인가
%% 
    % Set simulation option
    simulation_time = sim_options.simulation_time;
    transient_time = sim_options.transient_time;
    indices = [sim_options.constraints.index];
    constraints_num = sim_options.constraints_num;
    [rand_constraints, rand_constraint_seqs] = processOption( sim_options );
    %rand_constraints는 Latin hypercube sampling 결과에 rand_constraints( i, :
    %) * width + range(1)하여 constraint걸린 노드의 값을 원하는데로 만듬
    % 여기선 100*3 matrix
    % rand_constraint_seqs = zeros( size(rand_constraints, 1), size(rand_constraints, 2), simulation_time );
    %100*3*1000
    
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
layout = [6, 10];
num_plot = layout(1)*layout(2);
num_fig = ceil( length( node_names ) / num_plot );
node_counter = 1;

for i = 1 : num_fig
    figure( i );
    plot_counter = 1;
    for j = 1 : num_plot
        subplot( 6, 10, plot_counter );
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