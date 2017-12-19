
function [rand_constraints, rand_constraint_seqs] = processOption( sim_options )


    simulation_time = sim_options.simulation_time;
    transient_time = sim_options.transient_time;
    constraints = sim_options.constraints;
    constraints_num = sim_options.constraints_num;



%     simulation_time = 300;
%     transient_time = 100;
%     constraints = struct( ...
%         'index', {1, 2, 3}, ...
%         'range', {[0, 100], [20,60], [10,20]} ...
%     );
%     constraints_num = 10;



    indices = [constraints.index];
    ranges = {constraints.range};
    rand_constraints = lhsdesign( constraints_num, length( indices ) )';

    for i = 1 : length( indices )
        range = ranges{ i };
        width = range(2) - range(1) + 1;

        rand_constraints( i, : ) = rand_constraints( i, : ) * width + range(1);
    end;

    rand_constraints = floor( rand_constraints );


    block_count = simulation_time/100;

    rand_constraint_seqs = zeros( size(rand_constraints, 1), size(rand_constraints, 2), simulation_time );
    for i = 1 : length( indices )
        for j = 1 : constraints_num
            on_ratio = rand_constraints( i, j );
            for k = 0 : block_count-1
                rand_constraint_seqs( i, j, k*100+randperm( 100, on_ratio ) ) = 1;
            end;
        end;
    end;


end



