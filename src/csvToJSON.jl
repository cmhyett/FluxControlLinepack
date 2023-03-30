# this file transforms a scenario detailed via CSV files
#  into JSON scenario

function writeNetworkJSON(pipeFilepath::String,
                          nodeFilepath::String,
                          compressorFilepath::String,
                          destFilepath::String)
    pipes = CSV.File(pipeFilepath, comment="#");
    nodes = CSV.File(nodeFilepath, comment="#");
    compressors = CSV.File(compressorFilepath, comment="#");
    
    pipeDict = Dict{String, Any}();
    nodeDict = Dict{String, Any}();
    compressorDict = Dict{String, Any}();
    
    for i in 1:length(pipes)
        p = Dict("length"=>pipes.length[i],
                 "friction_factor"=>pipes.friction_factor[i],
                 "disc_seg"=>0,
                 "pipe_name"=>"p$(pipes.number[i])",
                 "pipe_id"=>pipes.number[i],
                 "to_node"=>pipes.end_node[i],
                 "from_node"=>pipes.start_node[i],
                 "area"=>(pipes.diameter[i]/2)^2*π,
                 "diameter"=>pipes.diameter[i])
        push!(pipeDict,
              "$(pipes.number[i])"=>p);
    end

    for i in 1:length(nodes)
        n = Dict("node_id"=>nodes.number[i],
                 "node_name"=>nodes.name[i],
                 "x_coord"=>nodes.lon[i],
                 "y_coord"=>nodes.lat[i],
                 "min_pressure"=>nodes.min_pressure[i],
                 "max_pressure"=>nodes.max_pressure[i],
                 "min_injection"=>nodes.min_injection[i],
                 "max_injection"=>nodes.max_injection[i],
                 "slack_bool"=>nodes.type[i])
        push!(nodeDict,
              "$(nodes.number[i])"=>n);
    end

    for i in 1:length(compressors)
        c = Dict("comp_id"=>compressors.number[i],
                 "comp_name"=>compressors.name[i],
                 "c_max"=>compressors.c_max[i],
                 "c_min"=>compressors.c_min[i],
                 "to_node"=>compressors.to_node[i],
                 "from_node"=>compressors.from_node[i],
                 "max_flow"=>compressors.max_flow[i],
                 "min_flow"=>compressors.min_flow[i],
                 "max_power"=>compressors.max_power[i]);
        push!(compressorDict,
              "$(compressors.number[i])"=>c);
    end

    nw = Dict("pipes"=>pipeDict,
              "nodes"=>nodeDict,
              "compressors"=>compressorDict);

    open(destFilepath, "w") do f
        JSON.print(f, nw, 4);
    end
    return destFilepath;
end

function writeICJSON(pipeFilepath::String,
                     nodeFilepath::String,
                     compressorFilepath::String,
                     destFilepath::String)
    pipes = CSV.File(pipeFilepath, comment="#");
    nodes = CSV.File(nodeFilepath, comment="#");
    compressors = CSV.File(compressorFilepath, comment="#");

    initCompressorFlow = Dict{String, Any}();
    initCompressorPIn = Dict{String, Any}();
    initCompressorPOut = Dict{String, Any}();
    initCRatio = Dict{String, Any}();
    for i in 1:length(compressors)
        number_string = "$(compressors.number[i])";
        push!(initCompressorFlow, number_string=>compressors.initial_compressor_flow[i]);
        push!(initCompressorPIn, number_string=>compressors.initial_compressor_pressure_in[i]);
        push!(initCompressorPOut, number_string=>compressors.initial_compressor_pressure_out[i]);
        push!(initCRatio, number_string=>compressors.initial_compressor_ratio[i]);
    end

    initNodeFlow = Dict{String, Any}();
    initNodePressure = Dict{String, Any}();
    for i in 1:length(nodes)
        nodeNum_string = "$(nodes.number[i])";
        push!(initNodeFlow, nodeNum_string=>nodes.initial_nodal_flow[i]);
        push!(initNodePressure, nodeNum_string=>nodes.initial_nodal_pressure[i]);
    end

    initPipeFlow = Dict{String, Any}();
    initPipePIn = Dict{String, Any}();
    initPipePOut = Dict{String, Any}();
    for i in 1:length(pipes)
        pipeNum_string = "$(pipes.number[i])";
        push!(initPipeFlow,
              pipeNum_string=>pipes.initial_pipe_flow[i]);
        push!(initPipePIn,
              pipeNum_string=>pipes.initial_pipe_pressure_in[i]);
        push!(initPipePOut,
              pipeNum_string=>pipes.initial_pipe_pressure_out[i]);
    end

    icData = Dict("initial_compressor_flow"=>initCompressorFlow,
                  "initial_compressor_pressure_in"=>initCompressorPIn,
                  "initial_compressor_pressure_out"=>initCompressorPOut,
                  "initial_compressor_ratio"=>initCRatio,
                  "initial_nodal_flow"=>initNodeFlow,
                  "initial_nodal_pressure"=>initNodePressure,
                  "initial_pipe_flow"=>initPipeFlow,
                  "initial_pipe_pressure_in"=>initPipePIn,
                  "initial_pipe_pressure_out"=>initPipePOut);

    open(destFilepath, "w") do file
        JSON.print(file, icData, 4);
    end
    return destFilepath;

end

function writeBCJSON( withdrawalNodeBCFilepath::String,
                      slackNodeBCFilepath::String,
                      compressorBCFilepath::String,
                      destFilepath::String)
    slackNodeBC = CSV.File(slackNodeBCFilepath, comment="#", types=String, header=false);
    withNodeBC = CSV.File(withdrawalNodeBCFilepath, comment="#", types=String, header=false);
    compressorBC = CSV.File(compressorBCFilepath, comment="#", types=String, header=false);

    boundaryCompressor = Dict{String, Any}();
    for i in 2:CSV.getcols(compressorBC)
        compDict = Dict{String, Any}();
        push!(compDict,
              "time"=>parse.(Float64,
                             CSV.getcolumn(compressorBC,1).column[2:end]));
        push!(compDict,
              "control_type"=>zeros(Int, CSV.getrows(compressorBC)-1));
        push!(compDict,
              "value"=>parse.(Float64, CSV.getcolumn(compressorBC,i).column[2:end]));
        push!(boundaryCompressor,
              compressorBC[1][i]=>compDict);
    end

    withNodeFlow = Dict{String, Any}();
    for i in 2:CSV.getcols(withNodeBC)
        nodeDict = Dict{String, Any}();
        push!(nodeDict,
              "time"=>parse.(Float64,
                             CSV.getcolumn(withNodeBC,1).column[2:end]));
        push!(nodeDict,
              "value"=>parse.(Float64, CSV.getcolumn(withNodeBC,i).column[2:end]));
        push!(withNodeFlow,
              withNodeBC[1][i]=>nodeDict);
    end

    slackNodeP = Dict{String, Any}();
    for i in 2:CSV.getcols(slackNodeBC)
        nodeDict = Dict{String, Any}();
        push!(nodeDict,
              "time"=>parse.(Float64,
                             CSV.getcolumn(slackNodeBC,1).column[2:end]));
        push!(nodeDict,
              "value"=>parse.(Float64, CSV.getcolumn(slackNodeBC,i).column[2:end]));
        push!(slackNodeP,
              slackNodeBC[1][i]=>nodeDict);
    end

    bcData = Dict("boundary_compressor"=>boundaryCompressor,
                  "boundary_nonslack_flow"=>withNodeFlow,
                  "boundary_pslack"=>slackNodeP);

    open(destFilepath, "w") do file
        JSON.print(file, bcData, 4);
    end
    return destFilepath;
end

function writeParamsJSON(paramsFilepath::String,
                         destFilepath::String)
    params = CSV.File(paramsFilepath, comment="#");
    paramsDict = Dict("Temperature (K):"=>params.temperature[1],
                      "Gas specific gravity (G):"=>params.gas_specific_gravity[1],
                      "Specific heat capacity ratio"=>params.specific_heat_capacity_ratio[1],
                      "units (SI=0, standard=1)"=>0,
                      "Initial time"=>params.initial_time[1],
                      "Final time"=>params.final_time[1],
                      "Discretization time step"=>params.discretization_time_step[1],
                      "Courant number (must be between 0 and 1, recommended value is 0.9)"=>params.courant_number[1],
                      "Output dt"=>params.output_dt[1],
                      "Output dx"=>params.output_dx[1],
                      "Save final state"=>params.save_final_state[1]);

    fullDict = Dict("simulation_params"=>paramsDict);

    open(destFilepath, "w") do file
        JSON.print(file, fullDict, 4);
    end
    return destFilepath;    
                      
end

function writeDisruptionsJSON(destFilepath::String)
    disruptions = Dict("node_disruptions"=>Dict(),
                       "pipe_disruptions"=>Dict());
    open(destFilepath, "w") do file
        JSON.print(file, Dict("disruption"=>disruptions), 4);
    end
    return destFilepath;
end
