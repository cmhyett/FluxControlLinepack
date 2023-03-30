using Pkg;
using GasTranSim;
using Plots;
using CSV, JSON;
using Tables;
using Serialization;
using Printf;
using Dierckx;
basepath = pwd() * "/";
include(basepath * "/src/csvToJSON.jl");
include(basepath * "/src/plotting.jl");

function runScenarios(scenNumArr::Array{Int,1})
    tsResults = [];
    for num in scenNumArr
        csvSpecPath = basepath * "examples/israelNetwork/reducedNetwork/csvSpecification/";
        jsonSpecPath = csvSpecPath * "../jsonSpec_$(num)/";
        mkpath(jsonSpecPath);
        writeNetworkJSON(csvSpecPath * "/network_data/pipes.csv",
                         csvSpecPath * "/network_data/nodes.csv",
                         csvSpecPath * "/network_data/compressors.csv",
                         jsonSpecPath * "/network.json");

        writeICJSON(csvSpecPath * "/initial_conditions/pipe_ic.csv",
                    csvSpecPath * "/initial_conditions/node_ic.csv",
                    csvSpecPath * "/initial_conditions/compressor_ic.csv",
                    jsonSpecPath * "/ic.json");

        writeBCJSON(csvSpecPath * "/boundary_condtions/withdrawal_node_bc_scen$(num).csv",
                    csvSpecPath * "/boundary_condtions/slack_node_bc.csv",
                    csvSpecPath * "/boundary_condtions/compressor_bc.csv",
                    jsonSpecPath * "/bc.json");

        writeParamsJSON(csvSpecPath * "/network_data/params_scen$(num).csv",
                        jsonSpecPath * "/params.json");

        writeDisruptionsJSON(jsonSpecPath * "/disruptions.json");

        folder = jsonSpecPath;

        ts = initialize_simulator(folder; eos=:full_cnga);
        run_simulator!(ts);
#        plot_pressure_profile(basepath * "/examples/israelNetwork/reducedNetwork/results/scenario$(num)/scen$(num).gif", ts);
        push!(tsResults, ts);
        println("finished running scenario $(num)");
    end
    return tsResults;
end

# sigma,mu percentage of average mass flux at node,
#  so mu=0, sigma=0.1 is N(0, 0.1*mean(ts.bc[:node][i]))
function monteCarlo(scenNum::Int, numRuns; mu=0.0, sigma=0.1)
    tsResults = Dict();
    csvSpecPath = basepath * "examples/israelNetwork/reducedNetwork/csvSpecification/";
    jsonSpecPath = csvSpecPath * "../jsonSpec_scen$(scenNum)/";
    mkpath(jsonSpecPath);
    writeNetworkJSON(csvSpecPath * "/network_data/pipes.csv",
                     csvSpecPath * "/network_data/nodes.csv",
                     csvSpecPath * "/network_data/compressors.csv",
                     jsonSpecPath * "/network.json");

    writeICJSON(csvSpecPath * "/initial_conditions/pipe_ic.csv",
                csvSpecPath * "/initial_conditions/node_ic.csv",
                csvSpecPath * "/initial_conditions/compressor_ic.csv",
                jsonSpecPath * "/ic.json");

    writeBCJSON(csvSpecPath * "/boundary_condtions/withdrawal_node_bc_scen$(scenNum).csv",
                csvSpecPath * "/boundary_condtions/slack_node_bc.csv",
                csvSpecPath * "/boundary_condtions/compressor_bc.csv",
                jsonSpecPath * "/bc.json");

    writeParamsJSON(csvSpecPath * "/network_data/params_scen$(scenNum).csv",
                    jsonSpecPath * "/params.json");

    writeDisruptionsJSON(jsonSpecPath * "/disruptions.json");

    folder = jsonSpecPath;

    function clip(minVal, maxVal, val)
        return max(min(maxVal, val), minVal)
    end
    
    Threads.@threads for i in 1:numRuns
        println("Starting run number $(i)");
        local ts = initialize_simulator(folder; eos=:simple_cnga);
        for key in keys(ts.boundary_conditions[:node])
            old_spline = ts.boundary_conditions[:node][key]["spl"];
            if (mean(old_spline.c) > 0) #demand node
                localMean = (rand()-0.5)*mu*mean(old_spline.c);
                localStd = mean(old_spline.c) * sigma;
                #randVals = clip.(localMean - 2*localStd, localMean + 2*localStd, randn(length(old_spline.c)) .* localStd .+ localMean);
                # if (i == 1) #do max perturbation on all nodes
                #     randVals = (ones(length(old_spline.c)) .- 0.5) .* localStd .+ (localMean);
                # elseif (i == 2) #do -max perturbation on all nodes
                #     randVals = (zeros(length(old_spline.c)) .- 0.5) .* localStd .+ (localMean);
                # else #random draw
                    randVals = (rand(length(old_spline.c)) .- 0.5) .* localStd .+ (localMean);
                # end
                new_vals = old_spline.c .+ randVals;
                ts.boundary_conditions[:node][key]["spl"] =
                    Spline1D(old_spline.t[2:end-1], new_vals, k=1);
            end
        end
        run_simulator!(ts);
        push!(tsResults, i=>ts);
        println("finished running run number: $(i)");
    end

    return tsResults;
end

