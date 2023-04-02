using Pkg;
using GasTranSim;
using Plots;
using CSV, JSON;
using Tables;
using Serialization;
using Printf;
using Dierckx;
using DifferentialEquations;


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
function monteCarlo(scenNum::Int, numRuns; sigma=0.5, numNeighbors=11)
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

    function moving_average(arr,n)
        @assert isodd(n);
        result = zeros(length(arr));
        n2 = floor(Int, n/2);
        for i in n2+1:length(arr)-n2
            result[i] = sum(arr[i-n2:i+n2])/n;
        end
        for i in 1:n2
            result[i] = (sum(arr[end-n2+i:end]) + sum(arr[1:i+n2]))/n;
        end
        for i in length(arr)-n2+1:length(arr)
            result[i] = (sum(arr[i-n2:end]) + sum(arr[1:i-length(arr)+n2]))/n;
        end
        return result;
    end

    demandNodes = [2,3,4,5,6,7,9,10,11];
    perturbedBCs = Dict();
    begin
        local ts = initialize_simulator(folder; eos=:full_cnga);
        for key in demandNodes
            println("key = $(key)");
            old_spline = ts.boundary_conditions[:node][key]["spl"];
            new_spline = Spline1D(old_spline.t[2:end-1], moving_average(old_spline.c, numNeighbors), k=2);
            t = old_spline.t[1:end-1];
            b = sigma;
            v = var(old_spline.(t))*(1/sigma);
            localSigma = mean(old_spline.c) * b;
            localTheta = 10;
            localSigma = sqrt(2*localTheta*sigma*mean(old_spline.c)); #var = sigma*var(old)
            u0 = new_spline(0);
            f(u,p,t) = localTheta*(new_spline(t)-u);
            g(u,p,t) = localSigma;
            prob = SDEProblem(f, g, u0, (t[1], t[end]));
            function prob_func(prob,i,repeat)
                return prob; #independent noise
                #return remake(prob, seed=i) #same noise for each node
            end
            eprob = EnsembleProblem(prob, prob_func=prob_func);
            sol = solve(eprob, EM(), dt=0.1, saveat=t, trajectories=numRuns);
            push!(perturbedBCs, key=>sol);
        end
    end
    if (perturbedBCs[6][1].retcode != :Success)
        println("setting up boundary conditions failed!");
        return 0;
    end
    Threads.@threads for i in 1:numRuns
        println("Starting run number $(i)");
        local ts = initialize_simulator(folder; eos=:full_cnga);
        for key in demandNodes
            ts.boundary_conditions[:node][key]["spl"].c .= perturbedBCs[key][i].u;
        end
        run_simulator!(ts);
        push!(tsResults, i=>ts);
        println("finished running run number: $(i)");
    end

    return tsResults;
end

