basepath, scenNum, nRuns, sigma = ARGS;
using Pkg;
Pkg.activate(basepath);
Pkg.instantiate();
using GasTranSim, Serialization, Plots;
include(basepath * "/src/csvToJSON.jl");
include(basepath * "./src/plotting.jl");
include(basepath * "./examples/runExample.jl");

nRuns = parse(Int,nRuns);
scenNum = parse(Int, scenNum);
sigma = parse(Float64, sigma);

outputPath = basepath * "/results/";
mkpath(outputPath);
sol = monteCarlo(scenNum, nRuns; sigma=sigma);
serialize(outputPath * "scen$(scenNum)_$(nRuns)runs_mean_width$(sigma).jls", sol);
p = plotMonteCarlo(sol);
savefig(p, outputPath * "scen$(scenNum)_$(nRuns)runs_mean_width$(sigma).png");
