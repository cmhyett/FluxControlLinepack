basepath = pwd() * "../";
using Pkg;
using GasTranSim, Serialization, Plots;
include(basepath * "./examples/runExample.jl");
include(basepath * "./src/plotting.jl");

scenNum, nRuns, mu, sigma = ARGS;
nRuns = parse(Int,nRuns);
scenNum = parse(Int, scenNum);
mu = parse(Float64, mu);
sigma = parse(Float64, sigma);

outputPath = basepath * "/results/";
mkpath(outputPath);
sol = monteCarlo(scenNum, nRuns; mu=mu, sigma=sigma);
serialize(outputPath * "scen$(scenNum)_$(nRuns)runs_mean$(mu)_width$(sigma).jls", sol);
p = plotMonteCarlo(sol);
savefig(p, outputPath * "scen$(scenNum)_$(nRuns)runs_mean$(mu)_width$(sigma).png");
