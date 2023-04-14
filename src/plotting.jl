using Plots;
using Statistics;
using Printf;

function generateExamplePlot()
    T = 3600.0;
    t = 0.0:1:T;
    numTsteps = length(t);

    pressure = [80 for i in t];
    perturbedPressure = zeros(numTsteps);
    for i in 1:numTsteps
        if (t[i] < 2400)
            perturbedPressure[i] = pressure[i];
        else
            perturbedPressure[i] = 1.1*pressure[i];
        end
    end

    plt = plot(t, perturbedPressure, label=:none, linecolor=2, linestyle=:dash, legend=:topleft)
    plot!(plt, t, pressure, xlabel = "seconds", ylabel="Pressure (bar)", label="nominal pressure", color=1);
    plot!(plt, [], [], color=2, linestyle=:dash, label="perturbed pressure")
    plot!(plt, ylims=(70,90))

    withdrawal = 31*(0.1 .* sin.(2π .* (t./T)) .+ 1);
    perturbedWithdrawal = zeros(numTsteps);
    for i in 1:numTsteps
        if (t[i] < 3000)
            perturbedWithdrawal[i] = withdrawal[i];
        else
            perturbedWithdrawal[i] = 1.3*withdrawal[i];
        end
    end
    rightYlims = (0.8*minimum(perturbedWithdrawal), 1.2*maximum(perturbedWithdrawal));
    plot!(twinx(plt), t, perturbedWithdrawal, color=4, linestyle=:dash, ylims=rightYlims, xticks=:none, yticks=:none, legend=:none);
    plot!(plt, right_margin=15Plots.mm)
    plot!(twinx(plt), t, withdrawal, color=3, legend=:none, ylabel="Withdrawal (kg/m/s)", ylims=rightYlims, xticks=:none)
    plot!(plt, [], [], color=3, label="nominal withdrawal");
    plot!(plt, [], [], color=4, linestyle=:dash, label="perturbed withdrawal");
    

    #plot annotations
    xstart = 2350;
    xend = xstart;
    ystart = 80;
    yend = perturbedPressure[end];
    plot!(plt, [xstart, xend], [ystart, yend], arrow=true, label=:none, color=:black);
    plot!(plt, [xend, xstart], [yend, ystart], arrow=true, label=:none, color=:black);
    annotate!(plt, (2250, ystart+0.5*(yend-ystart), "Δp", 10));

    xstart = 2400;
    xend = 3000;
    ystart = 80.5;
    yend = 80.5;
    plot!(plt, [xstart, xend], [ystart, yend], arrow=true, label=:none, color=:black);
    plot!(plt, [xend, xstart], [yend, ystart], arrow=true, label=:none, color=:black);
    annotate!(plt, (xstart+0.5*(xend-xstart), 81, "Δt", 10));
    
    return plt;
end

function calculateLinepack(ts)

    numTsteps = length(ts.sol["time_points"]);
    numPipes = length(ts.sol["pipes"]);
    integratedLinepack = zeros(numTsteps);
    linepack = zeros(numTsteps);
    pipesSol = [ts.sol["pipes"]["$(i)"] for i in 1:numPipes];
    for p in pipesSol
        linepack += p["in_flow"] - p["out_flow"];
    end
    return integrateArray(linepack, ts.sol["time_step"]);
end

function integrateArray(f, dx)
    sums = zeros(length(f)-1);
    sums[1] = (f[1]+f[2])*dx/2;
    for i in 2:length(f)-1
        sums[i] = sums[i-1] + (f[i]+f[i+1])*dx/2;
    end
    return sums;
end

function plotPressureAtNodes(ts)
    t = ts.sol["time_points"]/3600;
    plot(t, ts.sol["nodes"]["1"]["pressure"])
    plot!(t, ts.sol["nodes"]["2"]["pressure"])
    plot!(t, ts.sol["nodes"]["3"]["pressure"])
    plot!(t, ts.sol["nodes"]["4"]["pressure"])
    plot!(xlabel = "time (hours)", ylabel = "Pressure (Pa)")
end

function plotPressureBoundaryCondtions(dataPath; nodeNum=1)
    simulations = deserialize(dataPath);
    p = plot();
    for i in 1:length(simulations)
        ts = simulations[i].ts;
        t = ts.sol["time_points"];
        pr = ts.sol["nodes"]["$(nodeNum)"]["pressure"];
        plot!(p, t, pr);
    end
    plot!(p, legend=false, ylabel="Pressure (Pa)", xlabel="Time (s)", title="Pressure at node $(nodeNum)");
    return p;
end


function linepackPlot(ts; currentTime_hrs = nothing)
    t = [ts.sol["pipes"]["1"]["density_profile"][i].t for i in 1:length(ts.sol["pipes"]["1"]["density_profile"])] ./ 3600;
    linepack = zeros(length(t));
    northLinepack = zeros(length(t));
    southLinepack = zeros(length(t));
    
    southKeys = ["2", "3", "4", "5"];

    pressureLimits = Dict("1"=>58*1e5,
                          "2"=>53*1e5,
                          "3"=>53*1e5,
                          "4"=>53*1e5,
                          "5"=>51*1e5,
                          "6"=>50*1e5,
                          "7"=>55*1e5,
                          "8"=>55*1e5,
                          "9"=>50*1e5,
                          "10"=>45*1e5,
                          "11"=>45*1e5);
    pressureCrossings = Dict();
    
    for p in keys(ts.sol["pipes"])
        pipe = ts.sol["pipes"][p];
        avgPressure = [mean(pipe["density_profile"][i].val) for i in 1:length(pipe["density_profile"])];
        avgDensity = GasTranSim._pressure_to_density_full_cnga(avgPressure ./ ts.nominal_values[:pressure],
                                                      ts.nominal_values,
                                                      ts.params) .* ts.nominal_values[:density];
        len = ts.data["pipes"][p]["length"] * ts.nominal_values[:length];
        area = ts.data["pipes"][p]["area"] * ts.nominal_values[:area];
        linepack .+= avgDensity .* (len * area)
        if (p in southKeys)
            southLinepack .+= avgDensity .* (len * area);
        else
            northLinepack .+= avgDensity .* (len * area);
        end
    end
    numDaysx2 = floor(Int, t[end]/6);
    p = plot(xlabel="time (hrs)", ylabel = "Linepack (kg)", legend=:bottomleft, xticks=([6*i for i in 1:numDaysx2]));
    plot!(p, t, linepack, label="Total Linepack");
    plot!(p, t, northLinepack, label="North Linepack");
    plot!(p, t, southLinepack, label="South Linepack");

    if (currentTime_hrs != nothing)
        scatter!(p, [currentTime_hrs], [ylims(p)[1]*1.1], markershape=:utriangle, label=false, ylims=ylims(p));
    end

    # plot pressure at indicative nodes
    pressure1 = ts.sol["nodes"]["1"]["pressure"] ./ 1e5;
    pressure8 = ts.sol["nodes"]["8"]["pressure"] ./ 1e5;
    pressure6 = ts.sol["nodes"]["6"]["pressure"] ./ 1e5;

    rightYlims = (0.5*minimum([pressure1; pressure8; pressure6]), 1.5*maximum([pressure1; pressure8; pressure6]))

    plot!(twinx(p), t, pressure1, color=4, ylims=rightYlims, xticks=:none,  legend=:none, ylabel="Pressure (bar)", linestyle=:dash);
    plot!(p, [], [], color=4, label="Node 1 Pressure", linestyle=:dash);

    plot!(twinx(p), t, pressure8, color=5, ylims=rightYlims, xticks=:none,  legend=:none, linestyle=:dash, yticks=:none);
    plot!(p, [], [], color=5, label="Node 8 Pressure", linestyle=:dash);

    plot!(twinx(p), t, pressure6, color=6, ylims=rightYlims, xticks=:none,  legend=:none, linestyle=:dash, yticks=:none);
    plot!(p, [], [], color=6, label="Node 6 Pressure", linestyle=:dash);

    plot!(p, right_margin=5Plots.mm, size=(800,600))

    # plot pressure crossings
    startY, endY = ylims(p);
    for n in keys(ts.sol["nodes"])
        node = ts.sol["nodes"][n];
        if (n in keys(pressureLimits))
            firstCrossingTimeInd = findfirst(x->x<=pressureLimits[n], node["pressure"]);
            if (firstCrossingTimeInd != nothing)
                firstCrossTime = t[firstCrossingTimeInd];
                plot!(p, [firstCrossTime; firstCrossTime], [startY; endY-0.75*endY], label=false, color=:red, linestyle=:dashdot)
                plot!(p, annotations = ([firstCrossTime], [endY * (0.78 + 0.02*parse(Int, n))], Plots.text.([n], 10, :red, :center)));
            end
        end
    end

    return p;
end

function plotNetworkDiagram(ts)
    x = [];
    y = [];
    nodeNums = [];
    for k in keys(ts.data["nodes"])
        if (k != "12" && k != "13") #12 & 13 are fictitious
            push!(x, ts.data["nodes"][k]["x_coord"]);
            push!(y, ts.data["nodes"][k]["y_coord"]);
            push!(nodeNums, k)
        end
    end
    Δx = maximum(x) - minimum(x)
    Δy = maximum(y) - minimum(y)
    Δs = max(Δx,Δy)

    p = plot(legend = false, aspect_ratio=:equal, showaxis=:hide, xaxis=nothing,
             yaxis=nothing, grid=false)
    for k = keys(ts.data["pipes"])
        if (parse(Int,k) <=11)
            i = parse(Int64,k)
            from_node = ts.data["pipes"][k]["from_node"]
            to_node = ts.data["pipes"][k]["to_node"]
            from_x, to_x = ts.data["nodes"]["$from_node"]["x_coord"], ts.data["nodes"]["$to_node"]["x_coord"]
            from_y, to_y = ts.data["nodes"]["$from_node"]["y_coord"], ts.data["nodes"]["$to_node"]["y_coord"]
            N = ts.ref[:pipe][i]["num_discretization_points"]
            
            width = 7.0;
            color = RGB(55/255, 134/255, 230/255);
            plot!(p, [from_x; to_x], [from_y; to_y], linewidth = width,
                  color=color, label="")
        end
        if (k == "1" || k == "2") #double pipes
            width = 14.0;
            plot!(p, [from_x; to_x], [from_y; to_y], linewidth = width,
                  color=color, label="")
            plot!(p, [from_x; to_x], [from_y; to_y], linewidth = 3.0,
                  color=:white, label="")
        end
    end

    scatter!(p, x,y,
             xlim = (minimum(x)-0.05*Δs, maximum(x)+0.05*Δs),
             ylim = (minimum(y)-0.05*Δs, maximum(y)+0.05*Δs),
             markersize = 14.0,
             markeralpha= 1.0,
             markercolor= :lightgreen,
             markerstrokewidth = 0.0)
    plot!(p, annotations = (x, y, Plots.text.(nodeNums, :center)));
    return p;
end

function plot_pressure_profile(filename, ts; cm=cgrad(:roma, 10, categorical=true, scale=:exp), fps=15)
    x = [];
    y = [];
    nodeNums = [];
    for k in keys(ts.data["nodes"])
        if (k != "12" && k != "13") #12 & 13 are fictitious
            push!(x, ts.data["nodes"][k]["x_coord"]);
            push!(y, ts.data["nodes"][k]["y_coord"]);
            push!(nodeNums, k)
        end
    end
    t = [ts.sol["pipes"]["1"]["density_profile"][i].t for i in 1:length(ts.sol["pipes"]["1"]["density_profile"])];
    Δx = maximum(x) - minimum(x)
    Δy = maximum(y) - minimum(y)
    Δs = max(Δx,Δy)
    
    vmax = 0
    vmin = Inf
    for i in 1:length(t)
        for k = keys(ts.data["pipes"])
            dp = ts.sol["pipes"][k]["density_profile"];
            vmin = min(vmin, minimum(dp[i].val)) / 1e5;
            vmax = max(vmax, maximum(dp[i].val)) / 1e5;
        end
    end
    vmax = 80;
    vmin = 40;
    anim = @animate for num in 1:length(t)
        plot(legend = false, aspect_ratio=:equal, showaxis=:hide, xaxis=nothing,
            yaxis=nothing, grid=false)
        for k = keys(ts.data["pipes"])
            i = parse(Int64,k)
            from_node = ts.data["pipes"][k]["from_node"]
            to_node = ts.data["pipes"][k]["to_node"]
            from_x, to_x = ts.data["nodes"]["$from_node"]["x_coord"], ts.data["nodes"]["$to_node"]["x_coord"]
            from_y, to_y = ts.data["nodes"]["$from_node"]["y_coord"], ts.data["nodes"]["$to_node"]["y_coord"]
            N = ts.ref[:pipe][i]["num_discretization_points"]
            rho = ts.sol["pipes"][k]["density_profile"][num].val ./ 1e5;
            
            for n = 1:N
                dx = to_x - from_x
                dy = to_y - from_y
                x1 = dx * (n-1) / N + from_x
                x2 = dx * n / N + from_x
                y1 = dy * (n-1) / N + from_y
                y2 = dy * n / N + from_y
                plot!([x1; x2], [y1; y2], linewidth = 7.0,
                      color=cgrad(cm)[max(min((rho[n]-vmin)/(vmax-vmin), 1),0)] ,label="")
            end
        end

        p1 = scatter!(x,y,
                      xlim = (minimum(x)-0.03*Δs, maximum(x)+0.03*Δs),
                      ylim = (minimum(y)-0.03*Δs, maximum(y)+0.03*Δs),
                      markersize = 14.0,
                      markeralpha= 1.0,
                      markercolor= :lightgreen,
                      markerstrokewidth = 0.0)
        formattedTime = "time (hrs) = $(@sprintf("%.3f", t[num]/3600))";
        #p1 = plot!(annotations = (xlims(p1)[1], ylims(p1)[2], Plots.text(formattedTime, :left)))
        plot!(p1, xlabel=formattedTime);
        plot!(p1, annotations = (x, y, Plots.text.(nodeNums, :center)));
        h2 = scatter([0,0], [0,1], zcolor=[0,1], clims=(vmin,vmax),
                     xlims=(1,1.1), label="", c=cm, colorbar_title="Pressure (bar)", framestyle=:none, margin=0.0Plots.mm)
        p2 = linepackPlot(ts, currentTime_hrs = t[num]/3600);
        l = @layout [a{0.01w} a{0.4w} a{0.6w}]
        p = plot(h2, p1, p2, layout=l, link=:all, size=(1200,800), margin=10Plots.mm, right_margin=20Plots.mm);
        if (num == floor(Int, length(t)/2))
            savefig(p, filename*".png");
        end
        p;
    end
    gif(anim, filename, fps = fps)
end

# function plotMonteCarlo(runs; nodes=[1,6,8], quantiles=[0.125, 0.375, 0.5, 0.625, 0.875],
#                         plotLinepack=true) #runs = Dict(Int, transientsim)
#     @assert isodd(length(quantiles));
#     p = plot();
#     ks = Array([keys(runs)...]);
#     t_hrs = runs[ks[1]].sol["time_points"] ./ 3600;

#     linepack = zeros(length(ks), length(t_hrs))
#     for i in 1:length(ks)
#         local ts = runs[ks[i]];
#         for n in keys(ts.sol["pipes"])
#             pipe = ts.sol["pipes"][n];
#             avgPressure = [mean(pipe["density_profile"][i].val) for i in 1:length(pipe["density_profile"])];
#             avgDensity = GasTranSim._pressure_to_density_full_cnga(avgPressure ./ ts.nominal_values[:pressure],
#                                                                    ts.nominal_values,
#                                                                    ts.params) .* ts.nominal_values[:density];
#             len = ts.data["pipes"][n]["length"] * ts.nominal_values[:length];
#             area = ts.data["pipes"][n]["area"] * ts.nominal_values[:area];
#             linepack[i,:] .+= avgDensity .* (len * area)
#         end
#     end
    
#     for node in nodes
#         pressures = zeros(length(ks), length(t_hrs));
#         for i in 1:length(ks)
#             pressures[i,:] = runs[ks[i]].sol["nodes"]["$(node)"]["pressure"];
#         end
#         qs = [quantile(pressures[:,i], quantiles) for i in 1:size(pressures)[2]];
#         qs = hcat(qs...);
#         for i in 1:floor(Int, length(quantiles)/2)
#             plot!(p, t_hrs, qs[i,:], fillrange=qs[length(quantiles)-i+1,:], fillalpha=quantiles[i], c=node, 
#                   label=nothing, linealpha=0);
#         end
#         plot!(p, t_hrs, qs[floor(Int, length(quantiles)/2)+1,:], c=node, label="\$ P_{n$(node)}(t)\$", legend=:bottomleft);
#     end
#     hrsPerTick = 4
#     numTicks = floor(Int, t_hrs[end]/hrsPerTick);
#     plot!(p, xlabel="time (hrs)", ylabel = "Pressure (pa)", legend=:bottomleft, xticks=([hrsPerTick*i for i in 1:numTicks]));

#     addedPressureCrossings = false;
# #        addPressureCrossings(p, runs, nodes=[1,6,9], quantiles=[0.125, 0.375, 0.5, 0.625, 0.875]);

#     if (plotLinepack)
#         qs = [quantile(linepack[:,i], quantiles) for i in 1:size(linepack)[2]];
#         qs = hcat(qs...);
#         ylims = (0.8*minimum(linepack), 1.2*maximum(linepack));
#         for i in 1:floor(Int, length(quantiles)/2)
#             plot!(twinx(p), t_hrs, qs[i,:], fillrange=qs[length(quantiles)-i+1,:], fillalpha=quantiles[i], c=2, 
#                   label=nothing, linealpha=0, xticks=:none, yticks=:none, ylims=ylims);
#         end
#         plot!(twinx(p), t_hrs, qs[floor(Int, length(quantiles)/2)+1,:], c=2, label="Total Linepack", legend=false,
#               xticks=:none, ylims=ylims, linestyle=:dash, ylabel="Linepack (kg)");
#         plot!(p, right_margin=25Plots.mm);
#         plot!(p, [],[],c=2,linestyle=:dash,label="Linepack");
#     end
#     return p;

# end

function plotMonteCarlo(runs; nodes=[1,6,8], quantiles=[0.125, 0.375, 0.5, 0.625, 0.875],
                        plotLinepack=true, plims=(55,85), llims=(7.0e5,8.5e5),
                        addPCrossings=false, pCrossNodes=[6,9,1]) #runs = Dict(Int, transientsim)
    p = plotMonteCarloPressure(runs, nodes=nodes, quantiles=quantiles, plims=plims);
    if (addPCrossings)
        addPressureCrossings(p, runs; nodes=pCrossNodes)
    end
    if (plotLinepack)
        plotMonteCarloLinepack(runs, quantiles=quantiles, p=p, llims=llims);
    end
    return p;
end

function plotMonteCarloPressure(runs; nodes=[1,8,6], quantiles=[0.125, 0.375, 0.5, 0.625, 0.875],plims=(55,85)) #runs = Dict(Int, transientsim)
    pascalToBar = 1e-5;
    @assert isodd(length(quantiles));
    p = plot();
    ks = Array([keys(runs)...]);
    t_hrs = runs[ks[1]].sol["time_points"] ./ 3600;
    for node in nodes
        pressures = zeros(length(ks), length(t_hrs));
        for i in 1:length(ks)
            pressures[i,:] = runs[ks[i]].sol["nodes"]["$(node)"]["pressure"] .* pascalToBar;
        end
        qs = [quantile(pressures[:,i], quantiles) for i in 1:size(pressures)[2]];
        qs = hcat(qs...);
        for i in 1:floor(Int, length(quantiles)/2)
            plot!(p, t_hrs, qs[i,:], fillrange=qs[length(quantiles)-i+1,:], fillalpha=quantiles[i], c=node, 
                  label=nothing, linealpha=0);
        end
        plot!(p, t_hrs, qs[floor(Int, length(quantiles)/2)+1,:], c=node, label="\$ P_{$(node)}\$", legend=:bottomleft);
    end
    hrsPerTick = 4
    numTicks = floor(Int, t_hrs[end]/hrsPerTick);
    plot!(p, xlabel="time (hrs)", ylabel = "Pressure (bar)", legend=:bottomleft,
          xticks=([hrsPerTick*i for i in 1:numTicks]),
          ylims=plims);
    return p;
end

function plotMonteCarloLinepack(runs; quantiles=[0.125, 0.375, 0.5, 0.625, 0.875],p=nothing, llims=(7.0e5,8.5e5)) #runs = Dict(Int, transientsim)
    kgToMMBTU = 0.052;
    @assert isodd(length(quantiles));
    if (p == nothing)
        p = plot();
    end
    ks = Array([keys(runs)...]);
    t_hrs = runs[ks[1]].sol["time_points"] ./ 3600;
    linepack = zeros(length(ks), length(t_hrs))
    northLinepack = zeros(length(ks), length(t_hrs))
    southLinepack = zeros(length(ks), length(t_hrs))
    southKeys = ["2", "3", "4", "5"];
    for i in 1:length(ks)
        local ts = runs[ks[i]];
        for p in keys(ts.sol["pipes"])
            pipe = ts.sol["pipes"][p];
            avgPressure = [mean(pipe["density_profile"][i].val) for i in 1:length(pipe["density_profile"])];
            avgDensity = GasTranSim._pressure_to_density_full_cnga(avgPressure ./ ts.nominal_values[:pressure],
                                                                   ts.nominal_values,
                                                                   ts.params) .* ts.nominal_values[:density];
            len = ts.data["pipes"][p]["length"] * ts.nominal_values[:length];
            area = ts.data["pipes"][p]["area"] * ts.nominal_values[:area];
            linepack[i,:] .+= avgDensity .* (len * area)
            if (p in southKeys)
                southLinepack[i,:] .+= avgDensity .* (len * area);
            else
                northLinepack[i,:] .+= avgDensity .* (len * area);
            end
        end
    end
    linepack .*= kgToMMBTU; #convert units
    
    qs = [quantile(linepack[:,i], quantiles) for i in 1:size(linepack)[2]];
    qs = hcat(qs...);
    for i in 1:floor(Int, length(quantiles)/2)
        plot!(twinx(p), t_hrs, qs[i,:], fillrange=qs[length(quantiles)-i+1,:], fillalpha=quantiles[i]*1.5, c=2, 
              label=nothing, linealpha=0, xticks=:none, yticks=:none, ylims=llims);
        # plot!(p, t_hrs, qs[i,:], fillrange=qs[length(quantiles)-i+1,:], fillalpha=quantiles[i], c=2, 
        #       label=nothing, linealpha=0);
    end
    plot!(twinx(p), t_hrs, qs[floor(Int, length(quantiles)/2)+1,:], c=2, label="Total Linepack", legend=false,
          xticks=:none, ylims=llims, linestyle=:dash, ylabel="Linepack (MMBTU)");
    plot!(p, right_margin=20Plots.mm);
    plot!(p, [],[],c=2,linestyle=:dash,label="Linepack");
    
    return p;
end

function addPressureCrossings(p, runs; nodes=[1,6,11], quantiles=[0.125, 0.375, 0.5, 0.625, 0.875])
    @assert isodd(length(quantiles));
    pressureLimit = 50e5;
    ks = Array([keys(runs)...]);
    t_hrs = runs[ks[1]].sol["time_points"] ./ 3600;
    addedSomething = false;
    for node in nodes
        timeOfFirstCrossing = [];
        for i in 1:length(ks)
            firstCrossingTimeInd = findfirst(x->x<=pressureLimit, runs[ks[i]].sol["nodes"]["$(node)"]["pressure"]);
            if (firstCrossingTimeInd != nothing)
                push!(timeOfFirstCrossing, t_hrs[firstCrossingTimeInd]);
            else
                push!(timeOfFirstCrossing, Inf);
            end
        end
        if (length(timeOfFirstCrossing)>0)
            addedSomething = true;
            qs = quantile(timeOfFirstCrossing, quantiles);
            startY, endY = ylims(p);
            for i in 1:floor(Int, length(quantiles)/2)
                plot!(p, [qs[i], qs[length(quantiles)-i+1]], [startY, startY],
                      fillrange=0.8*[endY, endY], linealpha=0, fillalpha=quantiles[i], color=:red,
                      label=nothing);
            end
            medianXVal = qs[floor(Int, length(quantiles)/2)+1];
            plot!(p, [medianXVal, medianXVal], [startY, 0.8*endY],color=:red,label=nothing);
            plot!(p, annotations = ([medianXVal], [0.85*endY], Plots.text.(["\$ $(node) \$"], 12, :red, :center)));
            plot!(p, ylims=(startY,endY));
        end
    end
    return addedSomething;
end

function plotHistogramOfCrossingTimes(runs; node=6, nbins=20)
    pressureLimit = 50e5;
    ks = Array([keys(runs)...]);
    t_hrs = runs[ks[1]].sol["time_points"] ./ 3600;
    timeOfFirstCrossing = [];
    for i in 1:length(ks)
        firstCrossingTimeInd = findfirst(x->x<=pressureLimit, runs[ks[i]].sol["nodes"]["$(node)"]["pressure"]);
        if (firstCrossingTimeInd != nothing)
            push!(timeOfFirstCrossing, t_hrs[firstCrossingTimeInd]);
        else
            push!(timeOfFirstCrossing, Inf);
        end
    end
    hrsPerTick = 4
    numTicks = floor(Int, t_hrs[end]/hrsPerTick);

    h = histogram(timeOfFirstCrossing, nbins=nbins, mode=:pdf, xlims=(t_hrs[1],t_hrs[end]), legend=false,
                  xticks=([hrsPerTick*i for i in 1:numTicks]),
                  xlabel="time (hrs)", ylabel="# Pressure crossings");
    return h, timeOfFirstCrossing;
end
function plotScen1(ts)
    d = Dict(1=>ts);
    p = plotMonteCarlo(d, plotLinepack=true, plims=(55,83),llims=(6.5e5, 9.1e5))
    return p;
end

function plotScen2(runs)
    p = plotMonteCarlo(runs, plotLinepack=true, plims=(55,83),llims=(6.5e5, 9.1e5))
    return p;
end

function plotScen3(runs)
    p = plotMonteCarlo(runs, addPCrossings=true, pCrossNodes=[6,9,1],
                       plotLinepack=true, plims=(35,83), llims=(6.0e5,10.5e5));
    return p;
end

function plotScen4(runs)
    p = plotMonteCarlo(runs, addPCrossings=true, pCrossNodes=[6,9,1],
                       plotLinepack=true, plims=(35,83), llims=(6.0e5,10.5e5));
    return p;
end

function plotScen5(runs)
    p = plotMonteCarlo(runs, addPCrossings=true, pCrossNodes=[6,9,1],
                       plotLinepack=true, plims=(35,83), llims=(6.0e5,10.5e5));
    return p;
end

function plotScen6(runs)
    p = plotMonteCarlo(runs, addPCrossings=true, pCrossNodes=[6,9,1],
                       plotLinepack=true, plims=(35,83), llims=(6.0e5,10.5e5));
    return p;
end
