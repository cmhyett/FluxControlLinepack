using CSV;
include("../types/networkTypes.jl");

function parsePipes(inputFile::String)::Array{Pipe,1}
    file = CSV.File(inputFile);
    numPipes = length(file);
    output = [];
    for i in 1:numPipes
        push!(output, Pipe(file[i].id,
                           file[i].startNode,
                           file[i].endNode,
                           file[i].diameter,
                           file[i].length,
                           file[i].frictionFactor,
                           zeros(1,1),
                           zeros(1,1),
                           0,
                           0));
    end

    return output;
end

function parseNodes(staticInputFile::String,
                    dynamicInputFile::String,
                    pipes::Array{Pipe, 1})::Array{Junction,1}
    structureFile = CSV.File(staticInputFile);
    dynamicFile = CSV.File(dynamicInputFile, header=false);
    numTsteps = length(dynamicFile)-1;
    numJunctions = length(structureFile);
    output = [];
    
    for i in 1:numJunctions
        nodeId = structureFile.nodeId[i];
        nodeType = (structureFile.nodeType[i] == 1 ? SlackNode : WithdrawalNode);
        col = findfirst(x->x==nodeId, dynamicFile[1]);

        incidentPipes = [];
        for i in 1:length(pipes)
            if ((pipes[i].fromJunction == nodeId) ||
                (pipes[i].toJunction == nodeId))
                push!(incidentPipes, pipes[i].id);
            end
        end
        incidentPipes = convert(Array{Int, 1}, incidentPipes);
        
        if (structureFile.nodeType[i] == 1) #slackNode
            push!(output, Junction(nodeId,
                                   SlackNode,
                                   dynamicFile[col][2:end],
                                   zeros(numTsteps),
                                   [],
                                   incidentPipes));
        else # withdrawalNode
            push!(output, Junction(nodeId,
                                   WithdrawalNode,
                                   zeros(numTsteps),
                                   dynamicFile[col][2:end],
                                   [],
                                   incidentPipes));
        end            
    end
    return convert(Array{Junction, 1}, output);
end

function parseNetworkFromFiles(nwFilepath::String,
                                nodeFilepath::String,
                                pipeFilepath::String,
                                timeseriesFilepath::String)::Network
    nwFile = CSV.File(nwFilepath);
    pipes = parsePipes(pipeFilepath);
    nodes = parseNodes(nodeFilepath,
                       timeseriesFilepath,
                       pipes);

    @assert nwFile.units[1] == "si"; #hack
    
    return Network(nwFile.gasSpecificGravity[1],
                   nwFile.specificHeatCapacityRatio[1],
                   nwFile.temperature[1],
                   nwFile.R[1],
                   nwFile.compressibilityFactor[1],
                   nwFile.basePressure[1],
                   nwFile.baseLength[1],
                   nwFile.baseFlow[1],
                   si,
                   nwFile.soundSpeed[1],
                   pipes,
                   nodes,
                   nwFile.dt[1]);
end
