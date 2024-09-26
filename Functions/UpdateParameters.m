function [model] = UpdateParameters(model,parameters2Update,StartIndex)

%Updating Model parameters starting from an index. 
lengthOfParameters = length(parameters2Update);
for i = 1:lengthOfParameters
    index2Update = i+StartIndex-1;
    model.Parameters(index2Update).value = parameters2Update(i);
end 

end 