
%% Uniform distribution
% function : Generate uniform distribution between 
%            lower and upper, key is the number of elements in x
% Example  : x = generateUniform(3,-2,2) -> x = [-2,0,2]
function x = generateUniform(key,lower,upper)
    % 如果只返回一个数，就取上界和下届之间的中值
    if key == 1
        x = lower + (upper-lower)/2;
    else
        d = (upper-lower)/(key-1);
        x = lower:d:upper;
    end
end

