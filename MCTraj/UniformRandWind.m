function [ Wind ] = UniformRandWind(nTraj,paras)
Wind = (rand(1,nTraj)-0.5)*paras(1);
end

