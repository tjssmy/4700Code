function [ Wind ] = NormalRandWind(nTraj, paras)
Wind = randn(1, nTraj) * paras(1);
end
