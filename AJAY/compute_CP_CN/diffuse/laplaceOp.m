function [ value ] = laplaceOp( Un, x_index, y_index, D  )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
    j = y_index, i = x_index;
    value = (Un(j,i+1) - 2 * Un(j,i) + Un(j, i-1))/D.hx^2 + (Un(j+1,i) - 2 * Un(j,i) + Un(j-1, i))/D.hy^2;
    
end

