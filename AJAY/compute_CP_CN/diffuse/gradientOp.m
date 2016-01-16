function [ value ] = gradientOp( Un, x_index, y_index, D  )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
    j = y_index, i = x_index;
    value(1) = (Un(j, i+1) - Un(j, i-1))/2/D.hx;
    value(2) = (Un(j+1,i) - Un(j-1,i))/2/D.hy;
end

