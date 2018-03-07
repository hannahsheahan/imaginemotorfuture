function y=resample(x,frames)
    y = interp1(1:length(x),x,linspace(1,length(x),frames));
end
