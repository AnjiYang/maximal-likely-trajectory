function [f] = uexact(x,y,t)
f = normpdf(x,4.1,1/10)*normpdf(y,4.1,1/10); 
end
