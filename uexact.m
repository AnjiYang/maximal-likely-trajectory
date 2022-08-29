% 精确解函数-初始函数
function [f] = uexact(x,y,t)
f = normpdf(x,4.1,1/10)*normpdf(y,4.1,1/10); % 初始值为(1.5 1.5)，标准差为1/2,初始值反着放入，为什么？哪里有问题？
end