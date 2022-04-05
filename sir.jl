using PyPlot


N = 10000000
x = N - 10
y = 10

beta = 0.5/N
gamma = 0.2
dt = 0.001
tf = 100
t = 0:dt:tf
h = tf/length(t)

xt = zeros(length(t))
yt = zeros(length(t))
zt = zeros(length(t))
ynt = zeros(length(t))
R0 = zeros(length(t))

function dxdt(x,y)
 return (-beta*(x*y))
end
 
function dydt(x,y)
 return (beta*(x*y) - gamma * y)
end

for i in 1:length(t)

 k1 = h * dxdt(x, y) 
 l1 = h * dydt(x, y)
 k2 = h * dxdt(x+k1/2, y+l1/2)
 l2 = h * dydt(x+k1/2, y+l1/2)
 k3 = h * dxdt(x+k2/2, y+l2/2)
 l3 = h * dydt(x+k2/2, y+l2/2)
 k4 = h * dxdt(x+k3, y+l3)
 l4 = h * dydt(x+k3, y+l3)
 
 global x = x + (k1+2*k2+2*k3+k4)/6
 global y = y + (l1+2*l2+2*l3+l4)/6
 global yn = beta * x * y
 global r =  beta * x * N 

 xt[i] = x 
 ynt[i] = yn
 R0[i] = r
 yt[i] = y
 zt[i] = N - x - y

end

plot(t,xt)
plot(t,yt)
plot(t,zt)
plot(t,ynt, color = "red", ls="dashed" )
plot(t,R0, color = "black", ls="dashed" )
