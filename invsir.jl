using PyPlot
using RandomQuantum
using DelimitedFiles


dat = readdlm("COVID.dat")
d = dat[:,2]
z = dat[:,1]
N = 1260140

dt = 1
tf = length(z) 
t = 1:dt:tf
h = tf/length(t)
gamma = 0.08

function L2(d,s)
 return (sum((d .- s).^2) / sum(d.^2))
end


function direct(beta)

 xt = zeros(length(t))
 synt = zeros(length(t))
 zt = zeros(length(t))
 ynt = zeros(length(t))
 R0 = zeros(length(t))
 y = 1
 x = N - y
 
 function dxdt(x,y,b)
  return (-b*(x*y))
 end
 
 function dydt(x,y,b,g)
  return (b*(x*y) - g * y)
 end
 
 for i in 1:length(t)
  
  b = beta[i]
  g = gamma
  k1 = h * dxdt(x, y, b) 
  l1 = h * dydt(x, y, b, g)
  k2 = h * dxdt(x+k1/2, y+l1/2, b)
  l2 = h * dydt(x+k1/2, y+l1/2, b, g)
  k3 = h * dxdt(x+k2/2, y+l2/2, b)
  l3 = h * dydt(x+k2/2, y+l2/2, b, g)
  k4 = h * dxdt(x+k3, y+l3, b)
  l4 = h * dydt(x+k3, y+l3, b, g)
  
  x = x + (k1+2*k2+2*k3+k4)/6
  y = y + (l1+2*l2+2*l3+l4)/6

  synt[i] = y
 
 end
 return (synt)
end

bmax = 3.5
bmin = 0.02
beta = rand(length(t))
beta = ((bmax-bmin) .* beta .+ bmin)./N
ini = direct(beta)
e0 = L2(d, ini)
iter = 0

println(e0)
while true

 for i in 1:(length(beta)-1)

  global save = beta[i]

  global pert = ((2) * rand(1)[1] - 1)
  global beta[i] = beta[i] * N + pert
 
 if beta[i] > bmax || beta[i] < bmin 
   beta[i] = ((bmax-bmin) * rand(1)[1] + bmin)
  end
  
  beta[i] = beta[i]/N
  global syn = direct(beta)
  e = L2(d, syn)

  if e < e0 
   global e0 = e
  else
   beta[i] = save
  end
  
  #println(i)
  #plot(t,syn)
  #plot(t,d)
  #draw()
  #pause(0.0001)
  #clf()
  

 end

 global iter = iter + 1
 if iter == 1000 || e0 < 0.13
  break
 end


#println(iter)

end

println(e0)
println(iter)
plot(t,d)
plot(t,syn)
show()
