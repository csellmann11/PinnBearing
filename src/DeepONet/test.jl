using Lux

include("DeepONet.jl")

branch_nets = [[1,3,10,5],[1,3,10,5]]
trunc_net = [1,3,10,20]

y_data = rand(Float32,100)

net = DeepONet(branch_nets,trunc_net,y_data)

rng = Random.MersenneTwister(1234)

ps, st = Lux.setup(rng, net)


x = [rand(Float32,1,100),rand(Float32,1,1),rand(Float32,1,1)]

net(x,ps,st)

net.flags.training = false
@time net(x,ps,st)