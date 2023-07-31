abstract type Bearing end


mutable struct HD_Bearing <: Bearing
    d       :: Float32
    b       :: Float32
    psi     :: Float32
    rI      :: Float32
    om      :: Real    # om is either Float or Dual number
    eta     :: Float32
    circ    :: Float32

    c       :: Float32
    B       :: Float32
    alpha   :: Float32

    sysMat_dec :: Any
end

function HD_Bearing_Pfeil(;d = 0.1,b = 0.1,psi = 3.58e-03,om = 100,eta = 0.01)

    c = 150e-06
    rI = d/2
    circ = 2pi*rI
    B = b / rI

    alpha = 1/B^2

    sysMat_dec = nothing

    HD_Bearing(d,b,psi,rI,om,eta,circ,c,B,alpha,sysMat_dec)
end

function HD_Bearing_Schweizer(;d = 18e-03,b = 10e-03,psi = 3.58e-03,om = 100,eta = 6.4e-03)

    rI = d/2
    circ = 2pi*rI
    B = b / rI
    c = psi * rI
    alpha = 1/B^2

    sysMat_dec = nothing

    HD_Bearing(d,b,psi,rI,om,eta,circ,c,B,alpha,sysMat_dec)
end


mutable struct Foil_Bearing <: Bearing
    d       :: Float32
    b       :: Float32
    psi     :: Float32
    rI      :: Float32
    om      :: Float32
    eta     :: Float32
    circ    :: Float32

    c       :: Float32
    B       :: Float32
    alpha   :: Float32

    omega_max :: Float32

    Ksc     :: Float32
    p_atm   :: Float32
end

function Foil_Bearing_Leister(;d = 38.1e-03,b = 38.1e-03,psi = 3.58e-03,om = 100,eta = 1.85e-05,
    Eks = 214e09,ν = 0.3,omega_max = 400 * 2pi,t = 102e-06, lo = 1.778e-03, s0 = 4.572e-03 ,p_atm = 1.01325e05)

    rI = d/2
    circ = 2pi*rI
    B = b / rI

    c = 50e-06
    alpha = 1/B^2

    Ksc = Eks * t^3/(2 * s0 *  (1- ν^2) * lo^3) * c

    Foil_Bearing(d,b,psi,rI,om,eta,circ,c,B,alpha,omega_max,Ksc,p_atm)
    
end