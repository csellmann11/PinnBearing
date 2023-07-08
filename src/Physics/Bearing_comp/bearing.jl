mutable struct HD_Bearing
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
end

function HD_Bearing_Pfeil(;d = 0.1,b = 0.1,psi = 3.58e-03,om = 100,eta = 0.01)

    c = 150e-06
    rI = d/2
    circ = 2pi*rI
    B = b / rI

    alpha = 1/B^2

    HD_Bearing(d,b,psi,rI,om,eta,circ,c,B,alpha)
end

function HD_Bearing_Schweizer(;d = 18e-03,b = 10e-03,psi = 3.58e-03,om = 100,eta = 6.4e-03)

    rI = d/2
    circ = 2pi*rI
    B = b / rI
    c = psi * rI
    alpha = 1/B^2

    HD_Bearing(d,b,psi,rI,om,eta,circ,c,B,alpha)
end
