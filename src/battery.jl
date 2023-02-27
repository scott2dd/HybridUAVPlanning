
function get_OCV_func()
    SOC=100 .*[1,0.9503,0.9007,0.8510,0.8013,0.7517,0.7020,0.6524,0.6027,0.5530,0.5034,0.4537,0.4040,0.3543,0.3046,0.2550,0.2053,0.1556,0.1059,0.0563];
    OCV=[4.1617,4.0913,4.0749,4.0606,4.0153,3.9592,3.9164,3.8687,3.8163,3.7735,3.7317,3.6892,3.6396,3.5677,3.5208,3.4712,3.3860,3.2880,3.2037,3.0747];

    OCV_interp = LinearInterpolation(reverse(SOC), reverse(OCV))
    return OCV_interp, 0, 0
end
function get_OCV_func_LiPo()
    SOC=100 .*  [1.0,
    0.9668412021343039,
    0.934531323117733,
    0.9045785889184145,
    0.8724856321908501,
    0.8410823292377329,
    0.7019553954577109,
    0.5677038583756852,
    0.4366584105242912,
    0.30752419166017014,
    0.2511283417275686,
    0.19575747375201058,
    0.14133157690308176]
    
    OCV= [ 16.788396797180177,
    16.62169780731201,
    16.48680353164673,
    16.392482948303222,
    16.33413486480713,
    16.282863426208497,
    16.039941120147706,
    15.684522533416748,
    15.49577431678772,
    15.29106364250183,
    15.126658964157105,
    14.949568367004394,
    14.800385189056396]

    OCV_interp = LinearInterpolation(reverse(SOC), reverse(OCV))
    return OCV_interp, minimum(SOC), maximum(SOC)
end
function get_OCV_func(SOC, OCV)
    # SOC=100 .*[1,0.9503,0.9007,0.8510,0.8013,0.7517,0.7020,0.6524,0.6027,0.5530,0.5034,0.4537,0.4040,0.3543,0.3046,0.2550,0.2053,0.1556,0.1059,0.0563];
    # OCV=[4.1617,4.0913,4.0749,4.0606,4.0153,3.9592,3.9164,3.8687,3.8163,3.7735,3.7317,3.6892,3.6396,3.5677,3.5208,3.4712,3.3860,3.2880,3.2037,3.0747];

    OCV_interp = LinearInterpolation(reverse(SOC), reverse(OCV))
    return OCV_interp
end
function get_OCV_table()
    SOC=100 .*[1,0.9503,0.9007,0.8510,0.8013,0.7517,0.7020,0.6524,0.6027,0.5530,0.5034,0.4537,0.4040,0.3543,0.3046,0.2550,0.2053,0.1556,0.1059,0.0563];
    OCV=[4.1617,4.0913,4.0749,4.0606,4.0153,3.9592,3.9164,3.8687,3.8163,3.7735,3.7317,3.6892,3.6396,3.5677,3.5208,3.4712,3.3860,3.2880,3.2037,3.0747];


    for i in 1:length(SOC)
        println("$(round(SOC[i],digits=4)) & $(round(OCV[i],digits=4))  \\\\ \\hline"   )
    end
end
function get_one_by_V()
    SOC=100 .*[1,0.9503,0.9007,0.8510,0.8013,0.7517,0.7020,0.6524,0.6027,0.5530,0.5034,0.4537,0.4040,0.3543,0.3046,0.2550,0.2053,0.1556,0.1059,0.0563];
    OCV=1 ./[4.1617,4.0913,4.0749,4.0606,4.0153,3.9592,3.9164,3.8687,3.8163,3.7735,3.7317,3.6892,3.6396,3.5677,3.5208,3.4712,3.3860,3.2880,3.2037,3.0747];
    onebyV_interp = LinearInterpolation(reverse(SOC), reverse(OCV))
    return onebyV_interp
end



function get_one_by_Vsp()
    OCV, min, max = get_OCV_func()
    V(S,P) = (OCV(S) + sqrt( OCV(S)^2 - 4*P*70/1000 )  )/2

    S = 100 .*[1,0.9503,0.9007,0.8510,0.8013,0.7517,0.7020,0.6524,0.6027,0.5530,0.5034,0.4537,0.4040,0.3543,0.3046,0.2550,0.2053,0.1556,0.1059,0.0563]
    P = range(0, stop=20, length = length(S))
    S = reverse(S)
    Vdata = zeros(length(S), length(P))
    for i = 1:size(Vdata,1), j = 1:size(Vdata,2)
        Vdata[i,j] = 1/V(S[i], P[j])
    end
    
    itp = LinearInterpolation((S,P),Vdata)
    return itp

end
function get_one_by_Vsp(OCV, Svec, R, Pmax)
    V(S,P) = (OCV(S) + sqrt( OCV(S)^2 - 4*P*R )  )/2

    # S = 100 .*[1,0.9503,0.9007,0.8510,0.8013,0.7517,0.7020,0.6524,0.6027,0.5530,0.5034,0.4537,0.4040,0.3543,0.3046,0.2550,0.2053,0.1556,0.1059,0.0563];
    P = range(0, stop=Pmax, length = length(Svec))
    Svec = reverse(Svec)
    Vdata = zeros(length(Svec), length(P))
    for i = 1:size(Vdata,1), j = 1:size(Vdata,2)
        Vdata[i,j] = 1/V(Svec[i], P[j])
    end
    
    itp = LinearInterpolation((Svec,P),Vdata)
    return itp

end

function get_one_by_VspGLM(OCV, Svec, R, Pmax)
    V(S,P) = (OCV(S) + sqrt( OCV(S)^2 - 4*P*R )  )/2

    # S = 100 .*[1,0.9503,0.9007,0.8510,0.8013,0.7517,0.7020,0.6524,0.6027,0.5530,0.5034,0.4537,0.4040,0.3543,0.3046,0.2550,0.2053,0.1556,0.1059,0.0563];
    P = range(0, stop=Pmax, length = length(Svec))
    Svec = reverse(Svec)
    Vdata = zeros(length(Svec), length(P))
    for i = 1:size(Vdata,1), j = 1:size(Vdata,2)
        Vdata[i,j] = 1/V(Svec[i], P[j])
    end
    data = DataFrame(X=repeat(Svec, outer = length(P)), Y=repeat(P, outer=length(Svec)), Z =vec(Vdata))
    gm1 = lm(@formula(Z ~ X + Y), data)

    return gm1
end

function model(t0, tf; SOC0 = 90)
    tv = [t0]
    SOCv = [SOC0]
    SOC = SOC0
    Vv = [V0]
    V = V0
    Δ = (tf - t0)/500
    t = t0
    while t < tf

    end
end
function SOC_dyn(t, SOC)
    SOC += Δ*It/Q
    x2  += Δ*(I - 35*Dp*x2/Rp^2)
    return SOC, x2
end
# Make funciton for Eulers Method propogating
function riemman(V, tvec, S0, P, Cmax; N = 50)
    S = S0
    t0, tf = tvec[1], tvec[2]
    Δ =  (tf - t0)/N
    t = t0
    S_vec = Float64[]
    t_vec = Float64[]
    append!(S_vec,S)
    append!(t_vec, t0)
    while true
        S -= P/(Cmax*V(S))*Δ*100 #times 100 for "percent"
        t += Δ

        append!(t_vec, t)
        append!(S_vec, S)
        t >= tf && return S_vec, t_vec
    end
    return 0
end
#function for ohmic drop...
function riemman2(Vf,OCV, tvec, S0, P, Cmax, Req; N = 50)
    S = S0
    t0, tf = tvec[1], tvec[2]
    Δ =  (tf - t0)/N
    t = t0
    S_vec = Float64[]
    t_vec = Float64[]
    append!(S_vec,S)
    append!(t_vec, t0)
    while true
        V = Vf(OCV, S, Req, P)
        S -= P/(Cmax*V)*Δ*100 #times 100 for "percent"
        t += Δ
        append!(t_vec, t)
        append!(S_vec, S)
        t >= tf && return S_vec, t_vec
    end
    return 0
end

function V_model2(OCV, SOC, Req, P)
    V = (OCV(SOC) + sqrt(OCV(SOC)^2 - 4*Req*P))/2
    return V
end

function riemman7(OCV, tvec, S0, P, Cmax, Req, C; N = 50)
    S = S0
    t0, tf = tvec[1], tvec[2]
    Δ =  (tf - t0)/N
    t = t0
    τ = Req*C
    S_vec = Float64[]
    t_vec = Float64[]
    V_vec = Float64[]
    U_vec = Float64[]
    Vk = OCV(S)
    Uk = exp(-Δ/τ)*0 +Req*(1-exp(-Δ/τ))*P/Vk
    Vk1 = 0
    Uk1 = 0
    append!(S_vec,S)
    append!(t_vec, t0)
    append!(V_vec, Vk)
    append!(U_vec, Uk)

    while true
        Uk1 = exp(-Δ/τ)*Uk + Req*(1 - exp(-Δ/τ))*P/Vk
        Vk1 = OCV(S) - Uk1
        S -= P/(Cmax*Vk1)*Δ*100 #times 100 for "percent"
        t += Δ

        append!(t_vec, t)
        append!(S_vec, S)
        append!(V_vec, Vk1)
        append!(U_vec, Uk1)
        Uk = Uk1
        Vk = Vk1
        t >= tf && return S_vec, t_vec
    end
    return 0

end

function simps(xin::Vector, yin::Vector) #simpsons rule to get total charge used given current and time samples
    x,y = xin, yin
    n = length(y)-1

    n % 2 == 0 || (push!(x, x[end]  + (x[end] - x[end-1])/100); push!(y, 0); n += 1)
    length(x)-1 == n || error("`x` and `y` length must be equal")
    h = (x[end]-x[1])/n
    s = sum(y[1:2:n] + 4y[2:2:n] + y[3:2:n+1])
    return h/3 * s
end
