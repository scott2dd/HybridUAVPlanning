function array_to_list(X::Matrix{Float64})
    Xout = Vector{Float64}[]
    for i in 1:size(X,1)
        push!(Xout, X[i,:])
    end
    return Xout
end
function array_to_list(X::Matrix{Int64})
    Xout = Vector{Int64}[]
    for i in 1:size(X,1)
        push!(Xout, X[i,:])
    end
    return Xout
end
function dom2(X,Y)
    #returns true if arg1 >> arg2
    bool = false
    (X[1] <= Y[1] && X[2] <= Y[2]  ) && (bool = true) #&& X[3] <= Y[3]
    return bool 
end

function dom(X::Vector,Y::Vector)
    #returns true if X >> Y
    bool = false
    (X[1] <= Y[1] && X[2] >= Y[2]  && X[3] >= Y[3]) && (bool = true)
    return bool
end

function dom_min(X::Vector,Y::Vector)
    #returns true if X >> Y
    bool = false
    (X[1] <= Y[1] && X[2] <= Y[2]  && X[3] <= Y[3]) && (bool = true)  
    return bool 
end

function make_eff_list(N; maxes = [10,10,10])
    dim = length(maxes)
    M = zeros(N, )
end

function mergeXY(X::Vector{Vector{Float64}}, Y::Vector{Vector{Float64}})
    M = Vector{Float64}[]
    i,j = 1,1
    p,q = length(X), length(Y)
    while true
        if i>p
            M = vcat(M,  Y[j:end])
            return M
        end
        if j > q
            M = vcat(M, X[i:end])
            return M
        end
        Xi, Yj = X[i], Y[j]
        #3a
        if Xi[1] < Yj[1]
            while Xi[2] ≤ Yj[2] && Xi[1] < Yj[1]
                j += 1
                j > q && break 
                Yj = Y[j]
            end
            push!(M, Xi)
            i += 1
        #3b
        elseif Yj[1] < Xi[1]
            while Yj[2] ≤ Xi[2] && Yj[1] < Xi[1]
                i += 1
                i > p && break
                Xi = X[i]
            end
            push!(M, Yj)
            j += 1
        #3c
        elseif Xi[2] < Yj[2]
            while Yj[2] ≥ Xi[2] && Xi[1] == Yj[1]
                j += 1
                j > q && break
                Yj = Y[j]
            end
            push!(M, Xi)
            i +=1   
        #3d
        else
            while Xi[2] ≥ Yj[2] && Xi[1] == Yj[1]
                i += 1
                i > p && break
                Xi = X[i]
            end
            push!(M, Yj)
            j += 1
        end
    end
    return 0
end
function merge_sort3D(Temp, Perm)
    X, Y = Temp, Perm
    M = []
    i,j = 1,1
    p,q = length(X), length(Y)
    while true
        #2a
        if i>p
            M = vcat(M,  Y[j:end])
            return M
        #2b
        elseif j > q
            M = vcat(M, X[i:end])
            return M
        end
        Xi, Yj = X[i], Y[j]
        #3a
        if Xi[1] < Yj[1]
            while Xi[2] ≤ Yj[2] && Xi[1] < Yj[1]
                j += 1
                j > q && break 
                Yj = Y[j]
            end
            push!(M, Xi)
            i += 1
        #3b
        elseif Yj[1] < Xi[1]
            while Yj[2] ≤ Xi[2] && Yj[1] < Xi[1]
                i += 1
                i > p && break
                Xi = X[i]
            end
            push!(M, Yj)
            j += 1
        #3c
        elseif Xi[2] < Yj[2]
            while Yj[2] ≥ Xi[2] && Xi[1] == Yj[1]
                j += 1
                j > q && break
                Yj = Y[j]
            end
            push!(M, Xi)
            i +=1   
        #3d
        elseif Xi[2] > Yj[2]
            while Xi[2] ≥ Yj[2] && Xi[1] == Yj[1]
                i += 1
                i > p && break
                Xi = X[i]
            end
            push!(M, Yj)
            j += 1
            
        # extra for third dimension: if xi1 == yj1 && xi2 == yj2
        elseif Xi[3] < Yj[3]
            while Yj[3] ≥ Xi[3] && Xi[1] == Yj[1] && Xi[2] == Yj[2]
                j += 1
                j > q && break
                Yj = Y[j]
            end
            push!(M, Xi)
            i += 1     
        
        else 
            while Xi[3] ≥ Yj[3] && Xi[1] == Yj[1] && Xi[2] == Yj[2]
                i += 1
                i > p && break
                Xi = X[i]
            end
            push!(M, Yj)
            j += 1
        end
    end
    return 0
end

function brute_merge_sort(X, Y)  #need to adjust for equivalent labels....
    #for equiv labels: know that x == y cannot have duplicates in X or Y.. so just remove 1...
    M = []
    eq_keep = []
    for y in Y
        add_y = true
        for x in X
            if dom_min(x,y)
                add_y = false
                if x == y #if equal, may need to keep 1 
                    #add to y to eq_keep, add_y = true
                    if y ∉ eq_keep #if not in eq_keep, then we have not seen this label before so add it to eq_keep and set bool to true
                        push!(eq_keep, y)
                        add_y = true
                    else #else we have seen in beore, so set add_y to false 
                        add_y = false
                    end
                end
            end
        end
        add_y && push!(M, y) 
    end
    #we do not need to do equality corrections when looping through other "direction"
    #..because we only need to keep 1 of each equivalent set
    for x in X
        keep_x = true 
        for y in Y
            if dom_min(y,x)
                keep_x = false
            end
        end
        keep_x && push!(M, x)
    end
    return M
end


function merge_3D_drew(X, Y) #BROKEN LOGIC!!!!!!
    #need to eventually change this to X::heap and Y::heap
    #1) Sort [X;Y] by C
    M = vcat(X, Y)
    sort!(M, by = x -> x[1])
    println("X cat Y:")
    display(M)
    #2) remove dominated labels
    Bmin = Inf
    Qmin = Inf
    bool = ones(Bool, size(M,1))
    done = zeros(Bool, size(M,1)) #use this to skip in main loop after looping through equal C values
    for idx in 1:length(M)
        done[idx] == 1 && continue
        label = M[idx]
        Bi, Qi = label[2], label[3]
        idx2 = idx + 1
        if idx2 <= size(M,1) && M[idx2][1] == label[1]
            while true
                if idx2 == size(M,1) || M[idx2+1,1] != label[1]
                    break
                end
                idx2 += 1
            end
            println("here")
            temp = M[idx:idx2]
            booltemp = EFF_X(temp, dom_min)
            #set bool slice to bool2 vals
            bool[idx:idx2] = booltemp
            done[idx:idx2] .= 1

            #need to grab new Bmin Qmin...
        elseif Bi >= Bmin && Qi >= Qmin #if inneficient
            #then remove this
            bool[idx] = 0
        else #else is dominated
            Bmin = minimum([Bi, Bmin])
            Qmin = minimum([Qi, Qmin]) 
        end

    end
    display(bool)
    return M[bool]
end


function merge_3D_old(X, Y) #BROKEN LOGIC!!!! FAILS IF EQUAL B AND Q BUT PRIOR C EXISTS
    Mtemp = vcat(X, Y)
    sort!(Mtemp) #sorts by 1st el then 2nd el then 3rd el....
    boolv = ones(Bool, length(Mtemp))
    #now go through each element 
    prior = (Inf, Inf, Inf)
    for i in eachindex(Mtemp)
        label =  Mtemp[i]
        if prior[1] == label[1] 
            if prior[2] == label[2] #we sorted, so if equal 1's and 2's then 3piror <= 3i
                boolv[i] = 0
            else #if 2's are not equal, then can just compare 
                if label[3] < prior[3] #if this label is efficient then update _prior_
                    #1'equal, 2's not.  We sorted, so if _label_[3] < prior[3] then we are eff
                    prior = copy(label)
                else #else if 3rds are equal or _label_3 is more, then set _boolv_[i] to 0 and keep prior as is
                    boolv[i] = 0
                end
            end
        else #if not equal C's, loop through every prior _i_ and check if dommed
            for ii in 1:(i-1)
                boolv[ii] == 0 && continue #if already dommed don't need to compare
                if dom_min(Mtemp[ii], label)
                    boolv[i] = 0
                    break
                end
            end
            if boolv[i] ==  1 #if this label is effieicnet across all others, then this is new _prior_
                prior = copy(label)
            end
        end
    end
    return Mtemp[boolv]

end

function brute_EFF_merge(X, Y)
    M = vcat(X, Y)
    sort!(M)
    boolM = EFF_X(M, dom_min)
    M = M[boolM]
    return M
end
function EFF_X(X, dommy)
    Xout = Vector{Number}[]
    eq_remove = Int[]
    eq_keep = Int[]
    unique!(X)
    bool_vec = ones(Bool, length(X))
    hard_dom = zeros(Bool, length(X))
    for i in 1:length(X)
        x = X[i]
        boolx = true
        for j in 1:length(X)
            xx = X[j] 
            i == j && continue
            if dommy(xx, x) && x != xx #if dominated AND not equal, set as hard dom
                boolx = false
                bool_vec[i] = false
                hard_dom[i] = 1
            end
            if x == xx #if they are equal label values, may need to turn boolx back on
                if i ∉ eq_keep && i ∉ eq_remove && hard_dom[i] == 0 #if we have never seen i before AND not hard dommed  then keep Xi
                    boolx = true; bool_vec[i] = true
                    println("first: ",i)
                    push!(eq_keep, i)
                    push!(eq_remove, j)
                elseif i ∈ eq_keep && hard_dom == 0
                    boolx = true; bool_vec[i] = true
                    push!(eq_remove, j)
                else 
                    boolx = false
                    bool_vec[i] = false
                end
            end
        end
    end
    return bool_vec
end
function EFF_list(Γ::BinaryMinHeap{Tuple{Float64, Vector{Float64}}}, new::Vector{Float64})  #should not need this, but may be faster if we have less? like a merge sort?
    EFF_bool = true
    map_Γ = Γ.node_map
    for i in 1:length(map_Γ)  #γ in Γ
        map_Γ[i] == 0 && continue
        γ = Γ[i][2]
        if dom(γ, new)
            EFF_bool = false
            return EFF_bool
        end
    end
    return EFF_bool
end


