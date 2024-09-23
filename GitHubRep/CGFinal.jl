function Total(args;Eguess,Pguess) 
    #Arguments of the function come in the following array:
    #"NameShape,PosFile,d,PosName,dieletric,kvec,ElectricField,Error, Guess"
    #Ideally the PosFile is in the Dir Data/NameShape/PosName, but it's not required
    #The Pos vector should come in the shape [N,3]
    #d is the distance between gridpoints in the cubic grid
    #The Dielectric,kvec,EletricField,Error are assumed to be inputed as strings
    #The Dieletric should be in the shape (real,complex)
    #The EF should be in the shape [(e1,e1c),(e2,e2c),(e3,e3c)]
    #The Kvec should be in the shape [k1,k2,k3]
    #Error is compared to the modulus of (A'Ax- A'E_inc)
    
    #The Eguess,Pguess optional are a basis for interpolation for a initial 
    #guess for the conjugated gradients.
    #The Eguess is interpolated using Pguess to the Pos, and is then used as
    #the initial eletric field.
    
    #############################################
    delta(i,j) = ==(i,j) #KroneckerDelta
    
    Dir = "Data"         #This is the default directory. 
    
    function CheckDir(Dir)
        if isdir(Dir)
            cd(Dir)  
        else
            mkdir(Dir)
            cd(Dir)
        end
    end
    
    CheckDir(Dir)

    ######################################

    Shape  = args[1] 
    Pos = readdlm(args[2])  #Getting the Positions from the Pos file
    N = size(Pos)[1]
    d = parse(Float64,args[3])
    PosName = args[4]
    erl = args[5]  
    b = split(erl,",")
    er = [parse(Float64,b[1][2:end]),parse(Float64,b[2][1:(end-1)])]
    
    k_vecl = args[6]
    b = split(k_vecl,",")
    k_vec = [parse(Float64,b[1][2:end]),parse(Float64,b[2]),parse(Float64,b[3][1:(end-1)])]
    
    E_0l = args[7]
    b = split(E_0l,",")
    E_0 = [(parse(Float64,b[1][3:end]),parse(Float64,b[2][1:end-1])), (parse(Float64,b[3][2:end]),parse(Float64,b[4][1:end-1])),(parse(Float64,b[5][2:end]),parse(Float64,b[6][1:(end-2)]))]
    
    deltaEl = args[8]
    deltaE = parse(Float64,deltaEl)
    
    #println("Variables:",args)
    
    #########################################
    
    function Iden(N,er,k_vec,E_0,deltaE)  #Creates the Dir Name for the Results
        ID = "CGF_"
    
        ID *= "N_"*string(N)*"e"
        e1 = er[1]
        e2 = er[2]
        ID*= string(round(e1,digits=3)) *"_" * string(round(e2,digits=3))
        ID*= "_k"
        for i in k_vec
            ID = ID * string( round(i[1],digits=2))
        end
        ID *= "E"
        for i in E_0
            ID = ID * string( (round(i[1],digits=2),round(i[2],digits=2)))
        end
        str = @sprintf "%.0E" deltaE
        ID *= "d"*str
        return ID
    end

    Dir = Shape
    CheckDir(Shape)
    Dir = PosName
    CheckDir(PosName)
    
    ID = Iden(N,er,k_vec,E_0,deltaE)
    
    println(ID)
    
    mkpath(ID)
    cd(ID)
    TextoVar = "Var.txt"
    
    touch(TextoVar)
    f = open(TextoVar, "w")

    for x in args
        write(f,string(x)* "\n")
    end
    close(f)
    
##########################
    
    function A_Mult(Positions, knorm, Vector, d, er; conjugation = 0)
    
        #I will be assuming that Vector is in the shape (N_dips,3)
        #This will give the matrix multiplication in the same shape

        alpha = 3/(4 * pi)*d^3 * (er - 1)/(er + 2) #dieletric constant
        k = knorm                                  #omega/c
        N_dips = size(Positions)[1]
        result = zeros(Complex,N_dips,3)

        if conjugation == 0
            for i in 1:(N_dips-1)
                for j in (i+1):N_dips         #Since the matrix is symmetric

                    Rvec = Positions[i,:] - Positions[j,:]
                    r = norm(Rvec)

                    for l in 1:2, m in l+1:3   #symmetry of the 3x3 blocks
                        A = k^2 * Rvec[l]*Rvec[m]/r^3
                        B = (1/r^3 - im*k/r^2)*(-3*Rvec[l]*Rvec[m]/r^2)
                        Mterm = exp(im*k*r)*(A+B)
                        result[i,l] += Mterm * Vector[j,m]
                        result[i,m] += Mterm * Vector[j,l]
                        result[j,l] += Mterm * Vector[i,m]
                        result[j,m] += Mterm * Vector[i,l]
                    end

                    for l in 1:3   #diagonal of the 3x3 blocks
                        A = k^2 * (Rvec[l]^2 - 1)/r^3
                        B = (1/r^3 - im*k/r^2) * (1 -3*Rvec[l]^2/r^2)
                        Mterm = exp(im*k*r)*(A+B)
                        result[i,l] += Mterm*Vector[j,l]
                        result[j,l] += Mterm*Vector[i,l]
                    end
                end
            end

            for i in 1:N_dips, j in 1:3      #diagonals of the matrix
                result[i,j] += 1/alpha * Vector[i,j]
            end

            else  #conjugation of the previous matrix multiplication
            for i in 1:(N_dips-1)
                for j in (i+1):N_dips

                    Rvec = Positions[i,:] - Positions[j,:]
                    r = norm(Rvec)

                    for l in 1:2, m in l+1:3
                        A = k^2 * Rvec[l]*Rvec[m]/r^3
                        B = (1/r^3 + im*k/r^2)*(-3*Rvec[l]*Rvec[m]/r^2)
                        Mterm = exp(-im*k*r)*(A+B)
                        result[i,l] += Mterm * Vector[j,m]
                        result[i,m] += Mterm * Vector[j,l]
                        result[j,l] += Mterm * Vector[i,m]
                        result[j,m] += Mterm * Vector[i,l]
                    end

                    for l in 1:3
                        A = k^2 * (Rvec[l]^2 - 1)/r^3
                        B = (1/r^3 + im*k/r^2) * (1 -3*Rvec[l]^2/r^2)
                        Mterm = exp(-im*k*r)*(A+B)
                        result[i,l] += Mterm*Vector[j,l]
                        result[j,l] += Mterm*Vector[i,l]
                    end
                end
            end
            for i in 1:N_dips, j in 1:3
                result[i,j] += conj(1/alpha) * Vector[i,j]
            end
        end

        return result
    end
   

##############################
    function Incident_field(Posicoes,k_vec,E_0) 
        #Returns an Incident field from the E_0
        #It gives it in the shape (N,3)
        
        k = zeros(3)
        E = zeros(Complex,3)
        for i in 1:3
            k[i] = k_vec[i]
            E[i] = E_0[i][1] + im*E_0[i][2]
        end
            
        N_dips = size(Posicoes)[1]
        E_inc = Array{Complex}(undef,N_dips,3)
    
        for i in 1:N_dips
            r_vec = Posicoes[i,:]
            E_inc[i,:] .= E * exp(im*dot(r_vec,k)) 
        end

        E_inc
    end

     ############################
    
     function DDAMemoryLessCG(y,x0,error;Positions,d,e,knorm)
    #This function solves the system Ax = y

        function AMult(Vector)
            return A_Mult(Positions, knorm, Vector, d, e; conjugation= 0)
        end
        function CMult(Vector)
            return A_Mult(Positions, knorm, Vector, d, e; conjugation= 1)
        end

        z = CMult(y)
        
        g = z - CMult(AMult(x0))
        p = g
        w = AMult(x0)
        v = AMult(p)
        x=x0

        i = 0
        while norm(p) > error
            i+=1
            alpha = dot(g,g)/dot(v,v)
            x += alpha * p

            if i%10!=0
                w += alpha * v
            else
                println("Norm, #iteration:")
                println((norm(p),i))
                w = AMult(x)
                
            end
            ng = dot(g,g)
            g = z - CMult(w)
            beta = dot(g,g)/ng
            p = g + beta*p
            if i%10!=0
                v = AMult(g) + beta*v
            else
                v = AMult(p)
            end
        end
       return x
    end
    
    
    ######################
    function Retrieve_Complex(;Field = "E")    
    #returns E or Pol
    #Also is the function that takes the solution and splits it in
    #Real and Complex Parts for the python scripts
    
    TextoE = "E.txt" 
    TextoP = "P.txt"
    
    if Field == "E"
        TextoER = "ER.txt" 
        TextoEF = "EF.txt" 
    elseif Field == "P"
        TextoE = TextoP
        TextoER = "PR.txt" 
        TextoEF = "PF.txt"           
    elseif Field == "Ei"
        TextoE =  "EInc.txt"
        TextoER = "EIncR.txt"
        TextoEF = "EIncF.txt"        
    end
    
    EIn = readdlm(TextoE)
    N = size(EIn)[1]
    E = zeros(Complex,3*N)
    for i in 1:N, j in 1:3
        E[3*i-3+j]   += EIn[i,j*3-2]
        if EIn[i,j*3-1]== "+"
            E[3*i-3+j] += parse(Float64,EIn[i,j*3][1:end-2])*im
        else
            E[3*i-3+j] -= parse(Float64,EIn[i,j*3][1:end-2])*im
        end
    end
        
    ER = real.(E)
    EF = imag.(E)
    
    writedlm(TextoER,ER)
    writedlm(TextoEF,EF)
    return(E)
    end
    
    ########################
    
    x0 = zeros(N,3)
    Einc = Incident_field(Pos,k_vec,E_0)

    polarizations = DDAMemoryLessCG(Einc,x0,deltaE;Positions=Pos, d = d, e = er[1]+im*er[2], knorm = norm(k_vec))
    
    E = 4*π/(d^3*(er[1]+im*er[2]-1)) .*polarizations  #Clausius Mossoti Relation
    
    ##################################
    TextoE ="E.txt"
    TextoP = "P.txt"
    
    writedlm(TextoP, polarizations)
    writedlm(TextoE, E) 
    
    Retrieve_Complex(Field = "E")
    
end