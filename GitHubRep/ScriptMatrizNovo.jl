##########Imports para usar se pela primeira vez no terminal####################

#using Pkg
#Pkg.add("LinearAlgebra")
#Pkg.add("DelimitedFiles")
#Pkg.add("IterTools")
#Pkg.add("TimerOutputs")


##########PODE-SE USAR PARA K!=0, ESTÁ ATUALIZADO######################


using LinearAlgebra
using DelimitedFiles
using IterTools
using TimerOutputs
using Printf

function MatrixTotal()
#############################################
    δ(i,j) = ==(i,j)
    
    Dir = "Data"
    if isdir(Dir)
        cd(Dir)
        println("Enetering Directory")
    else
        mkdir(Dir)
        cd(Dir)
        println("Creating Directory")
    end

    
    ##############################
    
    Forma  = ARGS[1] 
    println("MShape=",Forma)
    Arestasl = ARGS[2] #lista
    b = split(Arestasl,",")
    Arestas = [parse(Float64,b[1][2:end]),parse(Float64,b[2]),parse(Float64,b[3][1:(end-1)])]
    println("Arestas = " * Arestasl)
     
    Nl = ARGS[3]
    N = parse(Int,Nl)
    println("N=" * Nl)
    
    erl = ARGS[4]  #(real,complexo)
    b = split(erl,",")
    er = [parse(Float64,b[1][2:end]),parse(Float64,b[2][1:(end-1)])]
    println("er=" , er)
    
    k_vecl = ARGS[5]
    b = split(k_vecl,",")
    k_vec = [parse(Float64,b[1][2:end]),parse(Float64,b[2]),parse(Float64,b[3][1:(end-1)])]
    println("k=", k_vec)
    
    E_0l = ARGS[6]
    b = split(E_0l,",")
    E_0 = [parse(Float64,b[1][2:end]), parse(Float64,b[2]), parse(Float64,b[3][1:(end-1)])]
    println("E=", E_0)
    
#######################################
    
    
    

    
#########################ID#################################

    function Iden(Forma,Arestas,N,er,k_vec,E_0)
        ID = "M"*Forma
        ID *= "N" 
        
        for i in Arestas
            ID = ID * string(floor(Int,i))
        end
        ID *= "D"* string(N)*"e"
        e1 = er[1]
        e2 = er[2]
        ID*= string(round(e1,digits=3)) *"_" * string(round(e2,digits=3))
        ID*= "k"
        for i in k_vec
            ID = ID * string(floor(Int,i))
        end
        ID *= "E"
        for i in E_0
            ID = ID * string(floor(Int,i))
        end
        ID
    end
    
    ####################################

    Erros_Matrizes = 0
    
    
    ###################VAR################################
    ID = Iden(Forma,Arestas,N,er,k_vec,E_0)
    
    mkpath(ID)
    cd(ID)
    TextoVar = "Var.txt"
    
    touch(TextoVar)
    f = open(TextoVar, "w")

    for x in ARGS
        write(f,string(x)* "\n")
    end
    close(f)
    
    ################Pos Placa#####################
    
    
    
    function PosicoesPlaca(N,Lx,Ly,Lz)
        #N é densidade dipolar linear

        N_x = floor(Int,Lx * (N))
        N_y = floor(Int,Ly * (N))
        N_z = floor(Int,Lz * (N))

        N_dips = (N_x)*(N_y)*(N_z)

        Rs = Array{Float64}(undef,N_dips,3)

        x,y,z = 0:(N_x-1),0:(N_y-1),0:(N_z-1)
        R = product(x,y,z)

        #print(size(R))

        j = 0
        for i in R
            j+=1
            Rs[j,:] .= i 
        end

        TextoRs = "Pos.txt"
        println(TextoRs)
        writedlm(TextoRs, Rs./N .+ 1/(2N) )

        Rs ./ N .+ 1/(2N) , 1/N
    end
    
    #######################################


    function MatrizA(Posicoes, k_vec, ϵ_r, d)

        N_dips = size(Posicoes)[1]
        k = norm(k_vec)

        N_M = 3 * N_dips
        α = 3/(4 * π)*d^3 * (ϵ_r - 1)/(ϵ_r + 2)
        A = zeros(Complex,N_M,N_M)

        for i in 1:(N_dips-1), j in (i+1):N_dips

            r_vec = Posicoes[i,:] - Posicoes[j,:]
            r = norm(r_vec)

            for l in 1:3, m in 1:3
                B = k^2 * (r_vec[l] * r_vec[m] - r^2 * δ(l,m) ) 
                C = 1/r^2 - im* k/r 
                D = r^2 * δ(l,m) - 3*r_vec[l] * r_vec[m]
                A[3*(i-1) + l, 3*(j-1) + m ] = exp(im * k * r)/r^3 * (B + C*D)
            end

        end

        for i in 1:N_M
            A[i,i] = 1/α
        end

        A_f = Symmetric(A)
        A_f
    end

##################Assume-se que o zero das posições está correcto#################################

    function Incident_field(Posicoes,k_vec,E_0)    
        N_dips = size(Posicoes)[1]
        N_M = 3 * N_dips
        E_inc = Array{Complex}(undef,N_M)
        for i in 1:N_dips
            r_vec = Posicoes[i,:]
            E_inc[(i*3 - 2) : i*3 ] .= E_0 * exp(im*dot(r_vec,k_vec)) 
        end

        TextoE =  "EInc.txt"
        writedlm(TextoE,E_inc)

        E_inc
    end


    function CalculoTotal(Posicoes, d, k_vec, ϵ_r,E_0,EM )
        
        er = ϵ_r[1] + ϵ_r[2]im
        Pos = Posicoes ./ d                   #adimensionalizou-se aqui
        α = 3/(4 * π) * (er - 1)/(er + 2)   #e aqui desaparece o d
        k = norm(k_vec)

        A = MatrizA(Pos, k_vec ,er,1)
        N_M = size(A)[1]
        E= Incident_field(Posicoes, k_vec,E_0)

        p  = A \ E    
        Ef = (1/α - 4* π/3) * p    


        TextoE = "E.txt"
        TextoA = "A.txt" 
        TextoP = "P.txt"
        writedlm(TextoE, Ef)
        writedlm(TextoA, A)
        writedlm(TextoP, p)

        EM= norm(A*p - E)
        TextoEM = "Erro.txt"
        writedlm(TextoEM,EM)

    end


    


    function Retrieve_Complex(;Field = "E")    #returns Pol,E

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
        println(EIn)
        N = size(EIn)[1]
        E = zeros(Complex,3*N)
        j=1
        for i in 1:N
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
        println("Complex Fields Retrieved")
    end


    function Complex_Matrix()
        TextoA = "A.txt" 
        TextoAR= "AR.txt"
        TextoAF= "AF.txt"
        AIn = readdlm(TextoA)
        N = size(AIn)[1]
        A = zeros(Complex,N,N)
        for i in 1:N, j in 1:N 
            A[i,j]   += AIn[i,j*3-2]

            if AIn[i,j*3-1]== "+"
                A[i,j] += parse(Float64,AIn[i,j*3][1:end-2])*im
            else
                A[i,j] -= parse(Float64,AIn[i,j*3][1:end-2])*im
            end        
        end
        AR = real.(A)
        AF = imag.(A)
        
        writedlm(TextoAR,ER)
        writedlm(TextoEF,EF)
        println("Matrix Retrieved")
    end



    to = TimerOutput()

    
    Rs,d = PosicoesPlaca(N,Arestas[1],Arestas[2],Arestas[3])

    @timeit to "Calculo N=$N" CalculoTotal(Rs,d,k_vec,er,E_0,Erros_Matrizes)

    Retrieve_Complex(Field = "E")
        
    TextoTime = "Tempos.txt"
    
    print(to)

    f = open(TextoTime,"w")
    print(f,to)
    close(f)

end

MatrixTotal()


