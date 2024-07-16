##########Imports para usar se pela primeira vez no terminal####################

#using Pkg
#Pkg.add("LinearAlgebra")
#Pkg.add("DelimitedFiles")
#Pkg.add("IterTools")
#Pkg.add("TimerOutputs")


##########SÓ USAR PARA APROXIMAÇÃO QUASE ESTÁTICA######################


using LinearAlgebra
using DelimitedFiles
using IterTools
using TimerOutputs
using Printf

###########################VARIÁVEIS####################################

function Total() 

    Forma  = ARGS[1]
    Arestasl = ARGS[2] #lista
    b = split(Arestasl,",")
    Arestas = [parse(Float64,b[1][2:end]),parse(Float64,b[2]),parse(Float64,b[3][1:(end-1)])]
    println(Arestas)
     
    Nl = ARGS[3]
    N = parse(Int,Nl)
    
    erl = ARGS[4]  #(real,complexo)
    b = split(erl,",")
    er = [parse(Float64,b[1][2:end]),parse(Float64,b[2][1:(end-1)])]
    println(er)
    
    k_vecl = ARGS[5]
    b = split(k_vecl,",")
    println(b)
    k_vec = [parse(Float64,b[1][2:end]),parse(Float64,b[2]),parse(Float64,b[3][1:(end-1)])]
    println(k_vec)
    
    E_0l = ARGS[6]
    println(E_0l)
    b = split(E_0l,",")
    println(b)
    E_0 = [parse(Float64,b[1][2:end]) parse(Float64,b[2]) parse(Float64,b[3][1:(end-1)])]
    println(E_0)
    
    deltaEl = ARGS[7]
    deltaE = parse(Float64,deltaEl)
    Lambdal = ARGS[8]
    Lambda = parse(Float64,Lambdal)
    
    Dir = "Data"
    cd(Dir)


    #########################ID#################################

    function Iden(Forma,Arestas,N,er,k_vec,E_0,deltaE,Lambda)
        ID = "F"*Forma
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
        str = @sprintf "%.0E" deltaE
        ID *= "d"*str*"l"*string(Lambda)
        ID
    end

    
    
    ###################VAR################################
    ID = Iden(Forma,Arestas,N,er,k_vec,E_0,deltaE,Lambda)
    
    mkpath(ID)
    cd(ID)
    TextoVar = "Var.txt"
    
    touch(TextoVar)
    f = open(TextoVar, "w")

    for x in ARGS
        write(f,string(x)* "\n")
    end
    close(f)
    
    
       #####################POS#######################

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
    
    
    ##################ER e EF############################

    
    

function PontoFixo_Max(Rs,d, E_0, δE, er; λ = 0.5)
    ϵ_r = er[1] + er[2]*im
    N_r = size(Rs)[1]

    α = 3*d^3/(4* π) * (ϵ_r -1)/(ϵ_r +2)   #polarização não adimensionalizada! Mudar isto!
    polarizacoes = zeros(Complex,N_r,3)
    for i in 1:N_r
        polarizacoes[i,:] = α * E_0 
    end
    Campo_eletrico = zeros(Complex,N_r,3)
    for i in 1:N_r
        Campo_eletrico[i,:] = E_0
    end

    Campo_2 = zeros(Complex,N_r,3)
    
    f = 0
    
    C = abs.(Campo_2- Campo_eletrico)
    max = C[argmax(C)]
    maximos = []
        
    print("OK1")
    
    while max > δE
        
        f+=1
        Campo_2 .= Campo_eletrico

        for i in 1:N_r
            polarizacoes[i,:] = λ .* α .* Campo_eletrico[i,:]  + (1-λ) .* polarizacoes[i,:]
        end

        for i in 1:N_r
            soma = zeros(Complex,3)
            for j in 1:N_r
                if j != i
                    rvec = Rs[i,:] - Rs[j,:]
                    rnorm = norm(rvec)
                    soma .+= 1/(rnorm^5) * ( 3 * (polarizacoes[j,1] *rvec[1] + polarizacoes[j,2] *rvec[2] + polarizacoes[j,3] *rvec[3]) .* rvec  - rnorm^2 .* polarizacoes[j,:]) 
                end
            end
            Campo_eletrico[i,:] = E_0[1,:] + soma
        end
        
        C = abs.(Campo_2- Campo_eletrico)
        max = C[argmax(C)]
        push!(maximos,max)
        
        if max > 5
            println("ERRO")
            TextoE = "E.txt"
            Campo_eletrico = zeros(Complex,N_r,3)
            writedlm(TextoE, Campo_eletrico)
            return()
            
        end
    end
    
    print("OK2")
    push!(maximos,f)
    
    Campo_eletrico = Campo_eletrico - 4 * π .* polarizacoes ./ (3*d^3)
    
    TextoE ="E.txt"
    TextoP = "P.txt"
    
    writedlm(TextoP, polarizacoes)
    writedlm(TextoE, Campo_eletrico)
    
    TextoM = "Max.txt"
    writedlm(TextoM,maximos)
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
    
    print(TextoER)
    
    writedlm(TextoER,ER)
    writedlm(TextoEF,EF)
    E
    print("OK3")
end

if isfile("Pos.txt")
    Rs = readdlm("Pos.txt")
    d = norm(Rs[2] -Rs[1])
        
    println(d)
    println(N)
else
    Rs,d = PosicoesPlaca(N, Arestas[1], Arestas[2], Arestas[3]) 
end
    
to = TimerOutput()

@timeit to "Calculo" PontoFixo_Max(Rs,d, E_0, deltaE, er; λ=Lambda)
            
Retrieve_Complex(Field = "E")


TextoTime = "Tempos.txt"

f = open(TextoTime,"w")
print(f,to)
close(f)


end 

Total()