"""
computing the asymptotic constants Qₚ and θ
Newton iterates for solving x+2x²
version 22-03-25
"""

function mₛ(x)         # the mantissas of BigFloats may be too large to print; we extract a few first values
    chb = 5; che = 7  # no. of chars to extract from beginning (mantissa) resp. end (exponent)
    w = string(x)
    if length(w)>15; w = w[1:chb] * "..." * w[end-che:end] end
    return w
end
function 𝑓(x::BigFloat)    # 𝑓 obtained by \itf, italic f
    return x + 2x^2
end
function 𝑓′(x::BigFloat)
    return 1 + 4.0x
end
function 𝑓″(x::BigFloat)
    return 4.0
end
myprec = 10^5 #113 for quadruple, 52 for double precision
setprecision(myprec)
println("############################","\n","epsBigFl=","\n",mₛ(eps(BigFloat)),"\n")
N = 17  # max no of Newton iterations
p = BigFloat(2.0)

x⃰ = BigFloat(0.0)  #the solution
θ⃰ = BigFloat(0.5) # conjecture!

x = zeros(BigFloat,N)
x[1] = 0.5
e, fx, f′x = zeros(BigFloat,N), zeros(BigFloat,N), zeros(BigFloat,N)

# #compute the Newton iterates; we store the values of f, f' in the vectors  fx, f'x
fx[1] = 𝑓(x[1])
f′x[1] = 𝑓′(x[1])
for k =2:N
    x[k] = x[k-1] - fx[k-1]/f′x[k-1]
    fx[k] = 𝑓(x[k])
    f′x[k] = 𝑓′(x[k])
end

Qₚ, Qₚ′, Q̂ₚ′, Qₚ″, Q̂ₚ″ = zeros(BigFloat,N), zeros(BigFloat,N), zeros(BigFloat,N), zeros(BigFloat,N), zeros(BigFloat,N)

println("time when computing e, s, Qₚ, Qₚ′, Qₚ″ and Q̂ₚ″ in parallel")
@time begin
    Threads.@threads for k=1:N
        e[k] = abs(x⃰-x[k])
        fx[k] = abs(fx[k])
        f′x[k] = abs(f′x[k])
    end
    s = zeros(BigFloat,N)
    Threads.@threads for k=2:N
        s[k] = abs(x[k]-x[k-1])
    end
    Threads.@threads for k=2:N
        Qₚ[k] = e[k]/e[k-1]^p
        Qₚ′[k] = s[k]/s[k-1]^p
        Qₚ″[k] = fx[k]/fx[k-1]^p*f′x[k]
        Q̂ₚ″[k] = fx[k]/fx[k-1]^p*f′x[k-1]^2/f′x[k]
    end
end

println("time when computing e, s, Qₚ, Qₚ′, Qₚ″ and Q̂ₚ″ in vectorized form")
@time begin
    e = abs.(x⃰.-x)
    s[2:N] = abs.(x[2:N] - x[1:N-1])
    fx = abs.(fx)
    f′x = abs.(f′x)
    Qₚ[2:N] = e[2:N]./e[1:N-1].^p
    Qₚ′[2:N] = s[2:N]./s[1:N-1].^p
    Qₚ″[2:N] = fx[2:N]./fx[1:N-1].^p.*f′x[2:N]
    Q̂ₚ″[2:N] = fx[2:N]./fx[1:N-1].^p.*f′x[1:N-1].^2 ./f′x[2:N]
end

println("##############################################################")

Qₚ⃰ = abs(𝑓″(x⃰)/(2𝑓′(x⃰)))
println("\nQₚ⃰= ",mₛ(Qₚ⃰),"\n")

# #print values
println("############################","\n","Qₚ","\n")
for k = N-7:N  println("k=",k-1, "  Qₚ=",mₛ(Qₚ[k]),"  Qₚ′=",mₛ(Qₚ′[k]),"  Qₚ″=",mₛ(Qₚ″[k]),"  Q̂p″=",mₛ(Q̂ₚ″[k]))  end
# #print errors
println("########","\n","errors in approximating Qₚ⃰: |Qₚ⃰-Qₚ|, |Qₚ⃰-Qₚ′|, |Qₚ⃰-Qₚ″|, |Qₚ⃰-Q̂p″|","\n")
for k = N-7:N  println("k=",k-1, "  er Qₚ=",mₛ(Qₚ⃰-Qₚ[k]),"  er Qₚ′=",mₛ(Qₚ⃰-Qₚ′[k]),"  er Qₚ″=",mₛ(Qₚ⃰-Qₚ″[k]),"  er Q̂p″=",mₛ(Qₚ⃰-Q̂ₚ″[k]))  end #one should insert 'abs'

println("############################","\n","errors eₖ=|x⃰-xₖ|","\n")
for k = N-7:N  println("k=",k-1,"   eₖ=",mₛ(e[k]))  end
println("############################","\n","corrections sₖ","\n")
for k=N-7:N  println("k=",k-1,"   sₖ=",mₛ(s[k])) end
println("############################","\n","nonlinear residuals fₖ","\n")
for k=N-7:N  println("k=",k-1,"   fₖ=",mₛ(fx[k])) end
println("############################","\n","errors eₖ approximated by fₖ/f'ₖ","\n")
for k=N-7:N  println("k=",k-1,"   fₖ/f'ₖ=",mₛ(fx[k]/f′x[k])) end

println("=============================","\n","θₖ's","\n")
# # it is not easy to put the hat on θ (and Greek symbols);
# # although possible (see the JuliaMono fonts at https://mono-math.netlify.app/#JuliaMono) this requires installation of some special programs
# # meanwhile we use the index "i" (improved) instead of the hat, and t̂ instead of theta hat
θ, θᵢ, θ′, θᵢ′, θ″, θᵢ″  = zeros(BigFloat,N), zeros(BigFloat,N), zeros(BigFloat,N), zeros(BigFloat,N), zeros(BigFloat,N), zeros(BigFloat,N)
Threads.@threads for i = 2:N
    θ[i] = (e[i])^(1/p^(i-1))
    θᵢ[i] = (e[i]*Qₚ[i])^(1/p^(i-1))
    θ′[i] = (s[i])^(1/p^(i-2))
    θᵢ′[i] = (s[i]*Qₚ′[i])^(1/p^(i-2))
    θ″[i] = (fx[i])^(1/p^(i-1))
    θᵢ″[i] = (fx[i]*Q̂ₚ″[i]/f′x[i])^(1/p^(i-1))
end
for i=N-7:N
    # #print value
    # println("i=",i-1," θ=",mₛ(θ[i]),"  θᵢ=",mₛ(θᵢ[i]),"  θ′=",mₛ(θ′[i]),"  θᵢ′=",mₛ(θᵢ′[i]),"  θ″=",mₛ(θ″[i]),"  θᵢ″=",mₛ(θᵢ″[i]))
    # #print error
    println("i=",i-1," er θ=",mₛ(θ⃰-θ[i])," er θᵢ=",mₛ(θ⃰-θᵢ[i])," er θ′=",mₛ(θ⃰-θ′[i])," er θᵢ′=",mₛ(θ⃰-θᵢ′[i])," er θ″=",mₛ(θ⃰-θ″[i])," er θᵢ″=",mₛ(θ⃰-θᵢ″[i]))
end;  println("--")

xₖ⃰ = zeros(BigFloat,N)   #extrapolated values using the asymptotic constants, e≈c⋅θ^p^k
xₖ⃰f = zeros(BigFloat,N) #extrapolated values using error estimation e≈f/f'
Threads.@threads for i = 2:N
    xₖ⃰[i] = x[i] - θᵢ″[i]^(p^(i-1))/Q̂ₚ″[i]  #/̄
    xₖ⃰f[i] = x[i] - fx[i]/f′x[i]
end
eₖ⃰f = abs.(x⃰.-xₖ⃰f) #extrapolated by f
for i = N-5:N   println("k=",i-1,"  er xₖ⃰f= ",mₛ(eₖ⃰f[i]))  end
eₑ = abs.(x⃰.-xₖ⃰) #extrapolated by c, θ
for i = N-5:N   println("i=",i-1,"  er xₖ⃰= ",mₛ(eₑ[i]))  end
