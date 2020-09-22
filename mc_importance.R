###3 Monte Carlo e Importance Sampling
##U1, U2.. ~ Unif(0,1) e X = {menor n | Sum(Ui) > 1}


######a) Estime E(X). Qual o erro padrão desta estimativa? 
#Seja g(U) = X  = {menor n | Sum(Ui) >1 }
#Eg(U) = integral(g(U)f(U)du)

##Monte Carlo:
#i) simular U1.. Un ~ Unif(0,1)

#ii) E(g(U)) = 1/10000 sum(g(U))

ests <- numeric(10000)
for (i in 1:10000){  ##para cada RE = 10000 sims, 
  sims <- runif(30, 0, 1) #simulamos draws de uma uniforme padrao
  sum <- 0
  j <- 0
  while(sum <= 1){    #calculamos quantos draws sao necessarios tal que a soma seja maior que 1
    sum <- sum + sims[j+1]
    j <- j+1
  }
  ests[i] <- j #guardamos para i in RE a quantia de runifs necessárias 
}
hist(ests)
EXsim <- (1/10000) * sum(ests) #calculamos o valor esperado 

####estimativa de X
EXsim
#[1] 2.7124

####erro padrão da est
sd_est <- (1/10000) * sqrt(sum((ests - EXsim)^2))
sd_est
#[1] 0.008718292


######b) Estimar P(X>=10). Vamos utilizar importance sampling, uma vez que f(U) ~ Unif(0,1). Para X>=10,
#sampling da uniforme vai nos dar pouquissimos valores que correspondem a este caso. 
#Assim, devemos escolher uma distribuição de importância que nos dê uma maior oportunidade de samplearmos de valores com 'n' alto
#Ou seja, devemos escolher uma distribuição com o mesmo suporte de U(0,1) mas que produza VAs pequenas;
#Vamos considerar entao dist Beta com alfa = 2 e beta = 5
#h(x) = P(X>=10) = E(I(sum(X1..X9) <= 1 ))


p1 = 2 #parametros da beta
p2 = 9
k = 9 #E(I(sum(X1..X9) <= 1))
RE = 10000

#vamos ver como a dist se comporta em relação a números menores que é o que nos interessa
beta = rbeta(RE, p1, p2)
hist(beta)


sim = matrix(nrow = RE, ncol = k)

#gerando realizações
for (i in 1:RE) {
  sim[i,] = rbeta(k, p1, p2) 
}

#h(x) é a indicadora de que as 9 variáveis somam menos que 9 
hx = rowSums(sim) <= 1

#f(x) é a densidade da unif padrao
fx = 1

#g(x) é a conjunta da beta p/ 1...9
fdp_beta = apply(sim, 2, function (x) { dbeta(x, p1, p2) })
gx = apply(fdp_beta, 1, prod)

#media
est_vec = (hx*fx)/gx
est = mean(est_vec)
est
# [1] 2.71119972633006e-06

#desvio
sd_est <- (1/RE) * sqrt(sum((est_vec - est)^2))
sd_est
#[1] 5.21503e-07



########### c) Estimar P(X = 10) e calcular erro padrão da estimativa:

#g(x) e f(x) são semelhantes a questão b) 
#porém h(X) = P(X=10) = E(I(X1..X9 < 1) & I(X1..X10 > 1))
#ou seja, estamos interessados no caso em que X1..X9 não ultrapassa 1, mas X1..X10 ultrapassa (assim, temos P(X=10))

p1 = 2 #parâmetros da beta
p2 = 9
k = 10 #E(I(X1..X9 < 1) & I(X1..X10 > 1))
RE = 10000


sim = matrix(nrow=RE, ncol=10)


for (i in 1:RE) {
  sim[i,] = rbeta(k, p1, p2) 
}

fx = 1

hx = (rowSums(sim[,1:9] < 1) & (rowSums(sim) > 1)) ##nova indicadora

gx = apply(dbeta(sim, p1,p2), 1, prod)

est_vec = (hx*fx)/gx
est = mean(est_vec)
est
#[1] 0.02623771


#desvio
sd_est <- (1/RE) * sqrt(sum((est_vec - est)^2))
sd_est
#[1] 0.01230647
















# parametros da beta escolhidos:
shape1 = 1.6
shape2 = 12

# exemplos de valores gerados com essa beta, sao pequenos:
rbeta(n=10, shape1=shape1, shape2=shape2)

# gera 10 realizacoes dessa beta, de tamanho RE:
betas = matrix(nrow=RE, ncol=10)
for (i in 1:10)
{
  betas[,i] = rbeta(n=RE, shape1=shape1, shape2=shape2 )
}

# densidade da uniforme:
fx = 1

# indicadora:
hx = (rowSums(betas[,1:9] < 1) & (rowSums(betas) > 1) )

gx = apply(dbeta(betas, shape1=shape1, shape2=shape2), 1, prod)

estimativa = hx * fx / gx

print(paste('P(X = 10) estimado:', mean(estimativa) ))
print(paste('Erro padrao estimado:', sqrt( var(estimativa) / RE )  ))




