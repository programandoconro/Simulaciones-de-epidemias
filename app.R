###################################################################################################  
# An R script to perform a stochastic epidemic simulation using an Agent Based model
    # (with homogenous mixing equivalent to an SIR model)
    #
    # R script Author: Sherry Towers
    # smtowers@asu.edu
    # Created: Feb 13, 2016
    #
    # Copyright Sherry Towers, 2016
  
####################################################################################################
  # A Shiny App to easily change the parameters of the model and visualize the results in the Browser:
    #
    # App developer: Rodrigo Díaz Lupanow
    # programandoconro@gmail.com
    # Created: March 13, 2019
    #
    # Copyright Rodrigo Díaz-Lupanow, 2019
####################################################################################################
  
    # This script is not guaranteed to be free of bugs and/or errors
    # This script can be freely used and shared as long as the authors and
    # copyright information in this header remain intact.

#####################################################################################################
    library (shiny) 

    ui <- fluidPage(           

    # Título de la app
    titlePanel(“Simulación de epidemias a partir de Modelos Basados en Agentes”),

    # La barra lateral con los sliders que nos interesan
    sidebarLayout(
    sidebarPanel(
    sliderInput(“bins”,
    “Tamaño poblacional:”,
    min = 1,
    max = 1000000,
    value = 10000),

    sliderInput(“bins2”,
    “Número de individuos infectados inicialmente:”,
    min = 1,
    max = 100,
    value = 10),
    sliderInput(“bins3”,
    “Número de simulaciones:”,
    min = 1,
    max = 100,
    value = 10),
    sliderInput(“bins4”,
    ” Período de recuperacion (días^{-1}):”,
    min = 0,
    max = 1,
    value = 1/3),
    sliderInput(“bins5”,
    “Tiempo final del experimento (días):”,
    min = 120,
    max = 1000,
    value = 120),
    sliderInput(“bins6”,
    “Fuerza hipotética del virus:”,
    min = 0,
    max = 10,
    value = 1.5,step = 0.1)
    ),

    # Para mostrar la gráfica final
    mainPanel(
    plotOutput(“distPlot”),

    #para introducir texto importante
    h5(“Basado en el código obtenido a partir de: Stochastic Epidemic Simulation using an Agent Based model (with homogenous mixing equivalent to an SIR model), por Sherry Towers, disponible en: http://sherrytowers.com/2016/02/26/simple-agent-based-disease-modelling-with-homogenous-mixing/&#8221;),
    h5(“Aplicación creada por Rodrigo Díaz para el seminario de la materia: Introducción a la Ecología de Ecosistemas, del Postgrado en Ecología. Prof. encargado: Carlos Méndez, Instituto Venezolano de Investigaciones Científicas (IVIC)”)

    )
    )
    )

    # Define el servidor para realizar la gráfica final
    server <- function(input, output) {

    #a continuación el script modificado

    SIR_agent = function(N # population size
    ,I_0 # initial number infected
    ,S_0 # initial number susceptible
    ,gamma # recovery rate in days^{-1}
    ,R0 # reproduction number
    ,tbeg # begin time of simulation
    ,tend # end time of simulation
    ,delta_t=0 # time step (if 0, then dynamic time step is implemented)
    ){

    # begin by reverse engineering the transmission rate, beta, from
    # R0 and gamma
    # beta is the number of contacts sufficient to transmit infection per unit time
    # Thus, in time delta_t, a susceptible person will contact a Poisson random
    # number of infected people, with mean beta*I*delta_t/N
    # The probability an infected person will recover in time
    # delta_t is p=1-exp(-gamma*delta_t)
    # (note that p does not depend on t! This is due to the
    # memoryless nature of the Exponential distribution)
    #########################
    beta = R0*gamma
    #####################
    # now set up the state vector of the population
    # vstate = 0 means susceptible
    # vstate = 1 means infected
    # vstate = 2 means recovered
    # randomly pick I_0 of the people to be infected
    #####################################
    vstate = rep(0,N)
    index_inf = sample(N,I_0) # randomly pick I_0 people from N
    vstate[index_inf] = 1 # these are the infected people
    ######################################
    # now begin the time steps
    ########################################
    t = tbeg
    S = S_0
    I = I_0
    vS = numeric(0) # this will be filled with the number of susceptibles in the population over time
    vI = numeric(0) # this will be filled with the number of infectives in the population over time
    vtime = numeric(0) # this will be filled with the time steps
    while (t<tend&I>0){ # continue the simulation until we have no more infectious people, or t>=tend
    S = length(vstate[vstate==0]) # count the number of susceptibles, based on the state vector
    I = length(vstate[vstate==1]) # count the number of infectives, based on the state vector
    vS = append(vS,S) # append this info to vectors that we will return from the function
    vI = append(vI,I)
    vtime = append(vtime,t)
    #cat(t,S,I,”\n”)
    deltat=delta_t
    if (delta_t==0){ # this is the calculation of a dynamic time step
    deltat = 1/(beta*I*S/N + gamma*I)
    }
    recover_prob = (1-exp(-gamma*deltat)) # the probability an infective recovers in this time step
    # sample Poisson random numbers of infected people contacted by each person
    avg_num_infected_people_contacted = beta*I*deltat/N
    vnum_infected_people_contacted = rpois(N,avg_num_infected_people_contacted)
    vprob = runif(N) # sample uniform random numbers
    vnewstate = vstate # copy the state vector to a temporary vector used for calculations

    # Infected people recover if the sampled uniform random
    # number is less than the recovery probability
    vnewstate[vstate==1&vprob<recover_prob] = 2
    # If a susceptible contacted at least one infective, they are infected
    vnewstate[vstate==0&vnum_infected_people_contacted>0] = 1
    vstate = vnewstate # update the state vector

    t = t + deltat # update the time
    }
    final = 0
    if (length(vS)>0) final = 1-min(vS)/N

    return(list(time=vtime,I=vI,S=vS,final_size=final))

    }####

    SIR_agent_erlang = function(N # population size
    ,I_0 # initial number infected
    ,S_0 # initial number susceptible
    ,gamma # recovery rate in days^{-1}
    ,k # shape parameter of the erlang distribution
    ,R0 # reproduction number
    ,tbeg # begin time of simulation
    ,tend # end time of simulation
    ,delta_t=0 # time step (if 0, then dynamic time step is implemented)
    ){
    cat(N,I_0,S_0,gamma,k,R0,tbeg,tend,delta_t,”\n”)
    ################################
    # begin by reverse engineering the transmission rate, beta, from
    # R0 and gamma
    #
    # beta is the number of contacts sufficient to transmit infection per unit time
    #
    # Thus, in time delta_t, a susceptible person will contact a Poisson random
    # number of infected people, with mean beta*I*delta_t/N
    #
    # The probability an infected person will recover in time
    # delta_t is p=1-exp(-gamma*delta_t)
    # (note that p does not depend on t! This is due to the
    # memoryless nature of the Exponential distribution)
    ################################
    beta = R0*gamma
    #################################
    # now set up the state vector of the population
    # vstate = 0 means susceptible
    # vstate = 1 means infected
    # vstate = 2 means recovered
    # randomly pick I_0 of the people to be infected
    ################################
    vstate = rep(0,N)
    index_inf = sample(N,sum(I_0)) # randomly pick I_0 people from N
    vstate[index_inf] = 1 # these are the infected people
    #################################
    # randomly set up the sojourn time spent so far in each state.
    # when an individual leaves a state, the sojourn time gets reset to 0
    ################################
    vsojourn = rep(0,N)
    ################################
    # now begin the time steps
    ################################
    t = tbeg
    S = S_0
    I = sum(I_0)
    vS = numeric(0) # this will be filled with the number of susceptibles in the population over time
    vI = numeric(0) # this will be filled with the number of infectives in the population over time
    vtime = numeric(0) # this will be filled with the time steps
    while (t<tend&I>0){ # continue the simulation until we have no more infectious people, or t>=tend
    S = length(vstate[vstate==0]) # count the number of susceptibles, based on the state vector
    I = length(vstate[vstate==1]) # count the number of infectives, based on the state vector
    vS = append(vS,S) # append this info to vectors that we will return from the function
    vI = append(vI,I)
    vtime = append(vtime,t)
    #cat(t,S,I,”\n”)

    deltat=delta_t
    if (delta_t==0){ # this is the calculation of a dynamic time step
    deltat = 1/(beta*I*S/N + gamma*I)
    }

    # the probability an infective recovers in this time step
    theta = 1/(gamma*k)
    a = 1-pgamma(vsojourn,shape=k,scale=theta)
    recover_prob = rep(1,N)
    l = which(a>1e-4)
    recover_prob[l] = (pgamma(vsojourn[l]+deltat,shape=k,scale=theta)-pgamma(vsojourn[l],shape=k,scale=theta))/a[l]

    # sample Poisson random numbers of infected people contacted by each person
    avg_num_infected_people_contacted = beta*I*deltat/N
    vnum_infected_people_contacted = rpois(N,avg_num_infected_people_contacted)
    vprob = runif(N) # sample uniform random numbers
    vnewstate = vstate # copy the state vector to a temporary vector used for calculations

    vsojourn = vsojourn+deltat
    # Infected people recover if the sampled uniform random number is less than the recovery probability
    vnewstate[vstate==1&vprob<recover_prob] = 2
    # reset the time they have spent in their new state to zero!
    vsojourn[vstate==1&vprob<recover_prob] = 0

    # If a susceptible contacted at least one infective, they are infected
    vnewstate[vstate==0&vnum_infected_people_contacted>0] = 1
    # reset the time they have spent in their new state to zero!
    vsojourn[vstate==0&vnum_infected_people_contacted>0] = 0

    vstate = vnewstate # update the state vector
    t = t + deltat # update the time
    }
    final = 0
    if (length(vS)>0) final = 1-min(vS)/N

    return(list(time=vtime,I=vI,S=vS,final_size=final))

    }
    output$distPlot <- renderPlot({
    set.seed(777)
    library(“deSolve”)
    nrealisations = input$bins3
    ################################
    # this is a function which, given a value of S,I and R at time t
    # calculates the time derivatives of S I and R
    # vparameters contains the parameters of the model, like the
    # recovery period, gamma, and the transmission rate, beta
    # this function gets passed to the deSolve package
    ################################
    SIRfunc=function(t, x, vparameters){
    S = x[1]
    I = x[2]
    R = x[3]
    if (I<0) I=0

    with(as.list(vparameters),{
    npop = S+I+R
    dS = -beta*S*I/npop
    dI = +beta*S*I/npop – gamma*I
    dR = +gamma*I
    out = c(dS,dI,dR)
    list(out)
    })
    }
    #############################
    # Set up initial conditions
    N = input$bins # population size
    I_0 = input$bins2 # number intially infected people in the population
    S_0 = N-I_0
    R_0 = 0 # assume no one has recovered at first

    delta_t = 0.1 # nominal time step
    tbeg = 0 # begin day
    tend = input$bins5 # end day
    gamma = input$bins4 # recovery period of influenza in days^{-1}
    R0 = input$bins6 # R0 of a hypothetical strain of pandemic influenza
    beta = R0*gamma # “reverse engineer” beta from R0 and gamma
    # first simulate the model with deterministic ODE’s, so that we have something
    # to compare our stochastic simulation to.
    ################################
    vt = seq(tbeg,tend,delta_t)
    vparameters = c(gamma=gamma,beta=beta)
    inits = c(S=S_0,I=I_0,R=R_0)
    sirmodel = as.data.frame(lsoda(inits, vt, SIRfunc, vparameters))
    # now plot the results of the deterministic model
    ################################
    par(mfrow=c(1,1)) # divides the page into two plotting areas

    plot(sirmodel$time,sirmodel$I/N,

    ylim=c(0,1.5*max(sirmodel$I/N)),

    type=”l”,col=1,lwd=5,xlab=”Tiempo (días)”,

    ylab=”Fracción infectada (prevalecencia)”,

    main=paste(“Pandemia de Influenza en una población de “,N,sep=””))
    cat(“The final size of epidemic from the deterministic model is “,max(sirmodel$R/N),”\n”)

    # now do several simulations using the agent based model, and overlay the
    # results on those from the deterministic model
    ################################
    vfinal = numeric(0) # we will fill this vector with the epidemic final size estimates from the simulations
    for (iter in 1:nrealisations){
    myagent = SIR_agent(N,I_0,S_0,gamma,R0,tbeg,tend,delta_t)
    lines(myagent$time,myagent$I/N,lwd=2,col=(iter+1),lty=3)
    cat(iter,nrealisations,”The final size of the epidemic from the agent based stochastic model is “,myagent$final_size,”\n”)
    vfinal = append(vfinal,myagent$final_size)
    }
    lines(sirmodel$time,sirmodel$I/N,lwd=5)
    cat(“The final size of epidemic from the deterministic model is “,max(sirmodel$R/N),”\n”)
    legend(“topright”,legend=c(“Determinismo”,”Simulaciones basadas en agentes”),col=c(1,2),lwd=3,lty=c(1,3),bty=”n”)

    })
    }

shinyApp(ui = ui, server = server)

