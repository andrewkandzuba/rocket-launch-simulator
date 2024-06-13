library(deSolve)
options(scipen=999)

# Define the differential equations
variable_mass_rocket <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    x <- state[1]  # x(t)
    vx <- state[2]  # dx/dt
    
    y <- state[3]  # y(t)
    vy <- state[4]  # dy/dt
    
    g = G * Me / (R + y)^2
    
    m0 = 0
    k  = 0
    
    F0 = 0
    mu = 0
    
    a0 = 90
    b  = 0
    
    # Find X(y)
    X_t = 0
    
    for(i in 1:nrow(t_gradients)) {
      tg = t_gradients[i,]
      h0 = as.double(tg["h0"])
      hk = as.double(tg["hk"])
      T0 = as.double(tg["T0"])
      dL = as.double(tg["dL"])
      P0 = as.double(tg["P0"])
      if(y >= h0 && y < hk ) {
        dHkm = (y - h0) * 0.001
        T_i = T0 + dL * dHkm

        if(dL != 0) {
          Pr_i = P0 * (T0/T_i) ^ (magic/dL)
        } else {
          Pr_i = P0 * exp(-magic * dHkm/T0)
        }
        
        rho_i = (Pr_i * M_)/(R * T_i)
        
        v_i = sqrt(vx^2 + vy^2)
        sound_vi =  (-1) * 0.0007 * T_i^2 + 0.9690 * T_i + 115.7923;
        mach_vi = v_i / sound_vi
        cx_i = Cx(mach_vi)
        
        X_t = cx_i * rho_i * (v_i^2) * S / 2
        
        atmosphere <<- rbind(atmosphere, list(
          t = t,
          h = y,
          v = v_i,
          T = T_i,
          Pr = Pr_i,
          rho = rho_i,
          sound_vi = sound_vi,
          mach_v = mach_vi,
          cx = cx_i,
          X = X_t,
          g = g
        ))

        break
      }
    }
    
    # Find P(t)

    t_t = 0
    
    stages <- data.frame(
      m0 = parameters["s.m0"],
      k =  parameters["s.k"],
      F0 = parameters["s.F0"],
      mu = parameters["s.mu"],
      Ts = parameters["s.Ts"],
      Te = parameters["s.Te"],
      b =  parameters["s.b"]
    )
    
    for(i in 1:nrow(stages)) {
      stage <- stages[i,]
      Ts = as.numeric(stage["s.Ts"])
      Te = as.numeric(stage["s.Te"])
      if (t >= Ts && t <= Te) {
        t_t = t - Ts
        m0  = as.double(stage["s.m0"])
        k   = as.double(stage["s.k"])
        F0  = as.double(stage["s.F0"])
        mu  = as.double(stage["s.mu"])
        b   = as.double(stage["s.b"])
      } 
    }

    a0 = tangent[nrow(tangent),]
    a = a0 - b
    tangent <<- rbind(tangent, a0 = c(a))

    F_t = F0 + mu * t_t
    m_t = m0 - k * t_t
    
    dx1 <- vx
    dy1 <- vy
    
    dx2 <- ((F_t - X_t) / m_t) * cos(a * (pi / (180)))
    dy2 <- ((F_t - X_t) / m_t) * sin(a * (pi / (180))) - g
    
    
    list(c(dx1, dx2, dy1, dy2))
  })
}

Cx <- function(m) {
  cx_table <- data.frame(
    M   =  c(0.00,	0.25,	0.50,	0.75,	1.0,	1.25,	1.50,	2.00,	3.00,	4.00,	5.00,	8.00,	9.00,	10.00,	11.00,	12.00,	13.00,	14.00,	15.00),
    Cx	 = c(0.78,	0.78,	0.78,	0.80,	1.1,	0.90,	0.85,	0.80,	0.75,	0.72,	0.70,	0.66,	0.66,	 0.66,	 0.66,	 0.66,	 0.66,	 0.66,	0.66)
  )
  
  k <- data.frame();
  for(i in 1:18) {
    k <- rbind(k, list(as.double((cx_table[i+1,]$Cx- cx_table[i,]$Cx)/(cx_table[i+1,]$M - cx_table[i,]$M))))
  }
  
  cx = 0.0;
  for(i in 1:nrow(k)) {
    if (cx_table[i,]$M <= m && m < cx_table[i+1,]$M) {
      cx = cx_table[i,]$Cx + as.double(k[i,]) * (m - cx_table[i,]$M)
    }
  }
  
  cx
}

# Atmospheric temperature
t_gradients = data.frame(
  h0 = c(     0,  11000,   20000,    32000,   47000,   51000,  71000),
  hk = c( 11000,  20000,   32000,    47000,   51000,   71000,  85000),
  T0 = c(288.15, 216.15,  216.15,   228.65,  270.00,  270.00, 216.00),
  dL = c( -6.50,  0.00,     1.00,     2.80,    0.00,   -2.80,  -2.00),
  P0 = c(1030.0,  229.8,    55.3,      8.7,     1.1,     0.6,   0.03)
)

g0 =    9.81
magic = 34.163
M_ =    2.89644
R =     8.31447
S =     3.141593 * (3.05/2)^2 ## Atlas-D Middle Squre
G =     6.67430 * 10^-11
Me =    5.972 * 10^24 
R =     6371000.0


TBurn = 315

# atmosphere model
atmosphere = data.frame()

# Tangent 
tangent = data.frame(a0 = c(90))

# Initial conditions: x(0) = 0, dx/dt(0) = 0
state <- c(x = 0, vx = 0, y = 0, vy = 0)

# Rocket stages
#           stage0            stage1,         stage2      
stages <- data.frame(
  m0  = c(   117500.00,      111479.89,   24817.23 ),
  k   = c(      704.7,           704.7,     120.00 ),
  F0  = c(  1587100.00,     1595354.69,  299317.19 ),
  mu  = c(      917.19,         917.19,     345.31 ),
  Ts  = c(          0,              11,        136 ),
  Te  = c(         10,             135,      TBurn ),
  b   = c(          0,             0.6,       0.11 )
)

# Parameters: F (thrust force), m0 (initial mass), k (mass loss rate), g (gravity)
parameters <- c(s = stages)

# Time sequence
times <- seq(0, TBurn, by = 1)

# Solve the differential equations
solution <- ode(y = state, times = times, func = variable_mass_rocket, parms = parameters, method="euler")


# Adjust the height by the Earth curve

for(i in 1:nrow(solution)) {
  row <- solution[i,]
  x = as.double(row["x"])
  y = as.double(row["y"])
  h = R/cos(atan(x/R)) - R
  row["y"] = y + h
  row["ys"] = -h
}

# Plot the results
plot(solution, xlab = "Time", ylab = "State Variables")
legend("right", legend = c("x(t)", "dx/dt(t)"), col = 1:2, lty = 1:2)

# Plot the trajectory 
plot(solution[,c("x","y")], xlab = "Distance (m)", ylab="Height (m)", type = "l")

# Plot atmosphere
plot(atmosphere[,c("T","h")], ylim = c(0, 85000), xlab = "Temperature (K)", ylab="Height (m)", type = "l", )
plot(atmosphere[,c("h", "Pr")], ylab = "Pressure (kPa)", xlab="Height (m)", type = "l")
plot(atmosphere[,c("h", "rho")], ylab = "Density (K)", xlab="Height (m)", type = "l")
plot(atmosphere[,c("h", "mach_v")], xlab = "Height (m)", ylab="Mach", type = "l")
plot(atmosphere[,c("t", "mach_v")], xlab = "Time (s)", ylab="Mach", type = "l")
plot(atmosphere[,c("mach_v", "cx")], xlab = "Mach", ylab="Cx", type = "l")
plot(atmosphere[,c("t", "X")], xlab = "Time (s)", ylab="Aerodynamic Force", type = "l")
