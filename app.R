## app.R — ARMA Simulation Lab
## Requisitos: shiny, stats, graphics, withr

library(shiny)
library(withr)

# ---------- 0) Utilidades ----------

# Tu función de correlograma (tal cual la pasaste)
correlograma <- function(
    x,
    max_lag = 36,
    orientation = c("vertical","horizontal"),
    cex = 0.8
){
  orientation <- match.arg(orientation)
  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE)
  par(
    las = 1,
    cex.axis = cex,
    cex.lab  = cex,
    cex.main = cex,
    mar = c(2.2, 3.2, 0.6, 0.6),
    oma = c(0.6, 0, 0, 0),
    mgp = c(1.4, 0.4, 0),
    tcl = -0.2
  )
  if (orientation == "horizontal") {
    par(mfrow = c(1, 2))
  } else {
    par(mfrow = c(2, 1))
  }
  acf(x, main = "", ylab = "ACF",  xlab = "Rezagos", lag.max = max_lag)
  pacf(x, main = "", ylab = "PACF", xlab = "Rezagos", lag.max = max_lag)
  par(mfrow = c(1, 1))
}

# Construye el vector de coeficientes del polinomio AR no estacional: 1 - phi1 z - ... - phip z^p
poly_coefs_ar <- function(phi) c(1, -as.numeric(phi))

# Polinomio MA no estacional: 1 + theta1 z + ... + thetaq z^q
poly_coefs_ma <- function(theta) c(1, as.numeric(theta))

# Polinomio estacional AR: 1 - Phi1 z^s - Phi2 z^{2s} - ...
poly_coefs_seasonal_ar <- function(Phi, s) {
  P <- length(Phi)
  if (P == 0) return(1)
  coefs <- rep(0, P * s + 1)
  coefs[1] <- 1
  for (i in seq_len(P)) coefs[i * s + 1] <- -Phi[i]
  coefs
}

# Polinomio estacional MA: 1 + Theta1 z^s + Theta2 z^{2s} + ...
poly_coefs_seasonal_ma <- function(Theta, s) {
  Q <- length(Theta)
  if (Q == 0) return(1)
  coefs <- rep(0, Q * s + 1)
  coefs[1] <- 1
  for (i in seq_len(Q)) coefs[i * s + 1] <- Theta[i]
  coefs
}

# Convolución de polinomios (producto)
poly_convolve <- function(a, b) as.numeric(convolve(a, rev(b), type = "open"))

# Raíces (en z) y módulo
poly_roots <- function(coefs) {
  # polyroot espera coeficientes en orden ascendente: c(a0, a1, a2, ...)
  roots <- polyroot(coefs)
  data.frame(root = roots, modulus = Mod(roots))
}

# Chequeo de estacionariedad/invertibilidad: todas las raíces del polinomio combinado fuera del círculo unidad
check_stability_invertibility <- function(phi, Phi, theta, Theta, s) {
  # Polinomios combinados = no estacional * estacional
  AR_poly <- poly_convolve(poly_coefs_ar(phi), poly_coefs_seasonal_ar(Phi, s))
  MA_poly <- poly_convolve(poly_coefs_ma(theta), poly_coefs_seasonal_ma(Theta, s))
  AR_roots <- poly_roots(AR_poly)
  MA_roots <- poly_roots(MA_poly)
  list(
    AR = list(coefs = AR_poly, roots = AR_roots, ok = all(AR_roots$modulus > 1)),
    MA = list(coefs = MA_poly, roots = MA_roots, ok = all(MA_roots$modulus > 1))
  )
}

# Simulación ARIMA estacionario (d=D=0); agregamos "nivel" como desplazamiento
simulate_series <- function(n, phi, theta, Phi, Theta, s, sigma2, level = 0) {
  model <- list(ar = phi, ma = theta, order = c(length(phi), 0, length(theta)),
                seasonal = list(order = c(length(Phi), 0, length(Theta)), period = s))
  y <- arima.sim(n = n, model = model, sd = sqrt(sigma2))
  as.numeric(y) + level
}

# Formateo amigable de un modelo
fmt_model <- function(phi, theta, Phi, Theta, s) {
  paste0(
    "ARMA(", length(phi), ",", length(theta), ")",
    if (length(Phi) + length(Theta) > 0) paste0(" x SARMA(", length(Phi), ",", length(Theta), ")[", s, "]") else "",
    if (length(phi) > 0) paste0("\n  phi = [", paste(round(phi, 3), collapse = ", "), "]") else "",
    if (length(theta) > 0) paste0("\n  theta = [", paste(round(theta, 3), collapse = ", "), "]") else "",
    if (length(Phi) > 0) paste0("\n  Phi = [", paste(round(Phi, 3), collapse = ", "), "]") else "",
    if (length(Theta) > 0) paste0("\n  Theta = [", paste(round(Theta, 3), collapse = ", "), "]") else ""
  )
}

# Genera coeficientes estables/invertibles a partir de raíces fuera del círculo unidad
coeffs_from_stable_roots <- function(p, q, P, Q, s,
                                     ar_mod_min = 1.2, ma_mod_min = 1.2) {
  # No estacional AR
  phi <- numeric(0)
  if (p > 0) {
    r_mod <- runif(p, ar_mod_min, ar_mod_min + 0.8)
    r_ang <- runif(p, -pi, pi)
    r     <- r_mod * exp(1i * r_ang)
    # Polinomio con raíces r: prod(1 - r z) = 1 + c1 z + c2 z^2 + ...
    # phi_i = -c_i  (porque 1 - phi1 z - ... = prod(1 - r z))
    c <- 1
    for (ri in r) c <- poly_convolve(c, c(1, -ri))
    c <- Re(c) # pequeños imaginarios numéricos
    phi <- -c[-1]
  }
  # No estacional MA
  theta <- numeric(0)
  if (q > 0) {
    r_mod <- runif(q, ma_mod_min, ma_mod_min + 0.8)
    r_ang <- runif(q, -pi, pi)
    r     <- r_mod * exp(1i * r_ang)
    # 1 + theta1 z + ... = prod(1 + r z) => theta_i = coeficientes sin el término constante
    c <- 1
    for (ri in r) c <- poly_convolve(c, c(1, ri))
    c <- Re(c)
    theta <- c[-1]
  }
  # Estacional AR (raíces para z^s)
  Phi <- numeric(0)
  if (P > 0) {
    r_mod <- runif(P, ar_mod_min, ar_mod_min + 0.8)
    r_ang <- runif(P, -pi, pi)
    r     <- r_mod * exp(1i * r_ang)
    # prod(1 - r z^s) -> coef en z^{ks} con -Phi_k
    # Los coeficientes Phi_k salen del polinomio en variable w = z^s, igual lógica que AR
    c <- 1
    for (ri in r) c <- poly_convolve(c, c(1, -ri))
    c <- Re(c)
    Phi <- -c[-1]
  }
  # Estacional MA
  Theta <- numeric(0)
  if (Q > 0) {
    r_mod <- runif(Q, ma_mod_min, ma_mod_min + 0.8)
    r_ang <- runif(Q, -pi, pi)
    r     <- r_mod * exp(1i * r_ang)
    c <- 1
    for (ri in r) c <- poly_convolve(c, c(1, ri))
    c <- Re(c)
    Theta <- c[-1]
  }
  list(phi = as.numeric(phi), theta = as.numeric(theta),
       Phi = as.numeric(Phi), Theta = as.numeric(Theta))
}

# Forecast under Quadratic and LINEX losses given mean and s.e.
linex_point <- function(mean, se, a) {
  # Si Y|info ~ N(mean, se^2), el minimizador de E[exp(a e) - a e - 1] es mean - 0.5*a*se^2
  mean - 0.5 * a * (se^2)
}

# Ajustar un ARIMA con órdenes dados (sin diferenciar)
fit_arima_safe <- function(y, p, q, P, Q, s) {
  stats::arima(y,
               order = c(p, 0, q),
               seasonal = list(order = c(P, 0, Q), period = s),
               include.mean = TRUE, method = "ML")
}

# ---------- 1) UI ----------
ui <- fluidPage(
  titlePanel("ARMA Simulation Lab"),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Control de semilla"),
      numericInput("seed", "Semilla (entero):", value = 1234, min = 0, step = 1),
      actionButton("set_seed", "Setear semilla"),
      tags$hr(),
      selectInput("freq", "Frecuencia:", choices = c("Mensual (12)" = 12, "Trimestral (4)" = 4), selected = 12),
      uiOutput("sidebar_tab_controls")
      # tabsetPanel(
      #   id = "tabs",
      #   type = "pills",
      #   tabPanel("Simular modelo",
      #            br(),
      #            numericInput("n1", "Tamaño de muestra (n):", value = 300, min = 50, step = 10),
      #            numericInput("mu1", "Nivel (desplazamiento):", value = 0, step = 0.1),
      #            numericInput("sigma2_1", "Varianza del residuo:", value = 1, min = 0.0001, step = 0.1),
      #            tags$hr(),
      #            sliderInput("p1", "Orden AR (p):", min = 0, max = 6, value = 1, step = 1),
      #            uiOutput("phi_inputs"),
      #            sliderInput("q1", "Orden MA (q):", min = 0, max = 6, value = 1, step = 1),
      #            uiOutput("theta_inputs"),
      #            sliderInput("P1", "Orden SAR (P):", min = 0, max = 2, value = 0, step = 1),
      #            uiOutput("Phi_inputs"),
      #            sliderInput("Q1", "Orden SMA (Q):", min = 0, max = 2, value = 0, step = 1),
      #            uiOutput("Theta_inputs"),
      #            actionButton("simulate1", "Simular")
      #   ),
      #   tabPanel("Adivinar el modelo",
      #            br(),
      #            numericInput("n2", "Tamaño de muestra (n):", value = 300, min = 50, step = 10),
      #            selectInput("level", "Nivel de dificultad:", choices = c("Fácil", "Medio", "Difícil")),
      #            actionButton("simulate2", "Simular (oculto)"),
      #            tags$hr(),
      #            h5("Proponer modelos a estimar"),
      #            fluidRow(
      #              column(6, sliderInput("p2", "p", min = 0, max = 6, value = 1, step = 1)),
      #              column(6, sliderInput("q2", "q", min = 0, max = 6, value = 1, step = 1))
      #            ),
      #            fluidRow(
      #              column(6, sliderInput("P2", "P (<=2)", min = 0, max = 2, value = 0, step = 1)),
      #              column(6, sliderInput("Q2", "Q (<=2)", min = 0, max = 2, value = 0, step = 1))
      #            ),
      #            actionButton("add_model", "Agregar propuesta"),
      #            actionButton("reveal", "Revelar verdadero modelo"),
      #            tags$hr(),
      #            tableOutput("proposals_table")
      #   ),
      #   tabPanel("Predicciones",
      #            br(),
      #            numericInput("n3", "Tamaño de muestra (n):", value = 300, min = 50, step = 10),
      #            numericInput("mu3", "Nivel (desplazamiento):", value = 0, step = 0.1),
      #            numericInput("sigma2_3", "Varianza del residuo:", value = 1, min = 0.0001, step = 0.1),
      #            sliderInput("p3", "p", min = 0, max = 6, value = 1, step = 1),
      #            uiOutput("phi3_inputs"),
      #            sliderInput("q3", "q", min = 0, max = 6, value = 1, step = 1),
      #            uiOutput("theta3_inputs"),
      #            sliderInput("P3", "P", min = 0, max = 2, value = 0, step = 1),
      #            uiOutput("Phi3_inputs"),
      #            sliderInput("Q3", "Q", min = 0, max = 2, value = 0, step = 1),
      #            uiOutput("Theta3_inputs"),
      #            numericInput("h", "Horizonte de proyección (h):", value = 8, min = 1, step = 1),
      #            checkboxGroupInput("losses", "Funciones de pérdida a graficar:",
      #                               choices = c("Cuadrática" = "quad", "LINEX" = "linex"),
      #                               selected = c("quad")),
      #            conditionalPanel(
      #              condition = "input.losses.includes('linex')",
      #              numericInput("a_linex", "Parámetro 'a' de LINEX:", value = 0.5, step = 0.1)
      #            ),
      #            actionButton("simulate3", "Simular y predecir")
      #   )
      # )
    ),
    
    
    mainPanel(
      width = 9,
      tabsetPanel(
        id = "outtabs",
        type = "tabs",
        tabPanel("Introducción",
                 tags$div(style = "max-width: 900px;",
                          h3("ARMA/ARIMA Simulation Lab — Piloto 2025"),
                          p(
                            "Esta aplicación es un piloto diseñado para acompañar el curso de ",
                            strong("Econometría de Series de Tiempo"),
                            " de la ",
                            strong("Universidad ORT Uruguay"),
                            " en el año ",
                            strong("2025"),
                            ". Su objetivo es ofrecer un entorno interactivo para experimentar con la",
                            " simulación, identificación y evaluación de modelos ARMA/SARMA y la generación de",
                            " predicciones bajo diferentes funciones de pérdida."
                          ),
                          h4("¿Qué podés hacer acá?"),
                          tags$ul(
                            tags$li(strong("Simular modelos:"),
                                    " elegí órdenes, coeficientes y varianza; observá la serie, su ACF/PACF y",
                                    " comprobá estacionariedad e invertibilidad vía raíces del polinomio."),
                            tags$li(strong("Adivinar el modelo:"),
                                    " genera una serie oculta (por niveles de dificultad), proponé órdenes y compará AIC/BIC; ",
                                    " luego revelá el modelo verdadero."),
                            tags$li(strong("Proyectar:"),
                                    " construí pronósticos multi-paso y compará la media condicional (pérdida cuadrática)",
                                    " con la predicción óptima bajo ",
                                    em("LINEX"),
                                    " (elegí el parámetro ",
                                    code("a"), ").")
                          ),
                          h4("Reproducibilidad"),
                          p(
                            "Podés fijar la semilla con el botón ",
                            code("Setear semilla"),
                            " para reproducir resultados. Si no fijás la semilla, cada simulación generará una realización distinta."
                          ),
                          h4("Sugerencia de uso en clase"),
                          tags$ul(
                            tags$li("Explorá primero “Simular modelo” para identificar patrones ACF/PACF."),
                            tags$li("Probá el desafío de “Adivinar el modelo” por niveles."),
                            tags$li("Cerrá con “Predicciones” comparando cuadrática vs. LINEX.")
                          ),
                          p(em("Versión piloto — comentarios y mejoras bienvenidos."))
                 )
        ),
        tabPanel("Simular modelo",
                 plotOutput("plot_series_1", height = 260),
                 plotOutput("acf_1", height = 220),
                 #fluidRow(
                #   column(6, plotOutput("acf_1", height = 220))
                   #column(6, plotOutput("pacf_1", height = 220))
                 #),
                 tags$hr(),
                 h4("Raíces y condición de estabilidad/invertibilidad"),
                 verbatimTextOutput("roots_box_1")
        ),
        tabPanel("Adivinar el modelo",
                 plotOutput("plot_series_2", height = 260),
                 plotOutput("acf_2", height = 220),
                 #fluidRow(
                #   column(6, plotOutput("acf_2", height = 220)),
                #   column(6, plotOutput("pacf_2", height = 220))
                 #),
                 tags$hr(),
                 h4("Modelos propuestos"),
                 tableOutput("proposals_table"), 
                 tags$hr(),
                 h4("Verdadero modelo"),
                 verbatimTextOutput("true_model_box")
        ),
        tabPanel("Predicciones",
                 plotOutput("plot_forecast", height = 320),
                 verbatimTextOutput("forecast_info")
        )
      )
    )
  )
)

# ---------- 2) SERVER ----------
server <- function(input, output, session) {
  
  # Frecuencia (s): 12 o 4
  s <- reactive(as.integer(input$freq))
  
  # Semilla controlada
  seed_val <- reactiveVal(1234)
  observeEvent(input$set_seed, {
    req(input$seed)
    seed_val(as.integer(input$seed))
    showNotification(paste("Semilla seteada a", seed_val()), type = "message")
  })
  
  output$sidebar_tab_controls <- renderUI({
    req(input$outtabs)
    switch(input$outtabs,
           "Simular modelo" = tagList(
             br(),
             numericInput("n1", "Tamaño de muestra (n):", value = 300, min = 50, step = 10),
             numericInput("mu1", "Nivel (desplazamiento):", value = 0, step = 0.1),
             numericInput("sigma2_1", "Varianza del residuo:", value = 1, min = 0.0001, step = 0.1),
             tags$hr(),
             sliderInput("p1", "Orden AR (p):", min = 0, max = 6, value = 1, step = 1),
             uiOutput("phi_inputs"),
             sliderInput("q1", "Orden MA (q):", min = 0, max = 6, value = 1, step = 1),
             uiOutput("theta_inputs"),
             sliderInput("P1", "Orden SAR (P):", min = 0, max = 2, value = 0, step = 1),
             uiOutput("Phi_inputs"),
             sliderInput("Q1", "Orden SMA (Q):", min = 0, max = 2, value = 0, step = 1),
             uiOutput("Theta_inputs"),
             actionButton("simulate1", "Simular")
           ),
           
           "Adivinar el modelo" = tagList(
             br(),
             numericInput("n2", "Tamaño de muestra (n):", value = 300, min = 50, step = 10),
             selectInput("level", "Nivel de dificultad:", choices = c("Fácil", "Medio", "Difícil")),
             actionButton("simulate2", "Simular (oculto)"),
             tags$hr(),
             h5("Proponer modelos a estimar"),
             fluidRow(
               column(6, sliderInput("p2", "p", min = 0, max = 6, value = 1, step = 1)),
               column(6, sliderInput("q2", "q", min = 0, max = 6, value = 1, step = 1))
             ),
             fluidRow(
               column(6, sliderInput("P2", "P (≤2)", min = 0, max = 2, value = 0, step = 1)),
               column(6, sliderInput("Q2", "Q (≤2)", min = 0, max = 2, value = 0, step = 1))
             ),
             actionButton("add_model", "Agregar propuesta"),
             actionButton("reveal", "Revelar verdadero modelo")
             # La tabla de propuestas queda en el mainPanel como ya la tenías
           ),
           
           "Predicciones" = tagList(
             br(),
             numericInput("n3", "Tamaño de muestra (n):", value = 300, min = 50, step = 10),
             numericInput("mu3", "Nivel (desplazamiento):", value = 0, step = 0.1),
             numericInput("sigma2_3", "Varianza del residuo:", value = 1, min = 0.0001, step = 0.1),
             sliderInput("p3", "p", min = 0, max = 6, value = 1, step = 1),
             uiOutput("phi3_inputs"),
             sliderInput("q3", "q", min = 0, max = 6, value = 1, step = 1),
             uiOutput("theta3_inputs"),
             sliderInput("P3", "P", min = 0, max = 2, value = 0, step = 1),
             uiOutput("Phi3_inputs"),
             sliderInput("Q3", "Q", min = 0, max = 2, value = 0, step = 1),
             uiOutput("Theta3_inputs"),
             numericInput("h", "Horizonte de proyección (h):", value = 8, min = 1, step = 1),
             checkboxGroupInput("losses", "Funciones de pérdida a graficar:",
                                choices = c("Cuadrática" = "quad", "LINEX" = "linex"),
                                selected = c("quad")),
             conditionalPanel(
               condition = "input.losses.includes('linex')",
               numericInput("a_linex", "Parámetro 'a' de LINEX:", value = 0.5, step = 0.1)
             ),
             actionButton("simulate3", "Simular y predecir")
           ),
           
           # default
           tagList(h5("Seleccioná una pestaña de salida para ver sus controles."))
    )
  })
  
  # --- UI dinámico para coeficientes (pestaña 1) ---
  output$phi_inputs <- renderUI({
    if (input$p1 > 0) {
      lapply(seq_len(input$p1), function(i)
        numericInput(paste0("phi1_", i), paste0("phi[", i, "]:"), value = ifelse(i==1, 0.5, 0), step = 0.05))
    }
  })
  output$theta_inputs <- renderUI({
    if (input$q1 > 0) {
      lapply(seq_len(input$q1), function(i)
        numericInput(paste0("theta1_", i), paste0("theta[", i, "]:"), value = ifelse(i==1, 0.3, 0), step = 0.05))
    }
  })
  output$Phi_inputs <- renderUI({
    if (input$P1 > 0) {
      lapply(seq_len(input$P1), function(i)
        numericInput(paste0("Phi1_", i), paste0("Phi[", i, "]:"), value = 0, step = 0.05))
    }
  })
  output$Theta_inputs <- renderUI({
    if (input$Q1 > 0) {
      lapply(seq_len(input$Q1), function(i)
        numericInput(paste0("Theta1_", i), paste0("Theta[", i, "]:"), value = 0, step = 0.05))
    }
  })
  
  # --- Simulación pestaña 1 ---
  sim1 <- eventReactive(input$simulate1, {
    #with_seed(seed_val(),
              {
      p <- input$p1; q <- input$q1; P <- input$P1; Q <- input$Q1
      phi   <- if (p>0) sapply(seq_len(p), function(i) input[[paste0("phi1_", i)]]) else numeric(0)
      theta <- if (q>0) sapply(seq_len(q), function(i) input[[paste0("theta1_", i)]]) else numeric(0)
      Phi   <- if (P>0) sapply(seq_len(P), function(i) input[[paste0("Phi1_", i)]]) else numeric(0)
      Theta <- if (Q>0) sapply(seq_len(Q), function(i) input[[paste0("Theta1_", i)]]) else numeric(0)
      chk <- check_stability_invertibility(phi, Phi, theta, Theta, s())
      y <- simulate_series(input$n1, phi, theta, Phi, Theta, s(),
                           sigma2 = input$sigma2_1, level = input$mu1)
      list(y = y, phi = phi, theta = theta, Phi = Phi, Theta = Theta, chk = chk)
    }
    #)
  })
  
  output$plot_series_1 <- renderPlot({
    req(sim1())
    plot(sim1()$y, type = "l", main = "Serie simulada", xlab = "t", ylab = "y_t")
  })
  output$acf_1 <- renderPlot({
    req(sim1())
    correlograma(sim1()$y, max_lag = 36, orientation = "horizontal")
  })
  # Reemplazar este bloque en app.R
  #output$pacf_1 <- renderPlot({
  #  req(sim1())
  #  y <- sim1()$y
  #  n <- length(y)
  #  m <- mean(y)
  #  s <- sd(y)
  #  se_mean <- s / sqrt(n)
  #  ci <- m + c(-1, 1) * 1.96 * se_mean
    
  #  plot(y, type = "l", xlab = "t", ylab = "y_t",
  #       main = "Serie + IC 95% de la media muestral")
  #  abline(h = m, lwd = 2)           # media muestral
  #  abline(h = ci, lty = 3)          # IC 95% de la media
  #  legend("topleft",
  #         legend = c("Serie", "Media muestral", "IC 95% (media)"),
  #         lty = c(1, 1, 3), lwd = c(1, 2, 1), bty = "n")
  #})
  output$roots_box_1 <- renderPrint({
    req(sim1())
    cat("Polinomio AR combinado (no estacional * estacional):\n")
    print(round(sim1()$chk$AR$roots, 4))
    cat("\n¿Estacionario? ->", ifelse(sim1()$chk$AR$ok, "SÍ", "NO"), "\n\n")
    cat("Polinomio MA combinado (no estacional * estacional):\n")
    print(round(sim1()$chk$MA$roots, 4))
    cat("\n¿Invertible? ->", ifelse(sim1()$chk$MA$ok, "SÍ", "NO"), "\n")
  })
  
  # --- Pestaña 2: Adivinar el modelo ---
  hidden_model <- reactiveVal(NULL)
  sim2_series  <- reactiveVal(NULL)
  revealed     <- reactiveVal(FALSE) 
  
  observeEvent(input$simulate2, {
    #with_seed(seed_val(), 
    {
      # 1) Resetear estado previo
      revealed(FALSE)  # ocultar "Verdadero modelo"
      proposals(data.frame(p=integer(), q=integer(), P=integer(), Q=integer(),
                           AIC=numeric(), BIC=numeric(), stringsAsFactors = FALSE))
      hidden_model(NULL)
      sim2_series(NULL)
      
      # 2) Selección de órdenes por nivel (igual que antes)
      lvl <- input$level
      s_  <- s()
      n   <- input$n2
      
      pick_orders <- function() {
        if (lvl == "Fácil") {
          choice <- sample(c("AR","MA","SAR","SMA"), size = 1)
          p <- q <- P <- Q <- 0
          if (choice == "AR")  p <- sample(1:2, 1)
          if (choice == "MA")  q <- sample(1:2, 1)
          if (choice == "SAR") P <- sample(1:2, 1)
          if (choice == "SMA") Q <- sample(1:2, 1)
        } else if (lvl == "Medio") {
          p <- sample(0:3, 1); q <- sample(0:3, 1)
          if (p == 0 && q == 0) p <- 1
          P <- Q <- 0
        } else {
          p <- sample(0:3, 1); q <- sample(0:3, 1)
          P <- sample(0:2, 1); Q <- sample(0:2, 1)
          if (p + q + P + Q == 0) p <- 1
        }
        list(p=p,q=q,P=P,Q=Q)
      }
      
      # 3) Generar modelo válido (estable/invertible)
      repeat {
        ord <- pick_orders()
        co  <- coeffs_from_stable_roots(ord$p, ord$q, ord$P, ord$Q, s_)
        chk <- check_stability_invertibility(co$phi, co$Phi, co$theta, co$Theta, s_)
        if (chk$AR$ok && chk$MA$ok) {
          hidden_model(list(
            phi = co$phi, theta = co$theta, Phi = co$Phi, Theta = co$Theta,
            p = ord$p, q = ord$q, P = ord$P, Q = ord$Q, s = s_
          ))
          break
        }
      }
      
      # 4) Simular nueva serie
      y <- simulate_series(n = n,
                           phi = hidden_model()$phi,
                           theta = hidden_model()$theta,
                           Phi = hidden_model()$Phi,
                           Theta = hidden_model()$Theta,
                           s = s_,
                           sigma2 = 1, level = 0)
      sim2_series(y)
      
      showNotification("Nueva serie simulada. Proponé modelos y luego revelá el verdadero.", type = "message")
    }
    #)
  })
  
  # Tabla de propuestas
  proposals <- reactiveVal(data.frame(
    p = integer(), q = integer(), P = integer(), Q = integer(),
    AIC = numeric(), BIC = numeric(), stringsAsFactors = FALSE
  ))
  
  observeEvent(input$add_model, {
    req(sim2_series())            # asegurate de haber simulado la serie oculta
    req(!is.null(hidden_model())) # por si alguien pulsa antes de tiempo
    
    p <- input$p2; q <- input$q2; P <- input$P2; Q <- input$Q2
    
    fit <- try(fit_arima_safe(sim2_series(), p, q, P, Q, s()), silent = TRUE)
    if (inherits(fit, "try-error")) {
      showNotification("Estimación falló (probá otros órdenes).", type = "error")
      return(invisible(NULL))
    }
    
    # actualizar tabla
    df <- proposals()
    df <- rbind(df, data.frame(
      p = p, q = q, P = P, Q = Q,
      AIC = AIC(fit),
      BIC = AIC(fit, k = log(length(sim2_series())))
    ))
    proposals(df)
    
    showNotification(sprintf("Agregada propuesta: (%d,%d) x (%d,%d)", p, q, P, Q), type = "message")
  }, ignoreInit = TRUE)
  
  output$proposals_table <- renderTable({
    df <- proposals()
    if (is.null(df) || nrow(df) == 0) return(data.frame(Mensaje="Sin propuestas todavía"))
    df
  }, digits = 3)
  
  output$plot_series_2 <- renderPlot({
    req(sim2_series())
    plot(sim2_series(), type = "l", main = "Serie simulada (oculta)", xlab = "t", ylab = "y_t")
  })
  output$acf_2 <- renderPlot({
    req(sim2_series())
    correlograma(sim2_series(), max_lag = 36, orientation = "horizontal")
  })
  #output$pacf_2 <- renderPlot({
  #  req(sim2_series())
  #  pacf(sim2_series(), main = "PACF", xlab = "Rezagos", lag.max = 36)
  #})
  
  # Calcula la info del verdadero modelo cuando se presiona "Revelar"
  output_true <- eventReactive(input$reveal, {
    req(hidden_model())
    m <- hidden_model()
    list(
      desc = fmt_model(m$phi, m$theta, m$Phi, m$Theta, m$s),
      chk  = check_stability_invertibility(m$phi, m$Phi, m$theta, m$Theta, m$s)
    )
  })
  
  # Cambia el estado de revelado cuando se hace click
  observeEvent(input$reveal, {
    revealed(TRUE)
  })
  
  # Mostrar/ocultar el verdadero modelo según 'revealed()'
  output$true_model_box <- renderPrint({
    if (!revealed()) {
      cat("Presioná 'Revelar verdadero modelo' para mostrarlo.")
      return(invisible(NULL))
    }
    info <- output_true()  # esto ya corre porque revealed() es TRUE
    cat(info$desc, "\n\n")
    cat("Raíces AR (combinadas):\n")
    print(round(info$chk$AR$roots, 4))
    cat("\nRaíces MA (combinadas):\n")
    print(round(info$chk$MA$roots, 4))
    cat("\nEstacionario:", ifelse(info$chk$AR$ok, "SÍ", "NO"),
        " | Invertible:", ifelse(info$chk$MA$ok, "SÍ", "NO"), "\n")
  })
  
  # --- Pestaña 3: Predicciones + pérdidas ---
  # UI dinámico coeficientes
  output$phi3_inputs <- renderUI({
    if (input$p3 > 0) {
      lapply(seq_len(input$p3), function(i)
        numericInput(paste0("phi3_", i), paste0("phi[", i, "]:"), value = ifelse(i==1, 0.5, 0), step = 0.05))
    }
  })
  output$theta3_inputs <- renderUI({
    if (input$q3 > 0) {
      lapply(seq_len(input$q3), function(i)
        numericInput(paste0("theta3_", i), paste0("theta[", i, "]:"), value = ifelse(i==1, 0.3, 0), step = 0.05))
    }
  })
  output$Phi3_inputs <- renderUI({
    if (input$P3 > 0) {
      lapply(seq_len(input$P3), function(i)
        numericInput(paste0("Phi3_", i), paste0("Phi[", i, "]:"), value = 0, step = 0.05))
    }
  })
  output$Theta3_inputs <- renderUI({
    if (input$Q3 > 0) {
      lapply(seq_len(input$Q3), function(i)
        numericInput(paste0("Theta3_", i), paste0("Theta[", i, "]:"), value = 0, step = 0.05))
    }
  })
  
  pred3 <- eventReactive(input$simulate3, {
    #with_seed(seed_val(), 
              {
      p <- input$p3; q <- input$q3; P <- input$P3; Q <- input$Q3
      phi   <- if (p>0) sapply(seq_len(p), function(i) input[[paste0("phi3_", i)]]) else numeric(0)
      theta <- if (q>0) sapply(seq_len(q), function(i) input[[paste0("theta3_", i)]]) else numeric(0)
      Phi   <- if (P>0) sapply(seq_len(P), function(i) input[[paste0("Phi3_", i)]]) else numeric(0)
      Theta <- if (Q>0) sapply(seq_len(Q), function(i) input[[paste0("Theta3_", i)]]) else numeric(0)
      
      chk <- check_stability_invertibility(phi, Phi, theta, Theta, s())
      if (!chk$AR$ok || !chk$MA$ok) {
        showNotification("El modelo elegido no es estacionario/invertible. Ajustá coeficientes.", type = "error")
        return(NULL)
      }
      
      y <- simulate_series(input$n3, phi, theta, Phi, Theta, s(),
                           sigma2 = input$sigma2_3, level = input$mu3)
      
      fit <- fit_arima_safe(y, p, q, P, Q, s())
      pr  <- predict(fit, n.ahead = input$h)
      mean_fc <- as.numeric(pr$pred)
      se_fc   <- as.numeric(pr$se)
      
      list(y = y, fit = fit, mean_fc = mean_fc, se_fc = se_fc,
           h = input$h, losses = input$losses, a = input$a_linex)
    }
    #)
  })
  
  output$plot_forecast <- renderPlot({
    req(pred3())
    y  <- pred3()$y
    h  <- pred3()$h
    mu <- pred3()$mean_fc
    se <- pred3()$se_fc
    losses <- pred3()$losses
    a <- pred3()$a
    
    # últimos 24 datos
    idx_start <- max(1, length(y) - 36 + 1)
    y_last <- y[idx_start:length(y)]
    x0 <- length(y_last)
    xf <- x0 + h
    x_fc <- (x0 + 1):xf
    
    # Plot base extendiendo eje x
    plot(y_last, type = "l",
         xlab = "t", ylab = "y_t",
         main = "Últimos 36 + proyección",
         xlim = c(1, xf))
    abline(v = x0, lty = 2)
    
    leg <- character(0); lty <- integer(0); lwd <- integer(0); colv <- character(0)
    
    if ("quad" %in% losses) {
      upper <- mu + 1.96 * se
      lower <- mu - 1.96 * se
      # Banda sombreada para el IC (gris claro, sin borde)
      polygon(
        x = c(x_fc, rev(x_fc)),
        y = c(upper, rev(lower)),
        col = rgb(0.5, 0.5, 0.5, 0.25), border = NA
      )
      # Línea de la media condicional (negra, dashed)
      lines(x_fc, mu, lwd = 2, lty = 2, col = "black")
      
      leg <- c(leg, "Media cond. (Cuadrática)", "IC 95%")
      lty <- c(lty, 2, NA)     # NA para el parche
      lwd <- c(lwd, 2, NA)
      colv <- c(colv, "black", rgb(0.5,0.5,0.5,0.25))
    }
    
    if ("linex" %in% losses) {
      mu_linex <- mu - 0.5 * a * (se^2)
      lines(x_fc, mu_linex, lwd = 2, lty = 1, col = "blue")
      
      leg <- c(leg, paste0("Pto. LINEX (a=", a, ")"))
      lty <- c(lty, 1)
      lwd <- c(lwd, 2)
      colv <- c(colv, "blue")
    }
    
    # Leyenda (el ítem de banda queda sin línea)
    legend("topleft", legend = leg, lty = lty, lwd = lwd, col = colv, bty = "n")
  })
  
  
  
  output$forecast_info <- renderPrint({
    req(pred3())
    cat("Modelo estimado para predecir:\n")
    print(pred3()$fit)
    if ("linex" %in% pred3()$losses) {
      cat("\nRegla LINEX usada:  m* = mean - 0.5 * a * se^2 (para Y|info ~ Normal)\n")
    }
  })
}

shinyApp(ui, server)

