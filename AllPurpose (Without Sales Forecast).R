#All Purpose Modelling
library(shiny)
library(readxl)
library(dplyr)
library(earth)
library(markovchain)
library(ggplot2)
library(ggrepel)
library(scales)
library(tidyr)
library(lubridate)
library(shinyjs)

suppressWarnings({
  if (requireNamespace("forecast", quietly = TRUE)) library(forecast)
  if (requireNamespace("extraDistr", quietly = TRUE)) library(extraDistr)
})

# Helper: get mode for imputation
get_mode <- function(x) {
  ux <- unique(x[!is.na(x)])
  ux[which.max(tabulate(match(x, ux)))]
}

# Convert Earth model to Excel-compatible formula
convert_to_excel_formula <- function(model) {
  terms <- model$bx
  coefs <- model$coefficients
  
  terms_str <- character(length = length(coefs))
  terms_str[1] <- as.character(round(coefs[1], 5))  # Intercept
  
  if (length(coefs) > 1) {
    for (i in 2:length(coefs)) {
      term_expr <- as.character(terms[i])
      coef_val <- round(coefs[i], 5)
      
      # Replace h(x - c) with MAX(0, x - c)
      term_expr <- gsub("h\\(([^-\\)]+)-([^)]+)\\)", "MAX(0, \\1 - \\2)", term_expr)
      term_expr <- gsub("h\\(([^-\\)]+)\\+([^)]+)\\)", "MAX(0, \\1 + \\2)", term_expr)  # just in case
      term_expr <- gsub("h\\(([^)]+)\\)", "MAX(0, \\1)", term_expr)
      
      terms_str[i] <- paste0("+(", coef_val, " * ", term_expr, ")")
    }
  }
  formula_excel <- paste(terms_str, collapse = "")
  return(formula_excel)
}

# MARS fitting function
fit_mars_model <- function(data, predictors, target, degree = 1, trace = 1, na_action = c("omit", "fail", "impute")) {
  na_action <- match.arg(na_action)
  
  if (!all(c(predictors, target) %in% names(data))) stop("Some predictor or target variables not found.")
  
  df <- data[, c(target, predictors), drop=FALSE]
  original_n <- nrow(df)
  
  # Handle missing values
  if (na_action == "omit") {
    df <- na.omit(df)
  } else if (na_action == "fail" && anyNA(df)) {
    stop("Missing values found and na_action='fail'.")
  } else if (na_action == "impute") {
    for (col in names(df)) {
      if (anyNA(df[[col]])) {
        if (is.numeric(df[[col]])) {
          df[[col]][is.na(df[[col]])] <- mean(df[[col]], na.rm = TRUE)
        } else {
          if (is.logical(df[[col]])) df[[col]] <- factor(df[[col]])
          df[[col]][is.na(df[[col]])] <- get_mode(df[[col]])
        }
      }
    }
  }
  
  formula_str <- paste(target, "~", paste(predictors, collapse = " + "))
  formula_obj <- as.formula(formula_str)
  
  model <- earth(formula = formula_obj, data = df, degree = degree, trace = trace)
  
  list(model = model, formula = formula_obj, summary = summary(model), equation = capture.output(print(model)), n_obs = nrow(df))
}

# Compute transition matrix from multiple variables joint state
compute_transition_matrix_multi <- function(data, vars) {
  # Replace NAs with "NA" strings
  joint_states <- apply(data[, vars, drop=FALSE], 1, function(row) {
    row[is.na(row)] <- "NA"
    paste(as.character(row), collapse = "__")
  })
  
  joint_states <- factor(joint_states)
  states <- levels(joint_states)
  n <- length(states)
  
  trans_mat <- matrix(0, nrow = n, ncol = n, dimnames = list(states, states))
  
  if (length(joint_states) < 2) stop("Need at least two rows for Markov chain transitions.")
  
  for (i in 1:(length(joint_states)-1)) {
    from <- as.character(joint_states[i])
    to <- as.character(joint_states[i+1])
    trans_mat[from, to] <- trans_mat[from, to] + 1
  }
  
  # Normalize rows, if zero row sum make absorbing
  row_sums <- rowSums(trans_mat)
  for (r in seq_len(n)) {
    if (row_sums[r] == 0) {
      trans_mat[r, r] <- 1
    } else {
      trans_mat[r, ] <- trans_mat[r, ] / row_sums[r]
    }
  }
  
  trans_mat
}

# Build a filter label from current selections
build_filter_label <- function(filter_count, input_list) {
  parts <- c()
  
  # Add regular filters
  if (!is.null(filter_count) && filter_count > 0) {
    for (i in seq_len(filter_count)) {
      ns <- paste0("filter_", i)
      col <- input_list[[paste0(ns, "_col")]]
      vals <- input_list[[paste0(ns, "_vals")]]
      if (!is.null(col) && !is.null(vals) && length(vals) > 0) {
        if (length(vals) <= 3) {
          parts <- c(parts, paste0(col, ": ", paste(vals, collapse = ", ")))
        } else {
          parts <- c(parts, paste0(col, ": ", length(vals), " values"))
        }
      }
    }
  }
  
  # Add data cleaning note
  parts <- c(parts, "Excl. NA/Zero")
  
  if (!length(parts)) "All data" else paste(parts, collapse = " | ")
}

# Summarize occurrence vs volume (sum/mean)
summarise_occurrence_volume <- function(df, group_col, volume_col, agg = c("sum", "mean"), top_n = NA_integer_) {
  agg <- match.arg(agg)
  stopifnot(group_col %in% names(df), volume_col %in% names(df))
  g <- rlang::sym(group_col)
  v <- rlang::sym(volume_col)
  
  out <- df |>
    dplyr::filter(!is.na(!!g)) |>
    dplyr::summarise(
      .by = !!g,
      occurrence = dplyr::n(),
      volume_sum = sum(!!v, na.rm = TRUE),
      volume_mean = mean(!!v, na.rm = TRUE)
    ) |>
    dplyr::mutate(volume = ifelse(agg == "sum", volume_sum, volume_mean))
  
  if (!is.na(top_n) && top_n > 0) {
    out <- out |>
      dplyr::arrange(dplyr::desc(occurrence)) |>
      dplyr::slice_head(n = top_n)
  }
  out
}

# Monte Carlo: sample from data or parametric distributions
simulate_mc <- function(n, dist, params, sample_vec = NULL) {
  set.seed(params$seed %||% 123L)
  switch(
    dist,
    "bootstrap" = {
      stopifnot(!is.null(sample_vec))
      sample(sample_vec, size = n, replace = TRUE)
    },
    "normal" = rnorm(n, mean = params$mean %||% 0, sd = params$sd %||% 1),
    "lognormal" = rlnorm(n, meanlog = params$meanlog %||% 0, sdlog = params$sdlog %||% 1),
    "triangular" = {
      if (!requireNamespace("extraDistr", quietly = TRUE)) {
        # Fallback triangular using uniform transformation
        u <- runif(n)
        a <- params$min %||% 0
        b <- params$max %||% 1
        c <- params$mode %||% 0.5
        
        fc <- (c - a) / (b - a)
        ifelse(u < fc,
               a + sqrt(u * (b - a) * (c - a)),
               b - sqrt((1 - u) * (b - a) * (b - c)))
      } else {
        extraDistr::rltri(n, a = params$min, b = params$max, c = params$mode)
      }
    },
    stop("Unsupported distribution")
  )
}

`%||%` <- function(x, y) if (is.null(x)) y else x

# Simple sales forecast using forecast::auto.arima if available; fallback to naive
simple_sales_forecast <- function(df, date_col, value_col, horizon) {
  stopifnot(date_col %in% names(df), value_col %in% names(df))
  d <- df |>
    dplyr::mutate(
      !!date_col := as.Date(.data[[date_col]])
    ) |>
    dplyr::group_by(.data[[date_col]]) |>
    dplyr::summarise(value = sum(.data[[value_col]], na.rm = TRUE), .groups = "drop") |>
    dplyr::arrange(.data[[date_col]])
  
  if (nrow(d) < 3) stop("Not enough observations for forecasting.")
  
  # infer frequency (daily -> 7 or 30.44? We'll use weekly if strong weekly, otherwise 1)
  # Keep it simple: use ts with frequency = 7 if nearly daily; else frequency = 1
  date_diff <- median(diff(d[[date_col]]))
  freq <- if (!is.na(date_diff) && date_diff <= 2) 7 else 1
  
  ts_obj <- stats::ts(d$value, frequency = freq)
  
  if ("forecast" %in% loadedNamespaces()) {
    fit <- tryCatch(forecast::auto.arima(ts_obj), error = function(e) NULL)
  } else {
    fit <- NULL
  }
  
  fc <- if (!is.null(fit)) {
    forecast::forecast(fit, h = horizon)
  } else {
    # fallback naive
    last_val <- tail(d$value, 1)
    mean_sd <- stats::sd(d$value, na.rm = TRUE)
    mean_sd <- ifelse(is.na(mean_sd) || mean_sd == 0, 1, mean_sd)
    list(
      mean = rep(last_val, horizon),
      lower = rep(last_val - 1.96 * mean_sd, horizon),
      upper = rep(last_val + 1.96 * mean_sd, horizon)
    )
  }
  
  future_dates <- seq.Date(from = max(d[[date_col]]) + 1, by = "day", length.out = horizon)
  
  if (!is.null(fit)) {
    tibble::tibble(
      date = future_dates,
      .mean = as.numeric(fc$mean),
      .lower = as.numeric(fc$lower[, 2]),
      .upper = as.numeric(fc$upper[, 2])
    ) |>
      dplyr::mutate(method = "auto.arima")
  } else {
    tibble::tibble(
      date = future_dates,
      .mean = fc$mean,
      .lower = fc$lower,
      .upper = fc$upper
    ) |>
      dplyr::mutate(method = "naive")
  }
}
# UI
ui <- fluidPage(
  useShinyjs(),
  titlePanel("All-Purpose Analytics Dashboard"),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("file", "Upload CSV or Excel File", accept = c(".csv", ".xlsx")),
      
      hr(),
      h4("Data Filters"),
      numericInput("filter_count", "Number of Filters to Apply:", value = 0, min = 0, max = 10),
      uiOutput("filters_ui"),
      
      hr(),
      h4("Analysis Selection"),
      radioButtons("operation", "Select Operation:",
                   choices = c(
                     "MARS Model" = "mars",
                     "Markov Chain" = "markov",
                     "GGPlot: Volume vs Occurrence" = "ggplot",
                     "Monte Carlo" = "mc",
                     "Sales Predictive" = "sales"
                   ),
                   selected = "mars"),
      
      conditionalPanel(
        condition = "input.operation == 'mars'",
        uiOutput("target_ui"),
        uiOutput("predictors_ui"),
        selectInput("na_action", "Handle Missing Values:", choices = c("impute", "omit", "fail"), selected = "impute"),
        numericInput("degree", "Model Degree (1 = linear, 2 = interactions)", value = 1, min = 1, max = 2),
        actionButton("run_model", "Run MARS Model")
      ),
      
      conditionalPanel(
        condition = "input.operation == 'markov'",
        uiOutput("state_vars_ui"),
        actionButton("run_markov", "Compute Markov Chain Steady-State")
      ),
      
      conditionalPanel(
        condition = "input.operation == 'ggplot'",
        uiOutput("gg_group_ui"),
        uiOutput("gg_volume_ui"),
        selectInput("gg_plot_type", "Plot Type", 
                    choices = c("Scatter (Volume vs Occurrence)" = "scatter", 
                                "Histogram" = "histogram"), 
                    selected = "scatter"),
        conditionalPanel(
          condition = "input.operation == 'ggplot' && input.gg_plot_type == 'scatter'",
          selectInput("gg_agg", "Volume Aggregation", 
                      choices = c("Sum (Total)" = "sum", "Average (Mean)" = "mean"), 
                      selected = "sum"),
          numericInput("gg_topn", "Top N Groups (optional)", value = NA, min = 1)
        ),
        conditionalPanel(
          condition = "input.operation == 'ggplot' && input.gg_plot_type == 'histogram'",
          numericInput("gg_topn", "Top N Plants (optional)", value = NA, min = 1),
          checkboxInput("hist_density", "Show Volume Values on Bars", value = TRUE),
          checkboxInput("hist_facet", "Sort by Volume (Descending)", value = TRUE)
        ),
        actionButton("run_ggplot", "Build Plot")
      ),
      
      conditionalPanel(
        condition = "input.operation == 'mc'",
        uiOutput("mc_source_ui"),
        selectInput("mc_dist", "Distribution",
                    c("bootstrap", "normal", "lognormal", "triangular"),
                    selected = "bootstrap"),
        numericInput("mc_n", "Iterations", value = 10000, min = 100, step = 1000),
        numericInput("mc_seed", "Seed", value = 123, min = 1),
        # Parameter inputs
        conditionalPanel(
          condition = "input.operation == 'mc' && input.mc_dist == 'normal'",
          numericInput("mc_mean", "Mean", value = 0),
          numericInput("mc_sd", "Standard Deviation", value = 1, min = 0.01)
        ),
        conditionalPanel(
          condition = "input.operation == 'mc' && input.mc_dist == 'lognormal'",
          numericInput("mc_meanlog", "Mean of log", value = 0),
          numericInput("mc_sdlog", "SD of log", value = 1, min = 0.01)
        ),
        conditionalPanel(
          condition = "input.operation == 'mc' && input.mc_dist == 'triangular'",
          numericInput("mc_min", "Minimum", value = 0),
          numericInput("mc_mode", "Mode", value = 0.5),
          numericInput("mc_max", "Maximum", value = 1)
        ),
        actionButton("run_mc", "Run Simulation")
      ),
      
      conditionalPanel(
        condition = "input.operation == 'sales'",
        uiOutput("sales_date_ui"),
        uiOutput("sales_value_ui"),
        numericInput("sales_horizon", "Forecast Horizon (days)", value = 30, min = 7),
        actionButton("run_sales", "Forecast")
      )
    ),
    
    mainPanel(
      verbatimTextOutput("filtered_info"),
      
      conditionalPanel(
        condition = "input.operation == 'mars'",
        verbatimTextOutput("model_summary"),
        verbatimTextOutput("model_equation"),
        h4("MARS Diagnostic Plots"),
        plotOutput("diagnostic_plots"),
        br(), br(), br(), br(),br(), br(), br(),br(), br(), br(), br(),br(), br(), br(),
        br(), br(), br(),br(), br(), br(),
        h4("Excel-Compatible Equation"),
        textAreaInput("excel_formula", label = NULL, value = "", rows = 3, width = "100%"),
        actionButton("copy_excel", "Copy to Clipboard", icon = icon("copy"))
      ),
      
      conditionalPanel(
        condition = "input.operation == 'markov'",
        h4("Transition Matrix"),
        tableOutput("transition_matrix"),
        h4("Steady-State Distribution"),
        verbatimTextOutput("steady_state")
      ),
      
      conditionalPanel(
        condition = "input.operation == 'ggplot'",
        plotOutput("occ_vol_plot", height = "480px")
      ),
      
      conditionalPanel(
        condition = "input.operation == 'mc'",
        plotOutput("mc_plot", height = "420px"),
        h4("Summary Statistics"),
        tableOutput("mc_summary")
      ),
      
      conditionalPanel(
        condition = "input.operation == 'sales'",
        plotOutput("sales_plot", height = "480px"),
        h4("Forecast Table"),
        tableOutput("sales_table")
      )
    )
  )
)

# Server
server <- function(input, output, session) {
  
  options(shiny.maxRequestSize = 1000 * 1024^2)
  
  # Reactive data load
  data_raw <- reactive({
    req(input$file)
    ext <- tools::file_ext(input$file$name)
    tryCatch({
      if (ext == "csv") {
        read.csv(input$file$datapath, stringsAsFactors = FALSE)
      } else if (ext == "xlsx") {
        readxl::read_excel(input$file$datapath)
      } else {
        stop("Unsupported file type.")
      }
    }, error = function(e) {
      showNotification(paste("Error loading file:", e$message), type = "error")
      NULL
    })
  })
  
  # Render filters UI once based on filter_count
  output$filters_ui <- renderUI({
    n <- input$filter_count
    df <- data_raw()
    if (is.null(df) || n == 0) return(NULL)
    
    filter_ui_list <- lapply(seq_len(n), function(i) {
      ns <- paste0("filter_", i)
      tagList(
        selectInput(paste0(ns, "_col"), paste0("Filter ", i, " - Select Column"), choices = names(df), selected = names(df)[1]),
        uiOutput(paste0(ns, "_vals_ui")),
        tags$hr()
      )
    })
    do.call(tagList, filter_ui_list)
  })
  
  # For each filter, dynamically update the values selectizeInput based on selected column
  observe({
    n <- input$filter_count
    df <- data_raw()
    if (is.null(df) || n == 0) return(NULL)
    
    for (i in seq_len(n)) {
      local({
        ii <- i
        ns <- paste0("filter_", ii)
        col_input <- input[[paste0(ns, "_col")]]
        if (is.null(col_input) || is.null(df)) return()
        
        vals <- unique(df[[col_input]])
        # Make sure values are character for consistent display
        vals <- as.character(vals)
        
        # Preserve previous selections if possible
        previous_sel <- isolate(input[[paste0(ns, "_vals")]])
        if (is.null(previous_sel)) {
          previous_sel <- vals
        } else {
          previous_sel <- intersect(previous_sel, vals)
          if (length(previous_sel) == 0) previous_sel <- vals
        }
        
        output[[paste0(ns, "_vals_ui")]] <- renderUI({
          selectizeInput(paste0(ns, "_vals"), "Select Values", choices = vals,
                         selected = previous_sel, multiple = TRUE,
                         options = list(plugins = list('remove_button')))
        })
      })
    }
  })
  
  # Filter data reactively based on all filters
  filtered_data <- reactive({
    df <- data_raw()
    req(df)
    n <- input$filter_count
    if (is.null(n) || n == 0) return(df)
    
    for (i in seq_len(n)) {
      ns <- paste0("filter_", i)
      col <- input[[paste0(ns, "_col")]]
      vals <- input[[paste0(ns, "_vals")]]
      if (!is.null(col) && !is.null(vals)) {
        df <- df %>% filter(.data[[col]] %in% vals)
      }
    }
    df
  })
  
  # Show filtered data info
  output$filtered_info <- renderPrint({
    df_raw <- data_raw()
    df_filtered <- filtered_data()
    if (!is.null(df_raw) && !is.null(df_filtered)) {
      cat("Data Info:\n")
      cat("Original rows:", nrow(df_raw), "\n")
      cat("Filtered rows:", nrow(df_filtered), "\n")
      if (input$filter_count > 0) {
        cat("Active filters:", build_filter_label(input$filter_count, input), "\n")
      }
    }
  })
  # Update target and predictors inputs choices based on filtered data
  output$target_ui <- renderUI({
    df <- filtered_data()
    if (is.null(df)) return(NULL)
    var_types <- sapply(df, function(col) class(col)[1])
    choices <- setNames(names(var_types), paste0(names(var_types), " (", var_types, ")"))
    selectizeInput("target", "Select Target Variable", choices = choices, options = list(placeholder = 'Type to search...'))
  })
  
  output$predictors_ui <- renderUI({
    df <- filtered_data()
    if (is.null(df)) return(NULL)
    var_types <- sapply(df, function(col) class(col)[1])
    choices <- setNames(names(var_types), paste0(names(var_types), " (", var_types, ")"))
    selectizeInput("predictors", "Select Predictor Variable(s)", choices = choices, multiple = TRUE, options = list(placeholder = 'Type to search...'))
  })
  
  # Show multiple categorical/logical vars for Markov chain states
  output$state_vars_ui <- renderUI({
    df <- filtered_data()
    if (is.null(df)) return(NULL)
    
    # Keep only logical, character, or factor columns
    valid_cols <- names(df)[sapply(df, function(col) {
      is.logical(col) || is.factor(col) || is.character(col)
    })]
    
    if (length(valid_cols) == 0) {
      return(tags$p("No categorical or logical variables found for Markov Chain analysis."))
    }
    
    var_types <- sapply(df[valid_cols], function(col) class(col)[1])
    choices <- setNames(valid_cols, paste0(valid_cols, " (", var_types, ")"))
    
    selectInput("state_vars", "Select State Variable(s) (multiple allowed):", choices = choices, multiple = TRUE)
  })
  
  # GGPlot selectors - Updated to include Plant
  output$gg_group_ui <- renderUI({
    df <- filtered_data()
    if (is.null(df)) return(NULL)
    
    # Get all potential grouping columns
    categorical_cols <- names(df)[sapply(df, function(x) {
      is.character(x) || is.factor(x) || is.logical(x)
    })]
    
    # Force include specific columns even if they're numeric or other types
    force_include <- c("Plant", "Fin.Month", "Fin.Month.Number", "Quarter", "Quarters", "Month", "Year")
    numeric_as_categorical <- force_include[force_include %in% names(df)]
    
    # Combine all valid options
    valid <- unique(c(categorical_cols, numeric_as_categorical))
    
    # Prioritize specific columns (in order of preference)
    priority_cols <- c("Plant", "Fin.Month", "Fin.Month.Number", "Quarter", "Quarters")
    available_priority <- priority_cols[priority_cols %in% valid]
    
    if (length(available_priority) > 0) {
      # Put priority columns first, then the rest
      valid <- c(available_priority, valid[!valid %in% available_priority])
    }
    
    if (input$gg_plot_type == "histogram") {
      # Required for histograms (bar chart style)
      if (!length(valid)) return(tags$p("No grouping columns available"))
      
      # Default selection priority: Plant > Fin.Month > Fin.Month.Number > others
      default_selection <- if("Plant" %in% valid) {
        "Plant"
      } else if("Fin.Month" %in% valid) {
        "Fin.Month"
      } else if("Fin.Month.Number" %in% valid) {
        "Fin.Month.Number"
      } else {
        valid[1]
      }
      
      selectInput("gg_group", "Group (X-axis)", 
                  choices = valid, selected = default_selection)
    } else {
      # Required for scatter
      if (!length(valid)) return(tags$p("No grouping columns available"))
      
      # Same priority for scatter plots
      default_selection <- if("Plant" %in% valid) {
        "Plant"
      } else if("Fin.Month" %in% valid) {
        "Fin.Month"
      } else if("Fin.Month.Number" %in% valid) {
        "Fin.Month.Number"
      } else {
        valid[1]
      }
      
      selectInput("gg_group", "Group (categorical)", 
                  choices = valid, selected = default_selection)
    }
  })
  
  output$gg_volume_ui <- renderUI({
    df <- filtered_data()
    if (is.null(df)) return(NULL)
    valid <- names(df)[sapply(df, is.numeric)]
    if (!length(valid)) return(tags$p("No numeric columns available for volume"))
    selectInput("gg_volume", "Volume (numeric)", choices = valid)
  })
  
  # Monte Carlo: choose source column for bootstrap
  output$mc_source_ui <- renderUI({
    df <- filtered_data()
    if (is.null(df)) return(NULL)
    valid <- names(df)[sapply(df, is.numeric)]
    selectInput("mc_source", "Sample Column (for bootstrap)", choices = c("(none)" = "", valid), selected = "")
  })
  
  # Sales predictive: date and value columns
  output$sales_date_ui <- renderUI({
    df <- filtered_data()
    if (is.null(df)) return(NULL)
    # date-like columns
    valid <- names(df)[sapply(df, function(x) inherits(x, "Date") || inherits(x, "POSIXt") || is.character(x))]
    selectInput("sales_date", "Date Column", choices = valid)
  })
  
  output$sales_value_ui <- renderUI({
    df <- filtered_data()
    if (is.null(df)) return(NULL)
    valid <- names(df)[sapply(df, is.numeric)]
    selectInput("sales_value", "Sales/Value Column", choices = valid)
  })
  
  # Run MARS model
  observeEvent(input$run_model, {
    req(filtered_data(), input$target, input$predictors)
    
    tryCatch({
      res <- fit_mars_model(
        data = filtered_data(),
        predictors = input$predictors,
        target = input$target,
        degree = input$degree,
        na_action = input$na_action
      )
      
      output$model_summary <- renderPrint({ res$summary })
      output$model_equation <- renderPrint({ res$equation })
      
      # Update Excel formula text area
      excel_formula <- convert_to_excel_formula(res$model)
      updateTextAreaInput(session, "excel_formula", value = excel_formula)
      
      # Diagnostic plots
      output$diagnostic_plots <- renderPlot({
        req(res$model)
        
        rsq <- round(res$model$rsq, 4)
        grsq <- round(res$model$grsq, 4)
        
        layout(matrix(1:4, nrow = 2, byrow = TRUE))
        par(mar = c(5, 5, 4, 2), oma = c(2, 2, 4, 2), cex.main = 1.2, cex.axis = 1, cex.lab = 1)
        
        plot(res$model, which = 1, main = "Plot 1: R² vs Terms")
        mtext("Improvement in R² and GRSq as more terms are added.", side = 1, line = 4, cex = 1)
        
        plot(res$model, which = 2, main = "Plot 2: CDF of Residuals")
        mtext("Smooth S-curve indicates good residual distribution.", side = 1, line = 4, cex = 1)
        
        plot(res$model, which = 3, main = "Plot 3: Residuals vs Fitted")
        mtext("Look for spread and patterns; should be random scatter.", side = 1, line = 4, cex = 1)
        
        plot(res$model, which = 4, main = "Plot 4: Q-Q Plot of Residuals")
        mtext(paste("Normality check of residuals (R² =", rsq, ", GRSq =", grsq, ")"), side = 1, line = 4, cex = 1)
        
        mtext("MARS Diagnostic Plots", outer = TRUE, line = 1, cex = 1.5, font = 2)
      }, height = 800)
      
    }, error = function(e) {
      showNotification(paste("MARS Model Error:", e$message), type = "error")
      output$model_summary <- renderPrint({ NULL })
      output$model_equation <- renderPrint({ NULL })
      output$diagnostic_plots <- renderPlot(NULL)
      updateTextAreaInput(session, "excel_formula", value = "")
    })
  })
  
  # Run Markov Chain
  observeEvent(input$run_markov, {
    req(filtered_data(), input$state_vars)
    
    df <- filtered_data()
    vars <- input$state_vars
    
    # Convert logical to character for consistency
    for (v in vars) {
      if (is.logical(df[[v]])) df[[v]] <- as.character(df[[v]])
      else df[[v]] <- as.character(df[[v]])
    }
    
    tryCatch({
      P <- compute_transition_matrix_multi(df, vars)
      mc <- new("markovchain", transitionMatrix = P, states = rownames(P), name = "FromData")
      steady <- steadyStates(mc)
      
      output$transition_matrix <- renderTable(round(P, 4), rownames = TRUE)
      output$steady_state <- renderPrint({
        cat("Steady-State Distribution:\n")
        print(round(steady, 5))
      })
      
    }, error = function(e) {
      showNotification(paste("Markov Chain Error:", e$message), type = "error")
      output$transition_matrix <- renderTable(NULL)
      output$steady_state <- renderPrint(NULL)
    })
  })
  
  # Run GGPlot
  observeEvent(input$run_ggplot, {
    req(filtered_data(), input$gg_volume, input$gg_plot_type)
    df <- filtered_data()
    
    # Remove NA and zero values from volume column
    df <- df %>% 
      filter(!is.na(.data[[input$gg_volume]]) & .data[[input$gg_volume]] > 0)
    
    # Validate we still have data
    if (nrow(df) == 0) {
      showNotification("No data remaining after filtering NA/zero values", type = "warning")
      return()
    }
    
    if (input$gg_plot_type == "histogram") {
      # Bar chart logic - Plant on X-axis, Volume sum on Y-axis
      req(input$gg_group)
      
      # Remove NA and zero values from volume column
      df <- df %>% 
        filter(!is.na(.data[[input$gg_volume]]) & .data[[input$gg_volume]] > 0)
      
      # Remove NA values from grouping variable
      if (!is.null(input$gg_group) && nzchar(input$gg_group)) {
        if (is.logical(df[[input$gg_group]])) df[[input$gg_group]] <- as.character(df[[input$gg_group]])
        df <- df %>% filter(!is.na(.data[[input$gg_group]]))
      }
      
      label_txt <- build_filter_label(input$filter_count, input)
      
      output$occ_vol_plot <- renderPlot({
        
        # Aggregate volume by group (Plant)
        plot_data <- df %>%
          group_by(.data[[input$gg_group]]) %>%
          summarise(
            total_volume = sum(.data[[input$gg_volume]], na.rm = TRUE),
            count = n(),
            .groups = "drop"
          ) %>%
          arrange(desc(total_volume))
        
        # Limit to top N if specified
        if (!is.na(input$gg_topn) && input$gg_topn > 0) {
          plot_data <- plot_data %>% slice_head(n = input$gg_topn)
        }
        
        # Create bar chart with Plant highlighting
        p <- ggplot(plot_data, aes(x = reorder(.data[[input$gg_group]], total_volume), y = total_volume)) +
          geom_col(aes(fill = .data[[input$gg_group]]), color = "white", alpha = 0.8, width = 0.7) +
          scale_fill_manual(values = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", 
                                       "#ffff33", "#a65628", "#f781bf", "#999999", "#66c2a5")) +
          coord_flip() +  # Flip for better plant name readability
          guides(fill = "none")  # Remove legend since x-axis shows plant names
        
        # Add value labels on bars if requested
        if (input$hist_density) {  # Reusing this checkbox for value labels
          p <- p + geom_text(aes(label = scales::comma(round(total_volume, 0))), 
                             hjust = -0.1, size = 3.5, color = "black", fontface = "bold")
        }
        
        # Smart zoom limits for Y-axis
        y_max <- max(plot_data$total_volume, na.rm = TRUE) * 1.15
        
        p + 
          scale_y_continuous(labels = scales::comma, limits = c(0, y_max), expand = c(0, 0)) +
          labs(
            title = paste("Total", input$gg_volume, "by", input$gg_group),
            x = input$gg_group,
            y = paste("Total", input$gg_volume),
            subtitle = label_txt
          ) +
          theme_minimal(base_size = 12) +
          theme(
            axis.text.y = element_text(size = 11, face = "bold"),  # Bold plant names
            axis.text.x = element_text(size = 10),
            plot.title = element_text(size = 14, face = "bold"),
            panel.grid.major.y = element_blank(),  # Remove horizontal grid lines
            panel.grid.minor = element_blank(),
            axis.ticks.y = element_line(color = "gray50")
          )
      })
      
    } else {
      # Scatter plot logic - Group on X-axis, Volume on Y-axis (same as histogram but as points)
      req(input$gg_group, input$gg_agg)
      
      # Remove NA and zero values from volume column
      df <- df %>% 
        filter(!is.na(.data[[input$gg_volume]]) & .data[[input$gg_volume]] > 0)
      
      # Remove NA values from grouping variable
      if (!is.null(input$gg_group) && nzchar(input$gg_group)) {
        if (is.logical(df[[input$gg_group]])) df[[input$gg_group]] <- as.character(df[[input$gg_group]])
        df <- df %>% filter(!is.na(.data[[input$gg_group]]))
      }
      
      label_txt <- build_filter_label(input$filter_count, input)
      
      output$occ_vol_plot <- renderPlot({
        
        # Aggregate volume by group - sum or average (same logic as histogram)
        if (input$gg_agg == "sum") {
          plot_data <- df %>%
            group_by(.data[[input$gg_group]]) %>%
            summarise(
              volume_agg = sum(.data[[input$gg_volume]], na.rm = TRUE),
              count = n(),
              .groups = "drop"
            )
          y_label <- paste("Total", input$gg_volume)
          title_prefix <- "Total"
        } else {
          plot_data <- df %>%
            group_by(.data[[input$gg_group]]) %>%
            summarise(
              volume_agg = mean(.data[[input$gg_volume]], na.rm = TRUE),
              count = n(),
              .groups = "drop"
            )
          y_label <- paste("Average", input$gg_volume)
          title_prefix <- "Average"
        }
        
        # Sort by volume (descending)
        plot_data <- plot_data %>% arrange(desc(volume_agg))
        
        # Limit to top N if specified
        if (!is.na(input$gg_topn) && input$gg_topn > 0) {
          plot_data <- plot_data %>% slice_head(n = input$gg_topn)
        }
        
        # Create scatter plot with same data as histogram
        ggplot(plot_data, aes(x = reorder(.data[[input$gg_group]], volume_agg), y = volume_agg)) +
          geom_point(aes(color = .data[[input$gg_group]]), size = 4, alpha = 0.8) +
          geom_text_repel(aes(label = .data[[input$gg_group]]), size = 3, fontface = "bold") +
          scale_color_manual(values = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", 
                                        "#ffff33", "#a65628", "#f781bf", "#999999", "#66c2a5",
                                        "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854", "#ffd92f")) +
          guides(color = "none") +  # Remove legend since labels show group names
          scale_x_discrete() +
          scale_y_continuous(labels = scales::comma) +
          labs(
            title = paste(title_prefix, input$gg_volume, "by", input$gg_group, "(Scatter)"),
            x = input$gg_group,
            y = y_label,
            subtitle = paste(label_txt, "|", toupper(input$gg_agg))
          ) +
          theme_minimal(base_size = 12) +
          theme(
            axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
            axis.text.y = element_text(size = 10),
            plot.title = element_text(size = 14, face = "bold"),
            panel.grid.minor = element_blank()
          )
      })
    }
  })
  
  # Run Monte Carlo
  observeEvent(input$run_mc, {
    req(filtered_data(), input$mc_n, input$mc_dist, input$mc_seed)
    
    df <- filtered_data()
    sample_vec <- NULL
    if (!is.null(input$mc_source) && nzchar(input$mc_source)) {
      validate(need(is.numeric(df[[input$mc_source]]), "Selected sample column must be numeric"))
      sample_vec <- df[[input$mc_source]]
    }
    
    params <- list(
      seed = input$mc_seed,
      mean = if (!is.null(input$mc_mean)) input$mc_mean else NULL,
      sd = if (!is.null(input$mc_sd)) input$mc_sd else NULL,
      meanlog = if (!is.null(input$mc_meanlog)) input$mc_meanlog else NULL,
      sdlog = if (!is.null(input$mc_sdlog)) input$mc_sdlog else NULL,
      min = if (!is.null(input$mc_min)) input$mc_min else NULL,
      mode = if (!is.null(input$mc_mode)) input$mc_mode else NULL,
      max = if (!is.null(input$mc_max)) input$mc_max else NULL
    )
    
    sim <- tryCatch({
      simulate_mc(input$mc_n, input$mc_dist, params, sample_vec)
    }, error = function(e) {
      showNotification(paste("Monte Carlo Error:", e$message), type = "error")
      return(NULL)
    })
    req(sim)
    
    qs <- stats::quantile(sim, probs = c(0.05, 0.25, 0.5, 0.75, 0.95), na.rm = TRUE)
    tbl <- data.frame(
      Metric = c("Mean", "SD", "P5", "P25", "P50", "P75", "P95"),
      Value = c(mean(sim), sd(sim), qs[1], qs[2], qs[3], qs[4], qs[5])
    )
    output$mc_summary <- renderTable({
      data.frame(Metric = tbl$Metric, Value = round(tbl$Value, 4))
    }, rownames = FALSE)
    
    output$mc_plot <- renderPlot({
      label_txt <- build_filter_label(input$filter_count, input)
      
      ggplot(data.frame(x = sim), aes(x = x)) +
        geom_histogram(bins = 50, fill = "#6baed6", color = "white") +
        geom_vline(xintercept = as.numeric(qs[c(1,3,5)]), linetype = c("dashed","solid","dashed"),
                   color = c("#fb6a4a","#ef3b2c","#fb6a4a")) +
        labs(
          title = paste("Monte Carlo -", input$mc_dist),
          subtitle = paste0("n=", input$mc_n, " | ", label_txt),
          x = "Value", y = "Count"
        ) +
        theme_minimal(base_size = 12)
    })
  })
  # Run Sales Forecast
  observeEvent(input$run_sales, {
    req(filtered_data(), input$sales_date, input$sales_value, input$sales_horizon)
    df <- filtered_data()
    
    # Coerce date column if needed
    if (is.character(df[[input$sales_date]])) {
      suppressWarnings({
        parsed <- as.Date(df[[input$sales_date]])
        if (all(is.na(parsed))) {
          parsed <- as.Date(lubridate::ymd(df[[input$sales_date]]))
        }
      })
      df[[input$sales_date]] <- parsed
    }
    validate(need(!all(is.na(df[[input$sales_date]])), "Could not parse the selected date column"))
    
    fc_df <- tryCatch({
      simple_sales_forecast(df, input$sales_date, input$sales_value, horizon = input$sales_horizon)
    }, error = function(e) {
      showNotification(paste("Forecast Error:", e$message), type = "error")
      return(NULL)
    })
    req(fc_df)
    
    # Build historical series (daily sum)
    hist_df <- df |>
      dplyr::mutate(date = as.Date(.data[[input$sales_date]])) |>
      dplyr::group_by(date) |>
      dplyr::summarise(value = sum(.data[[input$sales_value]], na.rm = TRUE), .groups = "drop") |>
      dplyr::arrange(date)
    
    output$sales_table <- renderTable({
      dplyr::transmute(fc_df, 
                       Date = date, 
                       Forecast = round(.mean, 3), 
                       Lower95 = round(.lower, 3), 
                       Upper95 = round(.upper, 3), 
                       Method = method)
    })
    
    output$sales_plot <- renderPlot({
      label_txt <- build_filter_label(input$filter_count, input)
      
      ggplot() +
        geom_line(data = hist_df, aes(x = date, y = value), color = "#2c7fb8", size = 1) +
        geom_ribbon(data = fc_df, aes(x = date, ymin = .lower, ymax = .upper), fill = "#9ecae1", alpha = 0.5) +
        geom_line(data = fc_df, aes(x = date, y = .mean), color = "#ef3b2c", size = 1) +
        labs(
          title = "Sales Forecast",
          subtitle = paste0("Horizon: ", input$sales_horizon, " days | ", label_txt),
          x = "Date", y = "Sales (aggregated daily)"
        ) +
        scale_y_continuous(labels = scales::comma) +
        theme_minimal(base_size = 12)
    })
  })
  
  # Clipboard copy handler
  observeEvent(input$copy_excel, {
    clip <- input$excel_formula
    session$sendCustomMessage("copyToClipboard", clip)
  })
}

# Simple JavaScript for clipboard copying (just a few lines)
jsCode <- "
Shiny.addCustomMessageHandler('copyToClipboard', function(message) {
  const el = document.createElement('textarea');
  el.value = message;
  document.body.appendChild(el);
  el.select();
  document.execCommand('copy');
  document.body.removeChild(el);
  alert('Excel formula copied to clipboard!');
});
"

# Final app
shinyApp(
  ui = tagList(
    tags$head(tags$script(HTML(jsCode))),
    ui
  ),
  server = server
)