# nolint start
# Funkcje pomocnicze – dostosowane pod Shiny:
# - compute_edge_weights: liczy w_ij = cnt_ij / cnt_i i mapuje do kolejności krawędzi w grafie
# - precompute_structures: prekomputacja ID/outgoing/weights (z grafu)
# - ic_run_fast: szybka IC z możliwością ustawienia max_iter i gotowym precompute
# - compute_rankings: liczy rankingi + zwraca pickery seedów
# - average_runs: średnia z wielu powtórzeń z parametrem ic_max_iter

library(shiny)
library(igraph)

# ------------- liczenie wag i mapowanie ----------------
compute_edge_weights <- function(df, g) {
  edges_cnt <- as.data.frame(table(df$src, df$dst), stringsAsFactors = FALSE)
  colnames(edges_cnt) <- c("src","dst","cnt_ij")
  edges_cnt <- subset(edges_cnt, cnt_ij > 0 & src != dst)
  
  send_cnt <- as.data.frame(table(df$src), stringsAsFactors = FALSE)
  colnames(send_cnt) <- c("src","cnt_i")
  
  edges_w <- merge(edges_cnt, send_cnt, by = "src")
  edges_w$weight <- edges_w$cnt_ij / edges_w$cnt_i
  
  uv_names <- ends(g, es = E(g), names = TRUE)
  key_g <- paste(uv_names[,1], uv_names[,2])
  key_w <- paste(edges_w$src, edges_w$dst)
  
  wmatch <- edges_w$weight[ match(key_g, key_w) ]
  wmatch[is.na(wmatch)] <- 0
  pmin(pmax(wmatch, 0), 1)
}

# ------------- prekomputacja struktur -------------------
precompute_structures <- function(g) {
  ne <- ecount(g); nv <- vcount(g)
  et <- ends(g, es = E(g), names = FALSE)
  from_id <- et[,1]; to_id <- et[,2]
  
  tmp_split <- split(seq_len(ne), from_id)
  out_eids <- vector("list", nv)
  for (i in seq_len(nv)) {
    idx <- tmp_split[[as.character(i)]]
    if (is.null(idx)) idx <- integer(0)
    out_eids[[i]] <- idx
  }
  
  list(
    nv = nv,
    ne = ne,
    from_id = from_id,
    to_id   = to_id,
    out_eids = out_eids,
    w = { w <- E(g)$weight; w[is.na(w)] <- 0; w }
  )
}

# ------------- szybki bieg IC (z precompute) ------------
ic_run_fast <- function(g, seeds, max_iter, precomp) {
  nv <- precomp$nv
  out_eids <- precomp$out_eids
  to_id <- precomp$to_id
  w <- precomp$w
  
  V(g)$activated   <- FALSE
  V(g)$to_activate <- FALSE
  E(g)$tried       <- FALSE
  
  seeds_id <- as.integer(seeds)
  V(g)$activated[seeds_id] <- TRUE
  frontier <- seeds_id
  
  active_counts <- c(sum(V(g)$activated))
  it <- 0
  
  repeat {
    it <- it + 1
    if (length(frontier) == 0 || it > max_iter) break
    
    V(g)$to_activate[] <- FALSE
    
    for (u in frontier) {
      eids <- out_eids[[u]]
      if (length(eids) == 0) next
      
      eids <- eids[!E(g)$tried[eids]]
      if (length(eids) == 0) next
      E(g)$tried[eids] <- TRUE
      
      vs <- to_id[eids]
      ww <- w[eids]
      
      keep <- !V(g)$activated[vs]
      if (any(keep)) {
        vs2 <- vs[keep]
        ww2 <- ww[keep]
        succ <- runif(length(ww2)) < ww2
        if (any(succ)) V(g)$to_activate[vs2[succ]] <- TRUE
      }
    }
    
    newly <- which(V(g)$to_activate)
    if (length(newly) == 0) break
    V(g)$activated[newly] <- TRUE
    frontier <- newly
    
    active_counts <- c(active_counts, sum(V(g)$activated))
  }
  
  active_counts
}

# ------------- rankingi + pickery -----------------------
compute_rankings <- function(g) {
  sc_out <- degree(g, mode="out")
  sc_bet <- betweenness(g, directed = TRUE)
  sc_clo <- closeness(g, mode="out")
  sc_pr  <- page_rank(g, directed = TRUE, weights = E(g)$weight)$vector
  
  ord_out <- order(sc_out, decreasing = TRUE)
  ord_bet <- order(sc_bet, decreasing = TRUE)
  ord_clo <- order(sc_clo, decreasing = TRUE)
  ord_pr  <- order(sc_pr,  decreasing = TRUE)
  
  make_picker_from_order <- function(ord) {
    force(ord)
    function(g, pct=0.05) {
      k <- max(1, ceiling(pct * vcount(g)))
      V(g)[ord[1:k]]
    }
  }
  
  list(
    pick_outdeg      = make_picker_from_order(ord_out),
    pick_betweenness = make_picker_from_order(ord_bet),
    pick_closeness   = make_picker_from_order(ord_clo),
    pick_pagerank    = make_picker_from_order(ord_pr),
    pick_random      = function(g, pct = 0.05, seed = NULL) {
      if (!is.null(seed)) set.seed(seed)
      k <- max(1, ceiling(pct * vcount(g)))
      sample(V(g), k)
    }
  )
}

# ------------- średnia z N biegów -----------------------
average_runs <- function(
    g, pick_fun, runs = 100, pct = 0.05,
    ic_max_iter = vcount(g), deterministic = TRUE, precomp
) {
  traces <- vector("list", runs)
  maxlen <- 0
  seeds_fixed <- if (deterministic) pick_fun(g, pct = pct) else NULL
  
  for (r in seq_len(runs)) {
    if (!deterministic && (r == 1L)) {
      s1 <- pick_fun(g, pct = pct)
      tr <- ic_run_fast(g, s1, max_iter = ic_max_iter, precomp = precomp)
    } else {
      seeds <- if (deterministic) seeds_fixed else pick_fun(g, pct = pct)
      tr <- ic_run_fast(g, seeds, max_iter = ic_max_iter, precomp = precomp)
    }
    traces[[r]] <- tr
    if (length(tr) > maxlen) maxlen <- length(tr)
  }
  
  mat <- matrix(NA_real_, nrow = runs, ncol = maxlen)
  for (r in seq_len(runs)) {
    tr <- traces[[r]]
    if (length(tr) < maxlen) tr <- c(tr, rep(tail(tr, 1), maxlen - length(tr)))
    mat[r, ] <- tr
  }
  colMeans(mat)
}

# Shiny: IC na grafie z regulacją prawdopodobieństw w_ij i liczbą iteracji

# ====== UI ======
ui <- fluidPage(
  titlePanel("Dyfuzja informacji (Independent Cascade)"),
  
  sidebarLayout(
    sidebarPanel(
      helpText("Dane domyślnie są wczytywane z URL udostępnionego w zadaniu."),
      
      # 1) Źródło danych
      textInput(
        "data_url",
        "URL danych:",
        value = "https://bergplace.org/share/out.radoslaw_email_email"
      ),
      fileInput(
        "datafile",
        "lub wgraj plik out.radoslaw_email_email:",
        accept = c(".txt", ".csv", "")
      ),
      
      # 2) Suwak do skalowania wag w_ij
      sliderInput(
        "wij_factor",
        "Skalowanie prawdopodobieństw w_ij:",
        min = 0.10, max = 2.00, value = 1.00, step = 0.05
      ),
      helpText("0.10 = 10% w_ij, 1.00 = 100% (oryginalne), 2.00 = 200%."),
      
      # 3) Suwak liczby iteracji IC
      sliderInput(
        "iters",
        "Maksymalna liczba iteracji IC:",
        min = 1, max = 50, value = 10, step = 1
      ),
      
      # 4) Parametry eksperymentu
      numericInput("runs", "Liczba powtórzeń (średnia):", value = 100, min = 1, step = 1),
      sliderInput("pct_seeds", "Odsetek seedów (0.5%–50%):", min = 0.005, max = 0.50, value = 0.05, step = 0.005),
      numericInput("random_seed", "Ziarno losowości:", value = 123, step = 1),
      
      # 5) Wybór strategii seedów
      checkboxGroupInput(
        "strategies",
        "Strategie wyboru seedów:",
        choices = c("Top outdegree", "Top betweenness", "Top closeness", "Losowe", "PageRank"),
        selected = c("Top outdegree", "Top betweenness", "Top closeness", "Losowe", "PageRank")
      ),
      
      actionButton("run", "Uruchom symulacje", class = "btn-primary")
    ),
    
    mainPanel(
      plotOutput("plot", height = 420),
      hr(),
      tableOutput("summary"),
      hr(),
      verbatimTextOutput("log")
    )
  )
)

# ====== SERVER ======
server <- function(input, output, session) {
  
  # 1) Reaktywne wczytanie danych
  df_reactive <- reactive({
    # priorytet ma plik wgrany przez użytkownika
    if (!is.null(input$datafile)) {
      read.table(input$datafile$datapath, skip = 2)[, 1:2]
    } else {
      # z URL (domyślne)
      df <- read.table(input$data_url, skip = 2)[, 1:2]
    }
  })
  
  # 2) Uruchomienie całego pipeline po kliknięciu
  observeEvent(input$run, {
    req(df_reactive())
    
    # Log buffer
    .logs <- character(0)
    add_log <- function(fmt, ...) {
      ts <- format(Sys.time(), "%H:%M:%S")
      .logs <<- c(.logs, sprintf("[%s] %s", ts, sprintf(fmt, ...)))
    }
    
    # --- Import i graf
    add_log("Start skryptu")
    df <- df_reactive()
    colnames(df) <- c("src", "dst")
    add_log("Węzły/łuki z pliku: %d wierszy", nrow(df))
    
    g <- simplify(
      graph_from_data_frame(df, directed = TRUE),
      remove.multiple = TRUE, remove.loops = TRUE
    )
    add_log("Graf: v=%d, e=%d", vcount(g), ecount(g))
    
    # --- Wagi w_ij i przypisanie do krawędzi
    wmatch <- compute_edge_weights(df, g)
    E(g)$weight <- pmin(pmax(wmatch, 0), 1)
    
    # --- Skalowanie w_ij suwakiem
    E(g)$weight <- pmin(E(g)$weight * input$wij_factor, 1)
    add_log("Skalowanie w_ij × %.2f (obcięte do [0,1])", input$wij_factor)
    
    # --- Prekomputacja (ID, listy krawędzi wychodzących, wektor wag)
    pre <- precompute_structures(g)
    
    # --- Rankingi seedów
    set.seed(input$random_seed)
    ranks <- compute_rankings(g)
    
    # --- Bieg strategii
    chosen <- input$strategies
    series_list <- list()
    
    run_strategy <- function(label, pick_fun, deterministic) {
      avg <- average_runs(
        g, pick_fun,
        runs = input$runs,
        pct = input$pct_seeds,
        ic_max_iter = input$iters,
        deterministic = deterministic,
        precomp = pre
      )
      c(avg[1], diff(avg))
    }
    
    if ("Top outdegree" %in% chosen)
      series_list[["Top outdegree"]]   <- run_strategy("Top outdegree",      ranks$pick_outdeg,      TRUE)
    if ("Top betweenness" %in% chosen)
      series_list[["Top betweenness"]] <- run_strategy("Top betweenness",    ranks$pick_betweenness, TRUE)
    if ("Top closeness" %in% chosen)
      series_list[["Top closeness"]]   <- run_strategy("Top closeness",      ranks$pick_closeness,   TRUE)
    if ("Losowe" %in% chosen)
      series_list[["Losowe"]]          <- run_strategy("Losowe",             ranks$pick_random,      FALSE)
    if ("PageRank" %in% chosen)
      series_list[["PageRank"]]        <- run_strategy("PageRank",           ranks$pick_pagerank,    TRUE)
    
    # --- Wykres
    maxlen <- max(sapply(series_list, length))
    pad_to <- function(x, n) if (length(x) < n) c(x, rep(0, n - length(x))) else x
    series_padded <- lapply(series_list, pad_to, n = maxlen)
    mat <- do.call(cbind, series_padded)
    colnames(mat) <- names(series_list)
    
    output$plot <- renderPlot({
      matplot(
        x = seq_len(nrow(mat)),
        y = mat,
        type = "o", lty = 1, pch = 1,
        xlab = "Iteracja",
        ylab = "Liczba aktywowanych węzłów w iteracji",
        main = "Dyfuzja informacji (Independent Cascade)"
      )
      legend("topright", legend = colnames(mat), lty = 1, pch = 1, col = seq_len(ncol(mat)))
    })
    
    # --- Podsumowanie
    output$summary <- renderTable({
      data.frame(
        strategia      = names(series_list),
        suma_aktywnych = sapply(series_list, sum),
        max_iteracji   = sapply(series_list, function(x) max(which(x > 0)))
      )
    }, rownames = FALSE)
    
    # --- Logi
    output$log <- renderText(paste(.logs, collapse = "\n"))
  })
}

shinyApp(ui, server)
# nolint end
