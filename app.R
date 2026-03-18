# RNAPeaks Shiny Application
APP_VERSION <- "1.2.0"

# Null coalescing operator (explicit definition for clarity)
`%||%` <- function(x, y) if (is.null(x)) y else x

# Load sample data
load("data/sample_bed.rda")
load("data/sample_se.mats.rda")

# GTF is loaded lazily to reduce startup memory
gtf_human_path_global <- "data/gtf_human.rda"

library(shiny)
library(bslib)
library(bsicons)

# THEME SETUP
light_theme <- bs_theme(
  version = 5,
bootswatch = "flatly",
  bg = "#f8fafc",
  fg = "#0f172a",
  info = "#4f46e5",
  secondary = "#64748b",
  success = "#10b981",
  warning = "#f59e0b",
  danger = "#ef4444",
  base_font = font_google("IBM Plex Sans"),
  heading_font = font_google("IBM Plex Sans"),
  font_scale = 0.98,
  radius = "1rem"
)


# UI COMPONENTS
navbar_ui <- function(current_page) {
  div(
    class = "main-navbar",
    tags$a(class = "navbar-brand", onclick = "Shiny.setInputValue('nav_to', 'home', {priority: 'event'})", "RNAPeaks"),
    div(
      class = "navbar-links",
      tags$a(class = paste("nav-link-item", if(current_page == "plot_gene") "active" else ""),
             onclick = "Shiny.setInputValue('nav_to', 'plot_gene', {priority: 'event'})", "Plot Gene"),
      tags$a(class = paste("nav-link-item", if(current_page == "plot_region") "active" else ""),
             onclick = "Shiny.setInputValue('nav_to', 'plot_region', {priority: 'event'})", "Plot Region"),
      tags$a(class = paste("nav-link-item", if(current_page == "splicing_map") "active" else ""),
             onclick = "Shiny.setInputValue('nav_to', 'splicing_map', {priority: 'event'})", "Splicing Map"),
      tags$a(class = paste("nav-link-item", if(current_page == "sequence_map") "active" else ""),
             onclick = "Shiny.setInputValue('nav_to', 'sequence_map', {priority: 'event'})", "Sequence Map")
    ),
    div()
  )
}

function_card <- function(id, icon, title, description) {
  div(class = "col-md-6 col-lg-3 mb-4",
    div(class = "card function-card h-100 shadow-sm",
      onclick = sprintf("Shiny.setInputValue('nav_to', '%s', {priority: 'event'})", id),
      div(class = "card-body text-center p-4",
        div(class = "card-icon", bs_icon(icon)),
        h5(class = "card-title fw-bold", title),
        p(class = "card-text text-muted", description)
      )
    )
  )
}

sidebar_section <- function(title, ...) {
  div(class = "sidebar-section", h6(title), ...)
}

alert_box <- function(type, title, message) {
  icon_name <- switch(type, info = "info-circle-fill", warning = "exclamation-triangle-fill",
                      error = "x-circle-fill", success = "check-circle-fill", "info-circle-fill")
  div(class = paste("alert-box", paste0("alert-", type)),
    div(class = "alert-box-icon", bs_icon(icon_name)),
    div(class = "alert-box-content",
      div(class = "alert-box-title", title),
      div(class = "alert-box-message", message)
    )
  )
}

# Helper function to format error messages for users
format_error_message <- function(e) {
  msg <- e$message
  # Make common errors more user-friendly
  if (grepl("cannot open", msg, ignore.case = TRUE)) {
    return("Unable to read the uploaded file. Please check the file format.")
  }
  if (grepl("subscript out of bounds", msg, ignore.case = TRUE)) {
    return("The data format appears incorrect. Please verify your input file structure.")
  }
  if (grepl("Gene .* not found|No gene found", msg, ignore.case = TRUE)) {
    return(msg)
  }
  if (grepl("BSgenome", msg, ignore.case = TRUE)) {
    return("BSgenome.Hsapiens.UCSC.hg38 package is required for sequence analysis. Please install it.")
  }
  # Return original message if no special handling
  msg
}


# HOME PAGE
home_page <- function() {
  tagList(
    div(class = "hero-section", h1("RNAPeaks"), p("Visualize RNA-Binding Protein Peaks on Gene Structures")),
    div(class = "container",
      h5(class = "text-center mb-4 text-muted", "Select a visualization tool to get started"),
      div(class = "row justify-content-center",
        function_card("plot_gene", "bar-chart-fill", "Plot Gene", "Visualize RBP peaks on a single gene structure"),
        function_card("plot_region", "map-fill", "Plot Region", "Visualize peaks across a genomic region"),
        function_card("splicing_map", "diagram-3-fill", "Splicing Map", "Analyze protein binding around splice junctions"),
        function_card("sequence_map", "file-earmark-text-fill", "Sequence Map", "Map sequence motif frequency")
      )
    )
  )
}


# PLOT GENE PAGE
plot_gene_page <- function() {
  tagList(
    navbar_ui("plot_gene"),
    div(class = "page-content",
      layout_sidebar(
        sidebar = sidebar(width = 320,
          # Generate button at top
          div(class = "generate-btn-section",
            actionButton("pg_generate", "Generate Plot", class = "btn btn-primary", icon = icon("play"))
          ),

          sidebar_section("Data Input",
            radioButtons("pg_bed_source", "BED File:", choices = c("Sample Data" = "sample", "Upload" = "upload"), selected = "sample", inline = TRUE),
            conditionalPanel("input.pg_bed_source == 'upload'",
              fileInput("pg_bed_file", NULL, accept = c(".bed", ".txt", ".tsv")),
              div(class = "upload-hint", "BED format: chr, start, end, name, score, strand")
            ),
            radioButtons("pg_gtf_source", "GTF Annotation:", choices = c("Human (bundled)" = "human", "Upload" = "upload"), selected = "human", inline = TRUE),
            conditionalPanel("input.pg_gtf_source == 'upload'",
              fileInput("pg_gtf_file", NULL, accept = c(".gtf", ".gtf.gz")),
              div(class = "upload-hint", "Ensembl/GENCODE GTF format")
            )
          ),

          sidebar_section("Gene",
            textInput("pg_gene_id", "Gene ID:", value = "GAPDH", placeholder = "e.g., GAPDH, TP53, BRCA1"),
            textInput("pg_tx_id", "Transcript ID (optional):", value = "", placeholder = "e.g., ENST00000229239")
          ),

          sidebar_section("Peak Options",
            selectInput("pg_order_by", "Order By:", choices = c("Count" = "Count", "Target Name" = "Target"), selected = "Count"),
            selectInput("pg_peak_col", "Peak Color:", choices = c("Purple" = "purple", "Blue" = "blue", "Red" = "red", "Green" = "darkgreen"), selected = "purple"),
            numericInput("pg_merge", "Merge Nearby Peaks (bp):", value = 0, min = 0, step = 10),
            numericInput("pg_max_proteins", "Max Proteins to Show:", value = 40, min = 1, max = 100, step = 1)
          ),

          sidebar_section("Styling",
            numericInput("pg_title_size", "Title Size:", value = 25, min = 8, max = 40, step = 1),
            numericInput("pg_label_size", "Label Size:", value = 5, min = 2, max = 12, step = 0.5),
            numericInput("pg_axis_breaks_n", "# Position Markers:", value = 5, min = 2, max = 20, step = 1),
            fluidRow(
              column(6, numericInput("pg_total_arrows", "Total Arrows:", value = 6, min = 1, max = 20, step = 1)),
              column(6, numericInput("pg_max_per_intron", "Max/Intron:", value = 2, min = 1, max = 10, step = 1))
            ),
            checkboxInput("pg_five_to_three", "Orient 5' to 3' (flip negative strand)", value = FALSE),
            checkboxInput("pg_show_junctions", "Show Exon/Intron Junction Lines", value = FALSE),
            conditionalPanel("input.pg_show_junctions",
              selectInput("pg_junction_color", "Junction Line Color:",
                choices = c("Gray" = "gray40", "Black" = "black", "Red" = "red", "Blue" = "blue", "Green" = "darkgreen", "Orange" = "orange"),
                selected = "gray40")
            )
          ),

          sidebar_section("Highlight Region",
            numericInput("pg_highlight_start", "Start Position:", value = NA, min = 1),
            numericInput("pg_highlight_stop", "End Position:", value = NA, min = 1),
            selectInput("pg_highlight_color", "Highlight Color:",
              choices = c("Red" = "red", "Green" = "green", "Yellow" = "yellow", "Pink" = "pink", "Blue" = "blue", "Orange" = "orange"),
              selected = "pink")
          )
        ),

        div(class = "plot-card",
          div(class = "plot-card-header",
            conditionalPanel("output.pg_plot_ready",
              div(class = "dropdown download-dropdown",
                tags$button(class = "btn btn-outline-secondary btn-sm dropdown-toggle", type = "button",
                  `data-bs-toggle` = "dropdown", bs_icon("download"), " Download"),
                tags$ul(class = "dropdown-menu dropdown-menu-end",
                  tags$li(downloadLink("pg_download_pdf", class = "dropdown-item", tagList(bs_icon("file-earmark-pdf"), " PDF"))),
                  tags$li(downloadLink("pg_download_csv", class = "dropdown-item", tagList(bs_icon("file-earmark-spreadsheet"), " CSV")))
                )
              )
            )
          ),
          div(class = "plot-card-body", uiOutput("pg_plot_message"), plotOutput("pg_plot", height = "600px"))
        )
      )
    )
  )
}


# PLOT REGION PAGE
plot_region_page <- function() {
  tagList(
    navbar_ui("plot_region"),
    div(class = "page-content",
      layout_sidebar(
        sidebar = sidebar(width = 320,
          div(class = "generate-btn-section",
            actionButton("pr_generate", "Generate Plot", class = "btn btn-primary", icon = icon("play"))
          ),

          sidebar_section("Data Input",
            radioButtons("pr_bed_source", "BED File:", choices = c("Sample Data" = "sample", "Upload" = "upload"), selected = "sample", inline = TRUE),
            conditionalPanel("input.pr_bed_source == 'upload'",
              fileInput("pr_bed_file", NULL, accept = c(".bed", ".txt", ".tsv"))
            ),
            radioButtons("pr_gtf_source", "GTF Annotation:", choices = c("Human (bundled)" = "human", "Upload" = "upload"), selected = "human", inline = TRUE),
            conditionalPanel("input.pr_gtf_source == 'upload'",
              fileInput("pr_gtf_file", NULL, accept = c(".gtf", ".gtf.gz"))
            )
          ),

          sidebar_section("Region",
            fluidRow(
              column(6, textInput("pr_chr", "Chr:", value = "12")),
              column(6, selectInput("pr_strand", "Strand:", choices = c("+" = "+", "-" = "-")))
            ),
            numericInput("pr_start", "Start:", value = 6534512, min = 1),
            numericInput("pr_end", "End:", value = 6538374, min = 1),
            textInput("pr_gene_id", "Gene ID (optional):", value = "", placeholder = "e.g., FTH1"),
            textInput("pr_tx_id", "Transcript ID (optional):", value = "", placeholder = "e.g., ENST00000620041")
          ),

          sidebar_section("Peak Options",
            selectInput("pr_order_by", "Order By:", choices = c("Count" = "Count", "Target Name" = "Target"), selected = "Count"),
            selectInput("pr_peak_col", "Peak Color:", choices = c("Purple" = "purple", "Blue" = "blue", "Red" = "red", "Green" = "darkgreen"), selected = "purple"),
            numericInput("pr_max_proteins", "Max Proteins to Show:", value = 40, min = 1, max = 100, step = 1),
            checkboxInput("pr_five_to_three", "Orient 5' to 3' (flip negative strand)", value = FALSE),
            checkboxInput("pr_show_junctions", "Show Exon/Intron Junction Lines", value = FALSE)
          ),

          sidebar_section("Styling",
            numericInput("pr_title_size", "Title Size:", value = 25, min = 8, max = 40, step = 1),
            numericInput("pr_label_size", "Label Size:", value = 5, min = 2, max = 12, step = 0.5),
            numericInput("pr_axis_breaks_n", "# Position Markers:", value = 5, min = 2, max = 20, step = 1),
            fluidRow(
              column(6, numericInput("pr_total_arrows", "Total Arrows:", value = 12, min = 1, max = 30, step = 1)),
              column(6, numericInput("pr_max_per_intron", "Max/Intron:", value = 5, min = 1, max = 15, step = 1))
            ),
            conditionalPanel("input.pr_show_junctions",
              selectInput("pr_junction_color", "Junction Line Color:",
                choices = c("Gray" = "gray40", "Black" = "black", "Red" = "red", "Blue" = "blue", "Green" = "darkgreen", "Orange" = "orange"),
                selected = "gray40")
            )
          ),

          sidebar_section("Highlight Region",
            numericInput("pr_highlight_start", "Start Position:", value = NA, min = 1),
            numericInput("pr_highlight_stop", "End Position:", value = NA, min = 1),
            selectInput("pr_highlight_color", "Highlight Color:",
              choices = c("Red" = "red", "Green" = "green", "Yellow" = "yellow", "Pink" = "pink", "Blue" = "blue", "Orange" = "orange"),
              selected = "pink")
          )
        ),

        div(class = "plot-card",
          div(class = "plot-card-header",
            conditionalPanel("output.pr_plot_ready",
              div(class = "dropdown download-dropdown",
                tags$button(class = "btn btn-outline-secondary btn-sm dropdown-toggle", type = "button",
                  `data-bs-toggle` = "dropdown", bs_icon("download"), " Download"),
                tags$ul(class = "dropdown-menu dropdown-menu-end",
                  tags$li(downloadLink("pr_download_pdf", class = "dropdown-item", tagList(bs_icon("file-earmark-pdf"), " PDF"))),
                  tags$li(downloadLink("pr_download_csv", class = "dropdown-item", tagList(bs_icon("file-earmark-spreadsheet"), " CSV")))
                )
              )
            )
          ),
          div(class = "plot-card-body", uiOutput("pr_plot_message"), plotOutput("pr_plot", height = "600px"))
        )
      )
    )
  )
}


# SPLICING MAP PAGE
splicing_map_page <- function() {
  tagList(
    navbar_ui("splicing_map"),
    div(class = "page-content",
      layout_sidebar(
        sidebar = sidebar(width = 320,
          div(class = "generate-btn-section",
            actionButton("sm_generate", "Generate Map", class = "btn btn-primary", icon = icon("play"))
          ),

          sidebar_section("Data Input",
            radioButtons("sm_bed_source", "BED File:", choices = c("Sample" = "sample", "Upload" = "upload"), selected = "sample", inline = TRUE),
            conditionalPanel("input.sm_bed_source == 'upload'",
              fileInput("sm_bed_file", NULL, accept = c(".bed", ".txt", ".tsv"))
            ),
            radioButtons("sm_mats_source", "SE.MATS File:", choices = c("Sample" = "sample", "Upload" = "upload"), selected = "sample", inline = TRUE),
            conditionalPanel("input.sm_mats_source == 'upload'",
              fileInput("sm_mats_file", NULL, accept = c(".txt", ".tsv", ".csv")),
              div(class = "upload-hint", "rMATS SE.MATS.JC.txt output")
            )
          ),

          sidebar_section("Region Width",
            fluidRow(
              column(6, numericInput("sm_exon_width", "Exon (bp):", value = 50, min = 10, max = 200)),
              column(6, numericInput("sm_intron_width", "Intron (bp):", value = 300, min = 50, max = 500))
            ),
            numericInput("sm_moving_avg", "Smoothing Window:", value = 50, min = 0, max = 100)
          ),

          sidebar_section("Event Filtering",
            numericInput("sm_p_value", "P-value Threshold:", value = 0.05, min = 0.001, max = 0.1, step = 0.01),
            numericInput("sm_retained_inc", "Retained IncLevel Diff:", value = 0.1, min = 0, max = 1, step = 0.05),
            numericInput("sm_excluded_inc", "Excluded IncLevel Diff:", value = -0.1, min = -1, max = 0, step = 0.05),
            checkboxGroupInput("sm_groups", "Groups:", choices = c("Retained", "Excluded", "Control"), selected = c("Retained", "Excluded", "Control"), inline = TRUE)
          ),

          sidebar_section("Control Sampling",
            numericInput("sm_control_mult", "Control Multiplier:", value = 2.0, min = 0.5, max = 5, step = 0.5),
            numericInput("sm_control_iter", "Iterations:", value = 20, min = 5, max = 100, step = 5)
          ),

          sidebar_section("Significance",
            numericInput("sm_z_threshold", "Z-score Threshold:", value = 1.96, min = 1, max = 3, step = 0.1),
            numericInput("sm_min_consec", "Min Consecutive:", value = 10, min = 1, max = 50, step = 1)
          ),

          sidebar_section("Appearance",
            textInput("sm_title", "Plot Title:", value = "", placeholder = "Leave blank for no title"),
            fluidRow(
              column(4, selectInput("sm_retained_col", "Retained:", choices = c("Blue" = "blue", "Navy" = "navy", "Teal" = "teal", "Purple" = "purple", "Red" = "red", "Black" = "black"), selected = "blue")),
              column(4, selectInput("sm_excluded_col", "Excluded:", choices = c("Red" = "red", "Orange" = "orange", "Pink" = "hotpink", "Blue" = "blue", "Black" = "black", "Green" = "darkgreen"), selected = "red")),
              column(4, selectInput("sm_control_col", "Control:", choices = c("Black" = "black", "Gray" = "gray50", "Dark Gray" = "gray30", "Blue" = "blue", "Green" = "darkgreen"), selected = "black"))
            ),
            fluidRow(
              column(6, numericInput("sm_line_width", "Line Width:", value = 0.8, min = 0.2, max = 3, step = 0.1)),
              column(6, numericInput("sm_axis_text_size", "Axis Text Size:", value = 11, min = 6, max = 20, step = 1))
            ),
            numericInput("sm_title_size", "Title Size:", value = 20, min = 8, max = 40, step = 1),
            selectInput("sm_exon_col", "Skipped Exon Color:",
              choices = c("Navy" = "navy", "Blue" = "blue", "Black" = "black", "Dark Green" = "darkgreen", "Purple" = "purple4", "Dark Gray" = "gray30"),
              selected = "navy")
          )
        ),

        div(class = "plot-card",
          div(class = "plot-card-header",
            conditionalPanel("output.sm_plot_ready",
              div(class = "dropdown download-dropdown",
                tags$button(class = "btn btn-outline-secondary btn-sm dropdown-toggle", type = "button",
                  `data-bs-toggle` = "dropdown", bs_icon("download"), " Download"),
                tags$ul(class = "dropdown-menu dropdown-menu-end",
                  tags$li(downloadLink("sm_download_pdf", class = "dropdown-item", tagList(bs_icon("file-earmark-pdf"), " PDF")))
                )
              )
            )
          ),
          div(class = "plot-card-body", uiOutput("sm_plot_ui"))
        )
      )
    )
  )
}


# SEQUENCE MAP PAGE
sequence_map_page <- function() {
  tagList(
    navbar_ui("sequence_map"),
    div(class = "page-content",
      layout_sidebar(
        sidebar = sidebar(width = 320,
          div(class = "generate-btn-section",
            actionButton("sqm_generate", "Generate Map", class = "btn btn-primary", icon = icon("play"))
          ),

          sidebar_section("Data Input",
            radioButtons("sqm_mats_source", "SE.MATS File:", choices = c("Sample" = "sample", "Upload" = "upload"), selected = "sample", inline = TRUE),
            conditionalPanel("input.sqm_mats_source == 'upload'",
              fileInput("sqm_mats_file", NULL, accept = c(".txt", ".tsv", ".csv")),
              div(class = "upload-hint", "rMATS SE.MATS.JC.txt output")
            )
          ),

          sidebar_section("Sequence",
            textInput("sqm_sequence", "Motif:", value = "YGCY", placeholder = "e.g., YGCY, GCAUG"),
            div(class = "upload-hint", "Supports IUPAC: Y=C/T, R=A/G, N=any")
          ),

          sidebar_section("Region Width",
            fluidRow(
              column(6, numericInput("sqm_exon_width", "Exon (bp):", value = 50, min = 10, max = 200)),
              column(6, numericInput("sqm_intron_width", "Intron (bp):", value = 250, min = 50, max = 500))
            ),
            numericInput("sqm_moving_avg", "Smoothing Window:", value = 40, min = 0, max = 100)
          ),

          sidebar_section("Event Filtering",
            numericInput("sqm_p_value", "P-value Threshold:", value = 0.05, min = 0.001, max = 0.1, step = 0.01),
            numericInput("sqm_retained_inc", "Retained IncLevel Diff:", value = 0.1, min = 0, max = 1, step = 0.05),
            numericInput("sqm_excluded_inc", "Excluded IncLevel Diff:", value = -0.1, min = -1, max = 0, step = 0.05),
            checkboxGroupInput("sqm_groups", "Groups:", choices = c("Retained", "Excluded", "Control"), selected = c("Retained", "Excluded", "Control"), inline = TRUE)
          ),

          sidebar_section("Control Sampling",
            numericInput("sqm_control_mult", "Control Multiplier:", value = 2.0, min = 0.5, max = 5, step = 0.5),
            numericInput("sqm_control_iter", "Iterations:", value = 20, min = 5, max = 100, step = 5)
          ),

          sidebar_section("Significance",
            numericInput("sqm_z_threshold", "Z-score Threshold:", value = 1.96, min = 1, max = 3, step = 0.1),
            numericInput("sqm_min_consec", "Min Consecutive:", value = 10, min = 1, max = 50, step = 1)
          ),

          sidebar_section("Appearance",
            textInput("sqm_title", "Plot Title:", value = "", placeholder = "Leave blank for no title"),
            fluidRow(
              column(4, selectInput("sqm_retained_col", "Retained:", choices = c("Blue" = "blue", "Navy" = "navy", "Teal" = "teal", "Purple" = "purple", "Red" = "red", "Black" = "black"), selected = "blue")),
              column(4, selectInput("sqm_excluded_col", "Excluded:", choices = c("Red" = "red", "Orange" = "orange", "Pink" = "hotpink", "Blue" = "blue", "Black" = "black", "Green" = "darkgreen"), selected = "red")),
              column(4, selectInput("sqm_control_col", "Control:", choices = c("Black" = "black", "Gray" = "gray50", "Dark Gray" = "gray30", "Blue" = "blue", "Green" = "darkgreen"), selected = "black"))
            ),
            fluidRow(
              column(6, numericInput("sqm_line_width", "Line Width:", value = 0.8, min = 0.2, max = 3, step = 0.1)),
              column(6, numericInput("sqm_axis_text_size", "Axis Text Size:", value = 11, min = 6, max = 20, step = 1))
            ),
            numericInput("sqm_title_size", "Title Size:", value = 20, min = 8, max = 40, step = 1),
            selectInput("sqm_exon_col", "Skipped Exon Color:",
              choices = c("Navy" = "navy", "Blue" = "blue", "Black" = "black", "Dark Green" = "darkgreen", "Purple" = "purple4", "Dark Gray" = "gray30"),
              selected = "navy")
          )
        ),

        div(class = "plot-card",
          div(class = "plot-card-header",
            conditionalPanel("output.sqm_plot_ready",
              div(class = "dropdown download-dropdown",
                tags$button(class = "btn btn-outline-secondary btn-sm dropdown-toggle", type = "button",
                  `data-bs-toggle` = "dropdown", bs_icon("download"), " Download"),
                tags$ul(class = "dropdown-menu dropdown-menu-end",
                  tags$li(downloadLink("sqm_download_pdf", class = "dropdown-item", tagList(bs_icon("file-earmark-pdf"), " PDF")))
                )
              )
            )
          ),
          div(class = "plot-card-body", uiOutput("sqm_plot_ui"))
        )
      )
    )
  )
}


# MAIN UI
ui <- page_fluid(
  theme = light_theme,
  tags$head(tags$link(rel = "stylesheet", type = "text/css", href = "styles.css")),
  uiOutput("main_content")
)


# SERVER
server <- function(input, output, session) {

  # NAVIGATION
  current_page <- reactiveVal("home")
  observeEvent(input$nav_to, { current_page(input$nav_to) })

  output$main_content <- renderUI({
    switch(current_page(),
      "home" = home_page(),
      "plot_gene" = plot_gene_page(),
      "plot_region" = plot_region_page(),
      "splicing_map" = splicing_map_page(),
      "sequence_map" = sequence_map_page(),
      home_page()
    )
  })

  # REACTIVE VALUES
  rv <- reactiveValues(
    pg_result = NULL, pg_generated = FALSE, pg_error = NULL,
    pr_result = NULL, pr_generated = FALSE, pr_error = NULL,
    sm_result = NULL, sm_generated = FALSE, sm_error = NULL,
    sqm_result = NULL, sqm_generated = FALSE, sqm_error = NULL,
    gtf_human_cached = NULL  # Lazy-loaded GTF cache
  )

  # Lazy load GTF only when needed
  get_gtf_human <- function() {
    if (is.null(rv$gtf_human_cached)) {
      showNotification("Loading GTF annotation (first time only)...", type = "message", duration = NULL, id = "gtf_loading")
      load(gtf_human_path_global, envir = environment())
      rv$gtf_human_cached <- gtf_human
      removeNotification("gtf_loading")
      showNotification("GTF loaded successfully!", type = "message", duration = 2)
    }
    rv$gtf_human_cached
  }

  # PLOT GENE
  pg_bed_data <- reactive({
    if (input$pg_bed_source == "sample") sample_bed
    else { req(input$pg_bed_file); read.table(input$pg_bed_file$datapath, header = FALSE, sep = "\t", stringsAsFactors = FALSE) }
  })

  pg_gtf_data <- reactive({
    if (input$pg_gtf_source == "human") get_gtf_human()
    else { req(input$pg_gtf_file); rtracklayer::import(input$pg_gtf_file$datapath) }
  })

  observeEvent(input$pg_generate, {
    rv$pg_error <- NULL
    rv$pg_generated <- FALSE

    # Validation
    if (input$pg_bed_source == "upload" && is.null(input$pg_bed_file)) {
      rv$pg_error <- list(title = "Missing BED File", message = "Please upload a BED file or select 'Sample Data'.")
      return()
    }
    if (is.null(input$pg_gene_id) || input$pg_gene_id == "") {
      rv$pg_error <- list(title = "Missing Gene ID", message = "Please enter a gene symbol or Ensembl ID.")
      return()
    }

    showNotification("Generating plot...", type = "message", duration = NULL, id = "pg_loading")

    tryCatch({
      # Handle optional TxID
      tx_id <- if (is.null(input$pg_tx_id) || input$pg_tx_id == "") NA else input$pg_tx_id

      # Handle optional highlight region
      highlight_start <- if (is.na(input$pg_highlight_start)) NULL else input$pg_highlight_start
      highlight_stop <- if (is.na(input$pg_highlight_stop)) NULL else input$pg_highlight_stop

      result <- PlotGene(
        bed = pg_bed_data(),
        geneID = input$pg_gene_id,
        gtf = pg_gtf_data(),
        species = "Human",
        TxID = tx_id,
        order_by = input$pg_order_by,
        peak_col = input$pg_peak_col,
        merge = input$pg_merge,
        title_size = input$pg_title_size,
        label_size = input$pg_label_size,
        max_proteins = input$pg_max_proteins,
        total_arrows = input$pg_total_arrows,
        max_per_intron = input$pg_max_per_intron,
        five_to_three = input$pg_five_to_three,
        show_junctions = isTRUE(input$pg_show_junctions),
        junction_color = if (is.null(input$pg_junction_color)) "gray40" else input$pg_junction_color,
        axis_breaks_n  = if (is.null(input$pg_axis_breaks_n) || is.na(input$pg_axis_breaks_n)) 5L else as.integer(input$pg_axis_breaks_n),
        highlighted_region_start = highlight_start,
        highlighted_region_stop = highlight_stop,
        highlighted_region_color = input$pg_highlight_color,
        RNA_Peaks_File_Path = tempfile(fileext = ".pdf"),
        Bed_File_Path = tempfile(fileext = ".csv")
      )
      rv$pg_result <- result
      rv$pg_generated <- TRUE
      removeNotification("pg_loading")
      showNotification("Plot generated successfully!", type = "message", duration = 3)
    }, error = function(e) {
      removeNotification("pg_loading")
      rv$pg_error <- list(title = "Generation Failed", message = e$message)
    })
  })

  output$pg_plot <- renderPlot({ req(rv$pg_result); rv$pg_result$plot })
  output$pg_plot_ready <- reactive({ rv$pg_generated })
  outputOptions(output, "pg_plot_ready", suspendWhenHidden = FALSE)

  output$pg_plot_message <- renderUI({
    if (!is.null(rv$pg_error)) {
      alert_box("error", rv$pg_error$title, rv$pg_error$message)
    } else if (!rv$pg_generated) {
      div(class = "text-center text-muted py-5",
        bs_icon("bar-chart", size = "4rem", class = "mb-3 opacity-25"),
        h5("Ready to visualize"),
        p("Enter a gene ID and click 'Generate Plot'")
      )
    }
  })

  output$pg_download_pdf <- downloadHandler(
    filename = function() paste0("RNAPeaks_Gene_", input$pg_gene_id, "_", Sys.Date(), ".pdf"),
    content = function(file) { req(rv$pg_result); ggplot2::ggsave(file, rv$pg_result$plot, width = 16, height = 12, device = "pdf") }
  )
  output$pg_download_csv <- downloadHandler(
    filename = function() paste0("RNAPeaks_Gene_", input$pg_gene_id, "_", Sys.Date(), ".csv"),
    content = function(file) { req(rv$pg_result); write.csv(rv$pg_result$csv, file, row.names = FALSE) }
  )


  # PLOT REGION
  pr_bed_data <- reactive({
    if (input$pr_bed_source == "sample") sample_bed
    else { req(input$pr_bed_file); read.table(input$pr_bed_file$datapath, header = FALSE, sep = "\t", stringsAsFactors = FALSE) }
  })

  pr_gtf_data <- reactive({
    if (input$pr_gtf_source == "human") get_gtf_human()
    else { req(input$pr_gtf_file); rtracklayer::import(input$pr_gtf_file$datapath) }
  })

  observeEvent(input$pr_generate, {
    rv$pr_error <- NULL
    rv$pr_generated <- FALSE

    if (input$pr_bed_source == "upload" && is.null(input$pr_bed_file)) {
      rv$pr_error <- list(title = "Missing BED File", message = "Please upload a BED file or select 'Sample Data'.")
      return()
    }
    if (input$pr_start >= input$pr_end) {
      rv$pr_error <- list(title = "Invalid Region", message = "Start position must be less than end position.")
      return()
    }

    showNotification("Generating plot...", type = "message", duration = NULL, id = "pr_loading")

    tryCatch({
      # Handle optional highlight region
      highlight_start <- if (is.na(input$pr_highlight_start)) NULL else input$pr_highlight_start
      highlight_stop <- if (is.na(input$pr_highlight_stop)) NULL else input$pr_highlight_stop

      pr_gene_id <- if (is.null(input$pr_gene_id) || input$pr_gene_id == "") NULL else input$pr_gene_id
      pr_tx_id   <- if (is.null(input$pr_tx_id)   || input$pr_tx_id   == "") NA   else input$pr_tx_id

      result <- PlotRegion(
        bed = pr_bed_data(),
        Chr = input$pr_chr,
        Start = input$pr_start,
        End = input$pr_end,
        Strand = input$pr_strand,
        gtf = pr_gtf_data(),
        geneID = pr_gene_id,
        TxID = pr_tx_id,
        species = "Human",
        order_by = input$pr_order_by,
        peak_col = input$pr_peak_col,
        title_size = input$pr_title_size,
        label_size = input$pr_label_size,
        max_proteins = input$pr_max_proteins,
        five_to_three = input$pr_five_to_three,
        show_junctions = isTRUE(input$pr_show_junctions),
        junction_color = if (is.null(input$pr_junction_color)) "gray40" else input$pr_junction_color,
        axis_breaks_n  = if (is.null(input$pr_axis_breaks_n) || is.na(input$pr_axis_breaks_n)) 5L else as.integer(input$pr_axis_breaks_n),
        total_arrows = input$pr_total_arrows,
        max_per_intron = input$pr_max_per_intron,
        highlighted_region_start = highlight_start,
        highlighted_region_stop = highlight_stop,
        highlighted_region_color = input$pr_highlight_color,
        RNA_Peaks_File_Path = tempfile(fileext = ".pdf"),
        Bed_File_Path = tempfile(fileext = ".csv")
      )
      rv$pr_result <- result
      rv$pr_generated <- TRUE
      removeNotification("pr_loading")
      showNotification("Plot generated successfully!", type = "message", duration = 3)
    }, error = function(e) {
      removeNotification("pr_loading")
      rv$pr_error <- list(title = "Generation Failed", message = e$message)
    })
  })

  output$pr_plot <- renderPlot({ req(rv$pr_result); rv$pr_result$plot })
  output$pr_plot_ready <- reactive({ rv$pr_generated })
  outputOptions(output, "pr_plot_ready", suspendWhenHidden = FALSE)

  output$pr_plot_message <- renderUI({
    if (!is.null(rv$pr_error)) {
      alert_box("error", rv$pr_error$title, rv$pr_error$message)
    } else if (!rv$pr_generated) {
      div(class = "text-center text-muted py-5",
        bs_icon("map", size = "4rem", class = "mb-3 opacity-25"),
        h5("Ready to visualize"),
        p("Enter region coordinates and click 'Generate Plot'")
      )
    }
  })

  output$pr_download_pdf <- downloadHandler(
    filename = function() paste0("RNAPeaks_Region_chr", input$pr_chr, "_", Sys.Date(), ".pdf"),
    content = function(file) { req(rv$pr_result); ggplot2::ggsave(file, rv$pr_result$plot, width = 16, height = 12, device = "pdf") }
  )
  output$pr_download_csv <- downloadHandler(
    filename = function() paste0("RNAPeaks_Region_chr", input$pr_chr, "_", Sys.Date(), ".csv"),
    content = function(file) { req(rv$pr_result); write.csv(rv$pr_result$csv, file, row.names = FALSE) }
  )

  # SPLICING MAP
  sm_bed_data <- reactive({
    if (input$sm_bed_source == "sample") sample_bed
    else { req(input$sm_bed_file); read.table(input$sm_bed_file$datapath, header = FALSE, sep = "\t", stringsAsFactors = FALSE) }
  })

  sm_mats_data <- reactive({
    if (input$sm_mats_source == "sample") sample_se.mats
    else { req(input$sm_mats_file); read.table(input$sm_mats_file$datapath, header = TRUE, sep = "\t", stringsAsFactors = FALSE) }
  })

  observeEvent(input$sm_generate, {
    rv$sm_error <- NULL
    rv$sm_generated <- FALSE

    if (input$sm_bed_source == "upload" && is.null(input$sm_bed_file)) {
      rv$sm_error <- list(title = "Missing BED File", message = "Please upload a BED file or select 'Sample'.")
      return()
    }
    if (input$sm_mats_source == "upload" && is.null(input$sm_mats_file)) {
      rv$sm_error <- list(title = "Missing SE.MATS File", message = "Please upload an SE.MATS file or select 'Sample'.")
      return()
    }

    tryCatch({
      withProgress(message = "Generating Splicing Map", value = 0, {
        last_progress <- 0

        incProgress(0.05, detail = "Loading data...")
        bed_data <- sm_bed_data()
        mats_data <- sm_mats_data()

        incProgress(0.05, detail = "Starting analysis...")
        last_progress <- 0.10

        result <- createSplicingMap(
          bed_file = bed_data,
          SEMATS = mats_data,
          WidthIntoExon = input$sm_exon_width,
          WidthIntoIntron = input$sm_intron_width,
          moving_average = input$sm_moving_avg,
          p_valueRetainedAndExclusion = input$sm_p_value,
          retained_IncLevelDifference = input$sm_retained_inc,
          exclusion_IncLevelDifference = input$sm_excluded_inc,
          groups = input$sm_groups,
          control_multiplier = input$sm_control_mult,
          control_iterations = input$sm_control_iter,
          z_threshold = input$sm_z_threshold,
          min_consecutive = input$sm_min_consec,
          title = input$sm_title,
          retained_col = input$sm_retained_col,
          excluded_col = input$sm_excluded_col,
          control_col = input$sm_control_col,
          line_width = input$sm_line_width,
          axis_text_size = input$sm_axis_text_size,
          title_size = input$sm_title_size,
          exon_col = input$sm_exon_col,
          progress_callback = function(current, total, detail = NULL) {
            # createSplicingMap reports 0..100
            target <- max(0, min(1, current / total))

            # Ensure monotonic progress and avoid negative increments
            if (target > last_progress) {
              incProgress(target - last_progress, detail = detail %||% "Processing...")
              last_progress <<- target
            }
          }
        )

        incProgress(0.4, detail = "Complete!")
        rv$sm_result <- result
        rv$sm_generated <- TRUE
      })
      showNotification("Splicing map generated!", type = "message", duration = 3)
    }, error = function(e) {
      rv$sm_error <- list(title = "Generation Failed", message = e$message)
    })
  })

  output$sm_plot <- renderPlot({ req(rv$sm_result); rv$sm_result })
  output$sm_plot_ready <- reactive({ rv$sm_generated })
  outputOptions(output, "sm_plot_ready", suspendWhenHidden = FALSE)

  output$sm_plot_ui <- renderUI({
    if (!is.null(rv$sm_error)) {
      alert_box("error", rv$sm_error$title, rv$sm_error$message)
    } else if (rv$sm_generated) {
      plotOutput("sm_plot", height = "600px")
    } else {
      div(class = "text-center text-muted py-5",
        bs_icon("diagram-3", size = "4rem", class = "mb-3 opacity-25"),
        h5("Ready to analyze"),
        p("Upload data and click 'Generate Map'")
      )
    }
  })

  output$sm_download_pdf <- downloadHandler(
    filename = function() paste0("RNAPeaks_SplicingMap_", Sys.Date(), ".pdf"),
    content = function(file) { req(rv$sm_result); ggplot2::ggsave(file, rv$sm_result, width = 12, height = 8, device = "pdf") }
  )

  # SEQUENCE MAP
  sqm_mats_data <- reactive({
    if (input$sqm_mats_source == "sample") sample_se.mats
    else { req(input$sqm_mats_file); read.table(input$sqm_mats_file$datapath, header = TRUE, sep = "\t", stringsAsFactors = FALSE) }
  })

  observeEvent(input$sqm_generate, {
    rv$sqm_error <- NULL
    rv$sqm_generated <- FALSE

    if (input$sqm_mats_source == "upload" && is.null(input$sqm_mats_file)) {
      rv$sqm_error <- list(title = "Missing SE.MATS File", message = "Please upload an SE.MATS file or select 'Sample'.")
      return()
    }
    if (is.null(input$sqm_sequence) || input$sqm_sequence == "") {
      rv$sqm_error <- list(title = "Missing Sequence", message = "Please enter a sequence motif (e.g., YGCY).")
      return()
    }

    tryCatch({
      withProgress(message = "Generating Sequence Map", value = 0, {
        incProgress(0.05, detail = "Loading data...")
        mats_data <- sqm_mats_data()

        incProgress(0.05, detail = "Starting analysis...")
        last_progress <- 0.10

        result <- createSequenceMap(
          SEMATS = mats_data,
          sequence = input$sqm_sequence,
          WidthIntoExon = input$sqm_exon_width,
          WidthIntoIntron = input$sqm_intron_width,
          moving_average = input$sqm_moving_avg,
          p_valueRetainedAndExclusion = input$sqm_p_value,
          retained_IncLevelDifference = input$sqm_retained_inc,
          exclusion_IncLevelDifference = input$sqm_excluded_inc,
          groups = input$sqm_groups,
          control_multiplier = input$sqm_control_mult,
          control_iterations = input$sqm_control_iter,
          z_threshold = input$sqm_z_threshold,
          min_consecutive = input$sqm_min_consec,
          title = input$sqm_title,
          retained_col = input$sqm_retained_col,
          excluded_col = input$sqm_excluded_col,
          control_col = input$sqm_control_col,
          line_width = input$sqm_line_width,
          axis_text_size = input$sqm_axis_text_size,
          title_size = input$sqm_title_size,
          exon_col = input$sqm_exon_col,
          progress_callback = function(current, total, detail = NULL) {
            # createSequenceMap reports 0..100
            target <- max(0, min(1, current / total))

            # Ensure monotonic progress and avoid negative increments
            if (target > last_progress) {
              incProgress(target - last_progress, detail = detail %||% "Processing...")
              last_progress <<- target
            }
          }
        )

        incProgress(0.4, detail = "Complete!")
        rv$sqm_result <- result
        rv$sqm_generated <- TRUE
      })
      showNotification("Sequence map generated!", type = "message", duration = 3)
    }, error = function(e) {
      rv$sqm_error <- list(title = "Generation Failed", message = e$message)
    })
  })

  output$sqm_plot <- renderPlot({ req(rv$sqm_result); rv$sqm_result })
  output$sqm_plot_ready <- reactive({ rv$sqm_generated })
  outputOptions(output, "sqm_plot_ready", suspendWhenHidden = FALSE)

  output$sqm_plot_ui <- renderUI({
    if (!is.null(rv$sqm_error)) {
      alert_box("error", rv$sqm_error$title, rv$sqm_error$message)
    } else if (rv$sqm_generated) {
      plotOutput("sqm_plot", height = "600px")
    } else {
      div(class = "text-center text-muted py-5",
        bs_icon("file-earmark-text", size = "4rem", class = "mb-3 opacity-25"),
        h5("Ready to analyze"),
        p("Enter a sequence motif and click 'Generate Map'")
      )
    }
  })

  output$sqm_download_pdf <- downloadHandler(
    filename = function() paste0("RNAPeaks_SequenceMap_", input$sqm_sequence, "_", Sys.Date(), ".pdf"),
    content = function(file) { req(rv$sqm_result); ggplot2::ggsave(file, rv$sqm_result, width = 12, height = 8, device = "pdf") }
  )
}

shinyApp(ui = ui, server = server)
