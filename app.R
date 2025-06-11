# input packages ---- 
library(shiny)
library(ggplot2)
library(ggsignif)
library(ggpubr)
library(tidyr)
library(dplyr)
library(readxl)
library(DT)
options(shiny.maxRequestSize=30*1024^2)
# 定义 UI ----
ui <- fluidPage(
  titlePanel("Plasma Protein Contamination Correction and Differential Expression Analysis"),
  tabsetPanel(id = "Step",
              ## Welcome Tab ----
              tabPanel("Welcome",
                       # 外层容器：确保内容整体居中
                       div(style = "display: flex; justify-content: center; width: 100%;",
                           # 内容区域：固定最大宽度并左对齐
                           div(style = "max-width: 900px; width: 100%; text-align: left; padding: 20px;",
                               
                               # 标题（居中显示）
                               h2("Welcome to Plasma Protein Analysis Tool", 
                                  style = "text-align: center; font-size: calc(20px + 1vw); margin-bottom: 25px;"),
                               
                               # 图片（居中显示但内容左对齐）
                               div(style = "text-align: center; margin: 20px 0;",
                                   img(src = "Welcome_2.png", 
                                       style = "max-width: 100%; height: auto; border-radius: 8px;")
                               ),
                               
                               # 正文（强制左对齐）
                               div(style = "text-align: left;",  # 显式覆盖可能继承的居中样式
                                   p("This interactive tool allows you to analyze and correct for contamination in plasma proteomics data."),
                                   
                                   h4("Key features:", style = "margin-top: 25px;"),
                                   tags$ul(
                                     style = "padding-left: 20px;",
                                     tags$li("Data quality assessment and visualization"),
                                     tags$li("Contamination marker evaluation"),
                                     tags$li("Correction for erythrocyte, platelet and coagulation contamination"),
                                     tags$li("Differential expression analysis"),
                                     tags$li("Interactive visualization of results")
                                   ),
                                   
                                   h4("How to use:", style = "margin-top: 25px;"),
                                   tags$ol(
                                     style = "padding-left: 20px;",
                                     tags$li("Upload your protein expression data and group information"),
                                     tags$li("Check data quality and select contamination markers"),
                                     tags$li("Run correction for selected contamination types"),
                                     tags$li("Perform differential expression analysis")
                                   )
                               ),
                               
                               # 页脚（居中显示）
                               hr(style = "margin: 30px 0; border-top: 1px solid #eee;"),
                               p(style = "text-align: center; font-style: italic;", 
                                 "For questions or feedback, please contact support@example.com")
                           )
                       )
              ),
              ## Step1 ----
              # 改为 data input
              tabPanel("Step 1: Data Input",
                       sidebarLayout(
                         sidebarPanel(h3("Step 1: Data Input",
                                         tags$span(
                                           id = 'span1',
                                           `data-toggle` = "tooltip",
                                           title = 'In this part, users can upload their own proteomics expression data and sample group data. 
                                           The example data can be found when users click "Load example data" below. 
                                           Detailed descriptions are provided in the "Help" part.',
                                           tags$span(class = "glyphicon glyphicon-question-sign")
                                         )),
                                      radioButtons("data_source", "Select Data Source",
                                                   choices = list("Load experimental data" = "experimental", 
                                                                  "Load example data" = "example"),
                                                   selected = "experimental"),
                                      conditionalPanel(
                                        condition = "input.data_source == 'experimental'",
                                        fileInput("data_file", "Upload Data File (CSV)", 
                                                  accept = ".csv"),
                                        fileInput("group_file", "Upload Group Info File (CSV)", 
                                                  accept = ".csv")
                                      ),
                                      h3("Grouping settings"),
                                      # 新增去除生物学差异选项
                                      checkboxInput("remove_biological_diff", 
                                                    "Remove proteins with biological differences between groups", 
                                                    value = TRUE),
                                      # 条件面板控制分组设置显示
                                      conditionalPanel(
                                        condition = "input.remove_biological_diff == true",
                                        selectInput("group1", "Group 1", choices = NULL),
                                        selectInput("group2", "Group 2", choices = NULL)
                                      ),
                                      sliderInput("cor_cutoff", "Correlation Cutoff",
                                                  min = 0.5, max = 0.99, value = 0.9, step = 0.01),
                                      actionButton("run_check", "Run Data Check")
                         ),
                         mainPanel(conditionalPanel(
                           condition = "output.data_loaded == false",
                           h4("No data uploaded. Please upload your data or use example data to explore.")
                         ),
                         conditionalPanel(
                           condition = "output.data_loaded == true",
                           h4("Data File Preview"),
                           DTOutput("data_table"), 
                           h4("Group Information Preview"),
                           DTOutput("group_table")
                         )))
              ),
              ## Step2 ----
              # 拆分为data evaluation和define markers
              tabPanel("Step 2: Check markers and contamination levels",
                       sidebarLayout(
                         sidebarPanel(
                           sliderInput("cor_cutoff_step2", "Correlation Cutoff",
                                       min = 0.5, max = 0.99, value = 0.9, step = 0.01),
                           actionButton("rerun_check_step2", "Re-run Data Check"),
                           h3("Step 2: Contamination Assessment"),
                           # 修改后的单选按钮组
                           tags$div(class = "form-group",
                                    tags$label(class = "control-label", "Correction Type"),
                                    checkboxInput("type_all", "All", value = FALSE),
                                    uiOutput("erythrocyte_checkbox"),
                                    uiOutput("platelet_checkbox"),
                                    uiOutput("coagulation_checkbox")),
                           sliderInput("constraint_factor", "constraint factor",
                                       min = 0.5, max = 1.5, value = 1.0, step = 0.01),
                           actionButton("run_correct", "Run Correction")),
                         mainPanel(h3("Data Quality Assessment"),
                                   h4("Contamination Summary"),
                                   verbatimTextOutput("contamination_summary"),
                                   tabsetPanel(
                                     tabPanel("Erythrocyte", DTOutput("erythrocyte_marker_table")),
                                     tabPanel("Coagulation", DTOutput("coagulation_marker_table")),
                                     tabPanel("Platelet", DTOutput("platelet_marker_table"))
                                   ),
                                   h4("Quality Control"),
                                   tabsetPanel(
                                     tabPanel("correlation", fluidRow(
                                       column(width = 5,
                                              plotOutput("correlation_p_pre_plot")
                                       ),
                                       column(width = 5,
                                              plotOutput("correlation_r_pre_plot")))),
                                     tabPanel("PCA", plotOutput("pca_pre_plot")),
                                     tabPanel("Heatmap", plotOutput("heatmap_pre_plot")),
                                     tabPanel("Boxplot", plotOutput("boxplot_pre_plot"))
                                   ),
                                   h4("Contamination Marker Expression"),
                                   tabsetPanel(
                                     tabPanel("Erythrocyte",
                                              DTOutput("data_marker_erythrocyte"),
                                              plotOutput("contamination_erythrocyte_plot")),
                                     tabPanel("Coagulation", 
                                              DTOutput("data_marker_coagulation"),
                                              plotOutput("contamination_coagulation_plot")),
                                     tabPanel("Platelet", 
                                              DTOutput("data_marker_platelet"),
                                              plotOutput("contamination_platelet_plot"))
                                   ),
                                   h4("污染marker的相关性"),
                                   tabsetPanel(
                                     tabPanel("Erythrocyte",
                                              downloadButton("download_cor_data_erythrocyte", "Download Corrected Data"),
                                              DTOutput("cor_erythrocyte_data"),
                                              plotOutput("cor_erythrocyte_plot")),
                                     tabPanel("Coagulation", 
                                              downloadButton("download_cor_data_coagulation", "Download Corrected Data"),
                                              DTOutput("cor_coagulation_data"),
                                              plotOutput("cor_coagulation_plot")),
                                     tabPanel("Platelet", 
                                              downloadButton("download_cor_data_platelet", "Download Corrected Data"),
                                              DTOutput("cor_platelet_data"),
                                              plotOutput("cor_platelet_plot"))
                                   ),
                                   h4("Contamination Levels"),
                                   tabsetPanel(
                                     tabPanel("CV Analysis", plotOutput("cv_pre_plot")),
                                     tabPanel("Erythrocyte", plotOutput("erythrocyte_marker_pre_plot")),
                                     tabPanel("Coagulation", plotOutput("coagulation_marker_pre_plot")),
                                     tabPanel("Platelet", plotOutput("platelet_marker_pre_plot"))
                                   )
                         )
                       )
              ),
              ## Step3 ----
              # correction
              tabPanel("Step 3: Correction Results",
                       sidebarLayout(
                         sidebarPanel(h3("Step 3: Correction Analysis"),
                                      sliderInput("constraint_factor_step2", "constraint factor",
                                                  min = 0.5, max = 1.5, value = 1.0, step = 0.01),
                                      actionButton("run_correct_step2", "Re-run Correction"),
                                      h3("Step 4: DE Analysis"),
                                      actionButton("run_de", "Run Differential Expression Analysis")
                         ),
                         mainPanel(h3("Correction Outcomes"), 
                                   h4("Post-correction QC"),
                                   tabsetPanel(
                                     tabPanel("correlation", fluidRow(
                                       column(width = 5,
                                              plotOutput("correlation_p_post_plot")
                                       ),
                                       column(width = 5,
                                              plotOutput("correlation_r_post_plot")))),
                                     tabPanel("PCA", plotOutput("pca_post_plot")),
                                     tabPanel("Heatmap", plotOutput("heatmap_post_plot")),
                                     tabPanel("Boxplot", plotOutput("boxplot_post_plot"))
                                   ),
                                   h4("Post-correction Contamination"),
                                   tabsetPanel(
                                     tabPanel("CV Analysis", plotOutput("cv_post_plot")),
                                     tabPanel("Erythrocyte", plotOutput("erythrocyte_marker_post_plot")),
                                     tabPanel("Coagulation", plotOutput("coagulation_marker_post_plot")),
                                     tabPanel("Platelet", plotOutput("platelet_marker_post_plot"))
                                   ),
                                   h4("Corrected Data Matrix"),
                                   downloadButton("download_data", "Download Corrected Data"),
                                   DTOutput("data_correct_table")
                         )
                       )),
              ## Step4 DE ----
              # DE & enrichment
              tabPanel("Step 4: Differential Expression",
                       sidebarLayout(
                         sidebarPanel(
                           h3("Step 4: DE Analysis"),
                           # Add conditional message when remove_biological_diff is FALSE
                           conditionalPanel(
                             condition = "input.remove_biological_diff == false",
                             div(style = "color: red; font-weight: bold;",
                                 "Please upload grouping file and check 'Remove proteins with biological differences' in Step 1 to perform differential expression analysis.")
                           )
                         ),
                         mainPanel(tabsetPanel(
                           tabPanel("Raw Data Results",
                                    downloadButton("download_de_pre", "Download Raw DE Results"),
                                    DTOutput("result_de_pre_table")),
                           tabPanel("Corrected Results",
                                    downloadButton("download_de_post", "Download Corrected DE Results"),
                                    DTOutput("result_de_post_table")
                           )),
                           tabsetPanel(
                             tabPanel("Venn Diagram", plotOutput("vn_plot")),
                             tabPanel("Volcano Plot (Raw)", plotOutput("volc_de_pre")),
                             tabPanel("Volcano Plot (Corrected)", plotOutput("volc_de_post"))
                           )
                         )
                       )
              ),
              ## 在tabsetPanel中添加User Manual选项卡（Step4之后）----
              tabPanel("User Manual",
                       div(style = "padding: 20px; max-width: 1000px; margin: 0 auto;",
                           h2("User Manual", style = "color: #2c3e50; border-bottom: 2px solid #2c3e50; padding-bottom: 10px;"),
                           
                           h3("1. 工具概述", style = "color: #34495e;"),
                           p("本工具专为血浆蛋白质组学数据分析设计，提供以下核心功能："),
                           tags$ul(
                             tags$li("数据质量评估与可视化"),
                             tags$li("红细胞、血小板和凝血系统污染评估"),
                             tags$li("污染校正算法（基于稳健回归模型）"),
                             tags$li("差异表达分析与富集分析"),
                             tags$li("交互式结果可视化与数据导出")
                           ),
                           
                           h3("2. 使用指南", style = "color: #34495e;"),
                           h4("2.1 数据输入", style = "color: #7f8c8d;"),
                           tags$ol(
                             tags$li(strong("数据源选择："),"可选择示例数据进行参考，或上传CSV文件（基因表达矩阵和分组信息）"),
                             tags$li(strong("文件格式要求："),
                                     tags$ul(
                                       tags$li("表达矩阵：第一列为蛋白质名称，列为样本，需进行缺失值填充，无需进行log2转化（软件会自行进行log2转化）"),
                                       tags$li("分组信息：包含id（匹配表达矩阵列名）和group列")
                                     )),
                             tags$li(strong("参数设置："),
                                     "选择比较组，设置相关系数阈值（默认0.9）")
                           ),
                           
                           h4("2.2 污染评估", style = "color: #7f8c8d;"),
                           tags$ol(
                             tags$li(strong("质量评估："),"查看PCA、热图、相关系数分布等质量控制图表"),
                             tags$li(strong("标记物选择："),
                                     tags$ul(
                                       tags$li("从污染种类列表中选择Cv值较高的污染panel"),
                                       tags$li("通过相关性分析和差异表达筛选有效标记物")
                                     )),
                             tags$li(strong("污染水平："),
                                     tags$ul("通过CV分布评估各污染类型的影响程度，通过各样本的污染物marker表达来评估各样本的污染程度大小"),
                                     tags$ul("如果污染panel的CV值并不显著高于其他正常蛋白，或污染panel的marker无高相关形况，则该数据集无显著污染存在"),
                                     tags$ul("如果污染panel的marker在两分组均有显著差异，则无法分辨差异是污染还是生物学差异导致，无法进行矫正"))
                           ),
                           
                           h4("2.3 数据校正", style = "color: #7f8c8d;"),
                           tags$ol(
                             tags$li(strong("校正类型："),"选择需要校正的污染类型（红细胞、血小板、凝血系统），请勿选择无可用marker的污染类型进行矫正"),
                             tags$li(strong("约束因子："),"通过滑块调整校正强度（推荐范围0.8-1.2，默认为1）"),
                             tags$li(strong("质量控制："),"比较校正前后的PCA、污染markerCV变化等质量指标")
                           ),
                           
                           h4("2.4 差异分析", style = "color: #7f8c8d;"),
                           tags$ol(
                             tags$li(strong("分析方法："),"基于limma的差异表达分析"),
                             tags$li(strong("结果解读："),
                                     tags$ul(
                                       tags$li("通过Venn图比较校正前后差异蛋白重叠情况"),
                                       tags$li("火山图展示显著性差异蛋白")
                                     )),
                             tags$li(strong("数据导出："),"支持CSV格式的结果下载")
                           ),
                           
                           h3("3. 注意事项", style = "color: #34495e;"),
                           tags$ul(
                             tags$li(strong("数据预处理："),"建议进缺失值处理后再上传"),
                             tags$li(strong("标记物验证："),"需确保选择的污染标记物在数据集中稳定表达"),
                             tags$li(strong("参数优化："),"建议通过CV分布，相关系数分布图和PCA结果等调整约束因子，在绝大多数情况下，默认值即可满足需求"),
                             tags$li(strong("结果验证："),"校正后应观察到污染标记物的CV值显著降低，同时数据集的相关系数分布中高相关分布显著减少"),
                             tags$li(strong("技术支持："),"遇到问题请联系support@proteomics.com")
                           ),
                           
                           h3("4. 常见问题", style = "color: #34495e;"),
                           tags$ul(
                             tags$li(strong("Q1: "),"校正后数据出现负值怎么办？",
                                     "A: 这是正常现象，可能存在过小值，因为程序会自动进行log2转化"),
                             tags$li(strong("Q2: "),"如何确定最佳相关系数阈值？",
                                     "A: 默认0.9即可满足绝大多数情况下的需求，如果符合要求的marker过少，则可以适当放宽"),
                             tags$li(strong("Q3: "),"校正后差异蛋白数量显著变化是否正常？",
                                     "A: 是预期现象，说明污染对结果有显著影响，需结合生物学意义解读，正常情况下，去除的差异蛋白是由污染导致的，会富集到污染相关通路。而新增的差异蛋白则是被污染掩盖的符合生物学预期的蛋白，会富集到生物学相关通路")
                           )
                       )
              )
  )
)



# 定义 Server 逻辑 ----
server <- function(input, output, session) {
  ## 数据输入 ----
  ### 初始化 reactive values ----
  result_check <- reactiveVal()
  result_correct <- reactiveVal()
  result_de_pre <- reactiveVal()
  result_de_post <- reactiveVal()
  ### cor_cutoff响应值 ----
  # 新增cor_cutoff响应式值
  cor_cutoff <- reactive({
    input$cor_cutoff
  })
  # 同步滑块值
  observeEvent(input$cor_cutoff, {
    updateSliderInput(session, "cor_cutoff_step2", value = input$cor_cutoff)
  })
  observeEvent(input$cor_cutoff_step2, {
    updateSliderInput(session, "cor_cutoff", value = input$cor_cutoff_step2)
  })
  observeEvent(input$constraint_factor, {
    updateSliderInput(session, "constraint_factor_step2", value = input$constraint_factor)
  })
  observeEvent(input$constraint_factor_step2, {
    updateSliderInput(session, "constraint_factor", value = input$constraint_factor_step2)
  })
  constraint_factor <- reactive({
    # 当输入不存在时使用默认值1.0
    if (is.null(input$constraint_factor)) 1.0 else input$constraint_factor
  })
  ### 新增反应式值存储用户选择 ----
  selected_markers <- reactiveValues(
    erythrocyte = NULL,
    coagulation = NULL,
    platelet = NULL
  )
  ### 渲染marker表格 ----
  output$erythrocyte_marker_table <- renderDT({
    req(result_check())
    df <- result_check()$marker_stats$stats_erythrocyte
    datatable(df, selection = list(mode = 'multiple'), 
              options = list(pageLength = 5))
  })
  # 类似处理其他两个表格
  output$coagulation_marker_table <- renderDT({
    req(result_check())
    df <- result_check()$marker_stats$stats_coagulation
    datatable(df, selection = list(mode = 'multiple'), 
              options = list(pageLength = 5))
  })
  output$platelet_marker_table <- renderDT({
    req(result_check())
    df <- result_check()$marker_stats$stats_platelet
    datatable(df, selection = list(mode = 'multiple'), 
              options = list(pageLength = 5))
  })
  ### 获取用户选择 ----
  observe({
    req(result_check())
    
    # Erythrocyte
    erythrocyte_data <- result_check()$marker_stats$stats_erythrocyte 
    selected_erythrocyte <- erythrocyte_data$key[input$erythrocyte_marker_table_rows_selected]
    selected_markers$erythrocyte <- if(length(selected_erythrocyte) > 0) selected_erythrocyte else NULL
    
    # 类似处理其他两个类型
    # coagulation
    coagulation_data <- result_check()$marker_stats$stats_coagulation 
    selected_coagulation <- coagulation_data$key[input$coagulation_marker_table_rows_selected]
    selected_markers$coagulation <- if(length(selected_coagulation) > 0) selected_coagulation else NULL
    # platelet
    platelet_data <- result_check()$marker_stats$stats_platelet
    selected_platelet <- platelet_data$key[input$platelet_marker_table_rows_selected]
    selected_markers$platelet <- if(length(selected_platelet) > 0) selected_platelet else NULL
  })
  ## 加载示例数据 ----
  example_data <- reactive({
    df <- read.csv("./tests/raw_data_aggr.csv")  # 示例数据路径
    rownames(df) <- df[, 1]
    df[, -1, drop = FALSE]
  })
  
  example_group <- reactive({
    df <- read.csv("./tests/group.csv")  # 示例分组路径
    rownames(df) <- df[, 1]
    df
  })
  ## 读取数据 ----
  data <- reactive({
    if (input$data_source == "example") {
      example_data()
    } else {
      req(input$data_file)
      df <- read.csv(input$data_file$datapath)
      rownames(df) <- df[, 1]
      df[, -1, drop = FALSE]
    }
  })
  ### data_group
  data_group <- reactive({
    if (input$data_source == "example") {
      example_group()
    } else {
      req(input$group_file)
      df <- read.csv(input$group_file$datapath)
      rownames(df) <- df[, 1]
      df
    }
  })
  # 数据加载状态判断
  output$data_loaded <- reactive({
    !is.null(data()) && !is.null(data_group())
  })
  outputOptions(output, "data_loaded", suspendWhenHidden = FALSE)
  
  ## 展示原始数据 ----
  output$data_table <- renderDT({
    req(data())
    datatable(data(), options = list(pageLength = 10))
  })
  
  output$group_table <- renderDT({
    req(data_group())
    datatable(data_group(), options = list(pageLength = 10))
  })
  
  ## 更新分组选择 ----
  observe({
    req(data_group())
    if (!"group" %in% colnames(data_group())) {
      showNotification("Group file is missing 'group' column!", type = "error")
      return()
    }
    
    # 只有当用户选择去除生物学差异时才更新分组选择
    if (input$remove_biological_diff) {
      updateSelectInput(session, "group1", choices = unique(data_group()$group))
      updateSelectInput(session, "group2", choices = unique(data_group()$group))
    }
  })
  ## 更新DE analysis可使用条件 ----
  # Disable/enable DE button based on remove_biological_diff checkbox
  
  ## 数据检查 ----
  ### 数据检查（封装函数）----
  run_data_check <- function() {
    source("./R/data_check.R", local = TRUE)
    showModal(modalDialog("Running data check, please wait...", footer = NULL))
    # 根据用户选择决定是否传递分组信息
    if (input$remove_biological_diff && !is.null(input$group1)) {
      group1 <- input$group1
      group2 <- input$group2
      data_group <- data_group()
      DE_filter <- TRUE
    } else {
      group1 <- NULL
      group2 <- NULL
      data_group <- NULL
      DE_filter <- FALSE
    }
    check_result <- data_check(
      data = data(),
      data_group = data_group,
      DE_filter = DE_filter,
      cutoff = cor_cutoff(),
      group1 = input$group1,
      group2 = input$group2,
      custom_erythrocyte = selected_markers$erythrocyte,
      custom_coagulation = selected_markers$coagulation,
      custom_platelet = selected_markers$platelet
    )
    result_check(check_result)
    
  }
  ### Step1检查按钮 ----
  observeEvent(input$run_check, {
    run_data_check()
    updateTabsetPanel(session, "Step", selected = "Step 2: Check markers and contamination levels")
  })
  
  ### Step2重新检查按钮 ----
  observeEvent(input$rerun_check_step2, {
    run_data_check()
  })
  ## 数据校正 ----
  ### 处理校正类型选择逻辑 ----
  # 当选择"All"时自动选中其他三个
  observeEvent(input$type_all, {
    if (input$type_all) {
      updateCheckboxInput(session, "type_erythrocyte", value = TRUE)
      updateCheckboxInput(session, "type_platelet", value = TRUE)
      updateCheckboxInput(session, "type_coagulation", value = TRUE)
    }
  })
  
  # 当三个子项全选时自动选中"All"
  observe({
    all_selected <- all(
      input$type_erythrocyte,
      input$type_platelet,
      input$type_coagulation
    )
    updateCheckboxInput(session, "type_all", value = all_selected)
  })
  
  # 获取最终选择的校正类型
  selected_types <- reactive({
    types <- c()
    if (input$type_erythrocyte) types <- c(types, "erythrocyte")
    if (input$type_platelet) types <- c(types, "platelet")
    if (input$type_coagulation) types <- c(types, "coagulation")
    types
  })
  
  
  output$erythrocyte_checkbox <- renderUI({
    req(result_check())
    markers <- result_check()$marker_list$erythrocyte
    label <- if (length(markers) == 0) "Erythrocyte (no contamination)" else "Erythrocyte"
    checkboxInput("type_erythrocyte", label, 
                  value = if (length(markers) > 0) TRUE else FALSE)
  })
  
  output$platelet_checkbox <- renderUI({
    req(result_check())
    markers <- result_check()$marker_list$platelet
    label <- if (length(markers) == 0) "Platelet (no contamination)" else "Platelet"
    checkboxInput("type_platelet", label, 
                  value = if (length(markers) > 0) TRUE else FALSE)
  })
  
  output$coagulation_checkbox <- renderUI({
    req(result_check())
    markers <- result_check()$marker_list$coagulation
    label <- if (length(markers) == 0) "Coagulation (no contamination)" else "Coagulation"
    checkboxInput("type_coagulation", label, 
                  value = if (length(markers) > 0) TRUE else FALSE)
  })
  ### 主校正函数 ----
  run_correction <- function(constraint) {
    req(result_check())
    req(constraint)  # 确保约束因子存在
    
    source("./R/data_correct.R", local = TRUE)
    showModal(modalDialog("Performing data correction, please wait...", footer = NULL))
    
    # 根据选择的类型动态设置参数
    correction_type <- selected_types()
    
    # 确保至少选择一个类型
    if (length(correction_type) == 0) {
      showNotification("Please select at least one correction type!", type = "error")
      return(NULL)
    }
    
    # 调用矫正函数
    correct_result <- tryCatch({
      data_correct(
        data = result_check(), 
        type = correction_type,
        constraint = constraint,
        erythrocyte_marker = result_check()$marker_list$erythrocyte,
        coagulation_marker = result_check()$marker_list$coagulation,
        platelet_marker = result_check()$marker_list$platelet
      )
    }, error = function(e) {
      showNotification(paste("Correction failed:", e$message), type = "error")
      return(NULL)
    })
    
    removeModal()
    return(correct_result)
  }
  
  ### Step2校正按钮 ----
  observeEvent(input$run_correct, {
    correct_result <- run_correction(input$constraint_factor)
    if (!is.null(correct_result)) {
      result_correct(correct_result)
      updateTabsetPanel(session, "Step", selected = "Step 3: Correction Results")
    }
  })
  
  ### Step3重新校正按钮 ----
  observeEvent(input$run_correct_step2, {
    correct_result <- run_correction(input$constraint_factor_step2)
    if (!is.null(correct_result)) {
      result_correct(correct_result)
      showNotification("Correction re-run successfully!", type = "message")
    }
  })
  
  ## 差异表达分析 ----
  observeEvent(input$run_de, {
    # First check if remove_biological_diff is checked
    if (!isTRUE(input$remove_biological_diff)) {
      showModal(modalDialog(
        title = "Warning",
        "Please check 'Remove proteins with biological differences between groups' in Step 1 to perform differential expression analysis.",
        easyClose = TRUE,
        footer = NULL
      ))
      return()  # Exit the observer without running DE analysis
    }
    
    # Proceed with DE analysis only if remove_biological_diff is checked
    source("./R/de_analysis_module.R")
    req(result_correct())
    showModal(modalDialog("Running differential expression analysis, please wait...", footer = NULL))
    
    # 计算校正前结果
    de_pre <- limma_proteomics_analysis(
      expr_matrix = log2(result_correct()$rawdata),
      group_matrix = result_correct()$group,
      compare = c(input$group1, input$group2),
      p_type = "raw"
    )
    # 计算校正后结果
    de_post <- limma_proteomics_analysis(
      expr_matrix = log2(result_correct()$correct_data),
      group_matrix = result_correct()$group,
      compare = c(input$group1, input$group2),
      p_type = "raw"
    )
    
    result_de_pre(de_pre)
    result_de_post(de_post)
    removeModal()
    updateTabsetPanel(session, "Step", selected = "Step 4: Differential Expression")
  })
  ## 显示矫正后矩阵 ----
  output$data_correct_table <- renderDT({
    req(result_correct())
    datatable(result_correct()$correct_data, options = list(pageLength = 10))  # 每页显示 10 行
  })
  ## QC ----
  ### pre ----
  #### correlation ----
  output$correlation_p_pre_plot <- renderPlot({
    source("./R/plot_stat_distribution.R")
    req(result_check())
    plot_stat_distribution(data = result_check()$rawdata,
                           type = "pearson", statistic = "pvalue", alpha = 0.05)
  },height = 400,width = 500)
  output$correlation_r_pre_plot <- renderPlot({
    source("./R/plot_stat_distribution.R")
    req(result_check())
    plot_stat_distribution(data = result_check()$rawdata,
                           type = "pearson", statistic = "correlation", alpha = 0.05)
  },height = 400,width = 500)
  #### pca ----
  output$pca_pre_plot <- renderPlot({
    showModal(modalDialog("Running QC PCA, please wait...", footer = NULL))
    source("./R/modules/QC_PCA.R")
    req(result_check())
    removeModal()
    return(QC_PCA(data = result_check()$rawdata,
                  data_group = result_check()$group))
  },height = 400,width = 500)
  #### heatmap ----
  output$heatmap_pre_plot <- renderPlot({
    source("./R/modules/QC_heatmap.R")
    req(result_check())
    QC_heatmap(data = result_check()$rawdata,
               data_group = result_check()$group)
  })
  #### boxplot ----
  output$boxplot_pre_plot <- renderPlot({
    source("./R/modules/QC_boxplot.R")
    req(result_check())
    QC_boxplot(data = result_check()$rawdata,
               data_group = result_check()$group)
  })
  ### post ----
  #### correlation ----
  output$correlation_p_post_plot <- renderPlot({
    source("./R/plot_stat_distribution.R")
    req(result_check())  # 确保数据存在
    req(result_correct())
    plot_stat_distribution(data = result_correct()$correct_data,
                           type = "pearson", statistic = "pvalue", alpha = 0.05)
  },height = 400,width = 500)
  output$correlation_r_post_plot <- renderPlot({
    source("./R/plot_stat_distribution.R")
    req(result_check())  # 确保数据存在
    req(result_correct())
    plot_stat_distribution(data = result_correct()$correct_data,
                           type = "pearson", statistic = "correlation", alpha = 0.05)
  },height = 400,width = 500)
  #### pca ----
  output$pca_post_plot <- renderPlot({
    source("./R/modules/QC_PCA.R")
    req(result_check())  # 确保数据存在
    req(result_correct())
    QC_PCA(data = result_correct()$correct_data,
           data_group = result_correct()$group)
  },height = 400,width = 500)
  #### heatmap ----
  output$heatmap_post_plot <- renderPlot({
    source("./R/modules/QC_heatmap.R")
    req(result_check())  # 确保数据存在
    req(result_correct())
    QC_heatmap(data = result_correct()$correct_data,
               data_group = result_correct()$group)
  })
  #### boxplot ----
  output$boxplot_post_plot <- renderPlot({
    source("./R/modules/QC_boxplot.R")
    req(result_check())  # 确保数据存在
    req(result_correct())
    QC_boxplot(data = result_correct()$correct_data,
               data_group = result_correct()$group)
  })
  ## 相关性分析及可视化 ----
  ### ery ----
  output$cor_erythrocyte_plot <- renderPlot({
    source("./R/plot_expression_correlation.R")
    req(result_check())
    result <- plot_expression_correlation(exprMatrix = result_check()$correlation$erythrocyte$r,
                                          displayNumbers = T,input_type = "correlation")
    return(result$plot)
  },height = 400,width = 800)
  # output$cor_erythrocyte_data <- renderDT({
  #   req(result_check())
  #   datatable(result_check()$correlation$erythrocyte$r, options = list(pageLength = 10))  # 每页显示 10 行
  # })
  ### coa ----
  output$cor_coagulation_plot <- renderPlot({
    source("./R/plot_expression_correlation.R")
    req(result_check())
    result <- plot_expression_correlation(exprMatrix = result_check()$correlation$coagulation$r,
                                          displayNumbers = T,input_type = "correlation")
    return(result$plot)
  },height = 400,width = 800)
  # output$cor_coagulation_data <- renderDT({
  #   req(result_check())
  #   datatable(result_check()$correlation$coagulation$r, options = list(pageLength = 10))  # 每页显示 10 行
  # })
  ### platelet ----
  output$cor_platelet_plot <- renderPlot({
    source("./R/plot_expression_correlation.R")
    req(result_check())
    result <- plot_expression_correlation(exprMatrix = result_check()$correlation$platelet$r,
                                          displayNumbers = T,input_type = "correlation")
    return(result$plot)
  },height = 400,width = 800)
  # output$cor_platelet_data <- renderDT({
  #   req(result_check())
  #   datatable(result_check()$correlation$platelet$r, options = list(pageLength = 10))  # 每页显示 10 行
  # })
  ## 显示缺失基因 ----
  output$contamination_summary <- renderPrint({
    req(result_check())
    
    # 获取各个污染面板的统计数据和可用marker列表
    stats_ery <- result_check()$marker_stats$stats_erythrocyte
    stats_coa <- result_check()$marker_stats$stats_coagulation
    stats_plt <- result_check()$marker_stats$stats_platelet
    marker_list <- result_check()$marker_list
    
    # 计算统计信息的函数
    calculate_stats <- function(stats, panel) {
      total <- nrow(stats)
      missing_pct <- sum(stats$exists == "NA") / total * 100
      
      exists_count <- sum(stats$exists == "Pass")
      de_pct <- ifelse(exists_count > 0,
                       sum(stats$DE == "Non-removable inter-sample heterogeneity" & stats$exists == "Pass", na.rm = TRUE) / exists_count * 100,
                       0)
      cor_pct <- ifelse(exists_count > 0,
                        sum(stats$correlation == "Low statistical correlation" & stats$exists == "Pass", na.rm = TRUE) / exists_count * 100,
                        0)
      
      available <- marker_list[[panel]]
      available_text <- if (length(available) > 0) {
        paste("Available markers:", paste(available, collapse = ", "))
      } else {
        "No available markers; no significant contamination detected, correction not required."
      }
      
      sprintf("In %s contamination: %.1f%% missing, %.1f%% non-removable biological variation, %.1f%% low correlation. %s",
              panel, missing_pct, de_pct, cor_pct, available_text)
    }
    
    # 生成各面板报告
    ery_report <- calculate_stats(stats_ery, "erythrocyte")
    coa_report <- calculate_stats(stats_coa, "coagulation")
    plt_report <- calculate_stats(stats_plt, "platelet")
    
    # 组合输出
    cat(paste(ery_report, coa_report, plt_report, sep = "\n\n"))
  })
  
  ## 污染marker表达情况可视化 ----
  ### CV ----
  output$cv_pre_plot <- renderPlot({
    source("./R/modules/get_cv.R")
    req(result_check())
    data <- result_check()
    # 获取各个marker列表
    marker_coagulation <- data$marker_list$coagulation
    marker_erythrocyte <- data$marker_list$erythrocyte
    marker_platelet <- data$marker_list$platelet
    
    # 初始化存储各类型CV数据的列表
    cv_list <- list()
    # 处理每个marker类型，仅当存在时计算CV
    if (length(marker_coagulation) > 0) {
      result_cv_coa <- get_cv(raw_data = data$rawdata, protein = marker_coagulation)
      result_cv_coa$type <- "coagulation"
      cv_list <- c(cv_list, list(result_cv_coa))
    } else {
      showNotification("Coagulation markers are empty.", type = "warning")
    }
    
    if (length(marker_erythrocyte) > 0) {
      result_cv_ery <- get_cv(raw_data = data$rawdata, protein = marker_erythrocyte)
      result_cv_ery$type <- "erythrocyte"
      cv_list <- c(cv_list, list(result_cv_ery))
    } else {
      showNotification("Erythrocyte markers are empty.", type = "warning")
    }
    
    if (length(marker_platelet) > 0) {
      result_cv_pla <- get_cv(raw_data = data$rawdata, protein = marker_platelet)
      result_cv_pla$type <- "platelet"
      cv_list <- c(cv_list, list(result_cv_pla))
    } else {
      showNotification("Platelet markers are empty.", type = "warning")
    }
    # 处理other protein类型
    all_markers <- c(marker_coagulation, marker_erythrocyte, marker_platelet)
    if (length(all_markers) > 0) {
      other_proteins <- data$rawdata[!rownames(data$rawdata) %in% all_markers, ]
    } else {
      other_proteins <- data$rawdata
    }
    
    if (nrow(other_proteins) > 0) {
      result_cv_other <- get_cv(raw_data = other_proteins)
      result_cv_other$type <- "other protein"
      cv_list <- c(cv_list, list(result_cv_other))
    } else {
      showNotification("No other proteins available.", type = "warning")
    }
    
    # 合并所有数据，跳过空元素
    result_cv <- do.call(rbind, cv_list)
    
    # 如果没有数据可绘制，返回空白图
    if (is.null(result_cv) || nrow(result_cv) == 0) {
      return(ggplot() + 
               geom_blank() + 
               labs(title = "No data available for CV analysis.") +
               theme_classic())
    }
    
    # 确定存在的类型并设置因子水平
    existing_types <- unique(result_cv$type)
    type_levels <- c("coagulation", "erythrocyte", "platelet", "other protein")
    existing_levels <- intersect(type_levels, existing_types)
    result_cv$type <- factor(result_cv$type, levels = existing_levels)
    
    # 生成动态比较列表（仅包含存在的类型）
    comparisons <- list()
    if ("other protein" %in% existing_levels) {
      for (t in setdiff(existing_levels, "other protein")) {
        comparisons <- c(comparisons, list(c("other protein", t)))
      }
    }
    
    # 定义颜色映射（仅包含存在的类型）
    color_values <- c(
      "coagulation" = "#E64B35FF",
      "erythrocyte" = "#F39B7FFF",
      "platelet" = "#7E6148FF",
      "other protein" = "#3C5488FF"
    )[existing_levels]
    
    # 绘制主图
    p <- ggplot(result_cv, aes(x = type, y = CV, fill = type)) +
      geom_violin() +
      geom_boxplot(fill = "white", width = 0.15) +
      scale_fill_manual(values = color_values) +
      theme_minimal(base_size = 14) +
      theme(
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 13),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 13)
      )
    
    # 添加统计检验（当有比较对时）
    if (length(comparisons) > 0) {
      p <- p + stat_compare_means(
        comparisons = comparisons,
        method = "wilcox.test",
        label = "p.format",
        hide.ns = TRUE
      )
    }
    
    removeModal()
    return(p)
    
  }, height = 400, width = 500)
  output$cv_post_plot <- renderPlot({
    source("./R/modules/get_cv.R")
    req(result_check())  # 确保数据存在
    req(result_correct())
    data <- result_correct()
    # 获取各个marker列表
    marker_coagulation <- data$marker_list$coagulation
    marker_erythrocyte <- data$marker_list$erythrocyte
    marker_platelet <- data$marker_list$platelet
    
    # 初始化存储各类型CV数据的列表
    cv_list <- list()
    
    # 处理每个marker类型，仅当存在时计算CV
    if (length(marker_coagulation) > 0) {
      result_cv_coa <- get_cv(raw_data = data$correct_data, protein = marker_coagulation)
      result_cv_coa$type <- "coagulation"
      cv_list <- c(cv_list, list(result_cv_coa))
    } else {
      showNotification("Post-correction: Coagulation markers are empty.", type = "warning")
    }
    
    if (length(marker_erythrocyte) > 0) {
      result_cv_ery <- get_cv(raw_data = data$correct_data, protein = marker_erythrocyte)
      result_cv_ery$type <- "erythrocyte"
      cv_list <- c(cv_list, list(result_cv_ery))
    } else {
      showNotification("Post-correction: Erythrocyte markers are empty.", type = "warning")
    }
    
    if (length(marker_platelet) > 0) {
      result_cv_pla <- get_cv(raw_data = data$correct_data, protein = marker_platelet)
      result_cv_pla$type <- "platelet"
      cv_list <- c(cv_list, list(result_cv_pla))
    } else {
      showNotification("Post-correction: Platelet markers are empty.", type = "warning")
    }
    
    # 处理other protein类型（动态排除所有有效marker）
    all_markers <- c(marker_coagulation, marker_erythrocyte, marker_platelet)
    other_proteins <- data$correct_data[!rownames(data$correct_data) %in% all_markers, ]
    
    if (nrow(other_proteins) > 0) {
      result_cv_other <- get_cv(raw_data = other_proteins)
      result_cv_other$type <- "other protein"
      cv_list <- c(cv_list, list(result_cv_other))
    } else {
      showNotification("Post-correction: No other proteins available.", type = "warning")
    }
    
    # 合并所有数据（自动跳过空元素）
    result_cv <- do.call(rbind, cv_list)
    
    # 无数据时返回空白图
    if (is.null(result_cv) || nrow(result_cv) == 0) {
      return(ggplot() + 
               geom_blank() + 
               labs(title = "No data available for CV analysis (Post-correction)") +
               theme_classic())
    }
    
    # 动态设置因子水平
    existing_types <- unique(result_cv$type)
    type_levels <- c("coagulation", "erythrocyte", "platelet", "other protein")
    existing_levels <- intersect(type_levels, existing_types)
    result_cv$type <- factor(result_cv$type, levels = existing_levels)
    
    # 生成动态比较对
    comparisons <- list()
    if ("other protein" %in% existing_levels) {
      for (t in setdiff(existing_levels, "other protein")) {
        comparisons <- c(comparisons, list(c("other protein", t)))
      }
    }
    
    # 动态颜色映射
    color_values <- c(
      "coagulation" = "#E64B35FF",
      "erythrocyte" = "#F39B7FFF",
      "platelet" = "#7E6148FF",
      "other protein" = "#3C5488FF"
    )[existing_levels]
    
    # 绘制主图
    p <- ggplot(result_cv, aes(x = type, y = CV, fill = type)) +
      geom_violin(trim = FALSE) +
      geom_boxplot(fill = "white",width = 0.15) +
      scale_fill_manual(values = color_values) +
      labs(y = "Coefficient of Variation", x = "") +
      theme_minimal(base_size = 14) +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "right",
        panel.grid.major.x = element_blank()
      )
    
    # 智能添加统计标注
    if (length(comparisons) > 0) {
      p <- p + ggpubr::stat_compare_means(
        comparisons = comparisons,
        method = "wilcox.test",
        label = "p.format",
        tip.length = 0.01,
        step.increase = 0.1
      )
    }
    
    p
  }, height = 400, width = 500)
  ### erythrocyte ----
  
  output$contamination_erythrocyte_plot <- renderPlot({
    req(result_check())
    if (is.null(result_check()$plot_contamination)) {
      return(NULL)  # 如果 plot_contamination 为 NULL，不绘制图表
    }
    plot(result_check()$plot_contamination$erythrocyte)
  })
  ### platelet ----
  output$contamination_platelet_plot <- renderPlot({
    req(result_check())
    if (is.null(result_check()$plot_contamination)) {
      return(NULL)  # 如果 plot_contamination 为 NULL，不绘制图表
    }
    plot(result_check()$plot_contamination$platelet)
  })
  ### coagulation ----
  output$contamination_coagulation_plot <- renderPlot({
    req(result_check())
    if (is.null(result_check()$plot_contamination)) {
      return(NULL)  # 如果 plot_contamination 为 NULL，不绘制图表
    }
    plot(result_check()$plot_contamination$coagulation)
  })
  ## 样本污染情况可视化 ----
  ### pre ----
  #### erythrocyte marker ----
  output$erythrocyte_marker_pre_plot <- renderPlot({
    req(result_check())
    if (is.null(result_check()$plot_marker$erythrocyte)) {
      return(NULL)  # 如果 erythrocyte_marker_plot 为 NULL，不绘制图表
    }
    plot(result_check()$plot_marker$erythrocyte)
  })
  
  #### platelet marker ----
  output$platelet_marker_pre_plot <- renderPlot({
    req(result_check())
    if (is.null(result_check()$plot_marker$platelet)) {
      return(NULL)  # 如果 platelet_marker_plot 为 NULL，不绘制图表
    }
    plot(result_check()$plot_marker$platelet)
  })
  
  #### coagulation marker ----
  output$coagulation_marker_pre_plot <- renderPlot({
    req(result_check())
    if (is.null(result_check()$plot_marker$coagulation)) {
      return(NULL)  # 如果 coagulation_marker_plot 为 NULL，不绘制图表
    }
    print(result_check()$plot_marker$coagulation)
  })
  ### post ----
  #### erythrocyte marker ----
  output$erythrocyte_marker_post_plot <- renderPlot({
    req(result_check())  # 确保数据存在
    req(result_correct())
    marker <- result_check()$marker_list$erythrocyte
    source("./R/plot_protein_by_sample.R")
    if (length(marker) > 0) {
      plot_protein_by_sample(data = result_correct()$correct_data[rownames(result_correct()$correct_data)%in%result_correct()$marker_list$erythrocyte,])
    }
    
  })
  
  #### platelet marker ----
  output$platelet_marker_post_plot <- renderPlot({
    req(result_check())  # 确保数据存在
    req(result_correct())
    marker <- result_check()$marker_list$platelet
    source("./R/plot_protein_by_sample.R")
    if (length(marker) > 0) {
      plot_protein_by_sample(data = result_correct()$correct_data[rownames(result_correct()$correct_data)%in%result_correct()$marker_list$platelet,])
    }
    
  })
  
  #### coagulation marker ----
  output$coagulation_marker_post_plot <- renderPlot({
    req(result_check())  # 确保数据存在
    req(result_correct())
    marker <- result_check()$marker_list$coagulation
    source("./R/plot_protein_by_sample.R")
    if (length(marker) > 0) {
      plot_protein_by_sample(data = result_correct()$correct_data[rownames(result_correct()$correct_data)%in%result_correct()$marker_list$coagulation,])
    }
    
  })
  
  # corrected_plot 显示校正后污染水平可视化 ----
  output$corrected_plot <- renderPlot({
    req(result_check())  # 确保数据存在
    req(result_correct())
    if (is.null(result_correct()$contamination_level)) {
      return(NULL)  # 如果 contamination_level 为 NULL，不绘制图表
    }
    plot_contamination(result_correct()$contamination_level, "校正后污染水平可视化", "corrected_plot")
  })
  
  
  # 显示差异表达结果 ----
  output$result_de_pre_table <- renderDT({
    req(result_check())  # 确保数据存在
    req(result_de_pre())
    datatable(result_de_pre(), options = list(pageLength = 10))  # 每页显示 10 行
  })
  output$result_de_post_table <- renderDT({
    req(result_check())  # 确保数据存在
    req(result_de_post())
    datatable(result_de_post(), options = list(pageLength = 10))  # 每页显示 10 行
  })
  # 绘制VN图 ----
  output$vn_plot <- renderPlot({
    req(result_check())  # 确保数据存在
    req(result_de_post())
    req(result_de_pre())
    venn_plot(
      set1 = result_de_post()[result_de_post()$significant == TRUE,"Protein"],
      set2 = result_de_pre()[result_de_pre()$significant == TRUE,"Protein"],
      categories = c("post", "pre"),
      title = "Two differences analysed results",
      colors = c("#ae6b81", "#6982b9"),
      alpha = 0.6,
      print.mode = "raw"
    )
  },height = 500,width = 500)
  # 绘制火山图 ----
  output$volc_de_pre <- renderPlot({
    req(result_check())  # 确保数据存在
    req(result_de_pre())
    plot_volc <- create_volcano_plot(result_de_pre(), 
                                     p_type = "raw", 
                                     p_cutoff = 0.05, 
                                     logFC_cutoff = 0,
                                     gene_col = "Protein",
                                     group_names = c(input$group1,input$group2),
                                     colors = c(Up = "#E64B35", Down = "#4DBBD5", Not = "grey80"))
    print(plot_volc)
  },height = 500,width = 600)
  output$volc_de_post <- renderPlot({
    req(result_check())  # 确保数据存在
    req(result_de_post())
    plot_volc <- create_volcano_plot(result_de_post(), 
                                     p_type = "raw", 
                                     p_cutoff = 0.05, 
                                     logFC_cutoff = 0,
                                     gene_col = "Protein",
                                     group_names = c(input$group1,input$group2),
                                     colors = c(Up = "#E64B35", Down = "#4DBBD5", Not = "grey80"))
    print(plot_volc)
  },height = 500,width = 600)
  # 下载校正数据 ----
  output$download_data <- downloadHandler(
    filename = function() {
      "corrected_data.csv"
    },
    content = function(file) {
      write.csv(result_correct()$correct_data, file)
    }
  )
  # 下载相关性结果 ----
  output$download_cor_data_erythrocyte <- downloadHandler(
    filename = function() {
      "correlation_matrix_erythrocyte.csv"
    },
    content = function(file) {
      req(result_check())
      write.csv(result_check()$correlation$erythrocyte$r, file)
    }
  )
  output$download_cor_data_coagulation <- downloadHandler(
    filename = function() {
      "correlation_matrix_coagulation.csv"
    },
    content = function(file) {
      req(result_check())
      write.csv(result_check()$correlation$coagulation$r, file)
    }
  )
  output$download_cor_data_platelet <- downloadHandler(
    filename = function() {
      "correlation_matrix_platelet.csv"
    },
    content = function(file) {
      req(result_check())
      write.csv(result_check()$correlation$platelet$r, file)
    }
  )
  # 下载差异表达结果 ----
  ## 矫正前 ----
  output$download_de_pre <- downloadHandler(
    filename = function() {
      "de_results.csv"
    },
    content = function(file) {
      write.csv(result_de_pre(), file)
    }
  )
  ## 矫正后 ----
  output$download_de_post <- downloadHandler(
    filename = function() {
      "de_results.csv"
    },
    content = function(file) {
      write.csv(result_de_post(), file)
    }
  )
}
# 运行 Shiny 应用 ----
shinyApp(ui = ui, server = server)
