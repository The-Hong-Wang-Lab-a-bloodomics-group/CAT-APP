# input packages ---- 
library(shiny)
library(ggplot2)
library(ggsignif)
library(ggpubr)
library(tidyr)
library(dplyr)
library(readxl)
library(DT)
library(RSQLite)
options(shiny.maxRequestSize=30*1024^2)

# 定义 UI ----
ui <- fluidPage(
  tags$head(
    #确保 Umami 脚本被加载
    tags$script(defer = NA, src = "//umami.bloodecosystem.com/script.js",
                `data-website-id` = "a8952081-8633-4fd2-a59f-4b32d27a8227"),
  tags$style(HTML("
      /* 全局字体调整 */
      body {
        font-size: 16px;  /* 基础字体增大 */
        line-height: 1.6; /* 增加行高提升可读性 */
      }
      
      /* 标题层级优化 */
      h1, h2, h3, h4, h5 {
        font-weight: 600; /* 标题加粗 */
        margin-top: 20px;
        margin-bottom: 15px;
      }
      
      h1 { font-size: 2.0rem; } /* 约32px */
      h2 { font-size: 1.8rem; } /* 约28px */
      h3 { font-size: 1.6rem; } /* 约25px */
      h4 { font-size: 1.4rem; } /* 约22px */
      h5 { font-size: 1.2rem; } /* 约19px */
      
      /* 面板标题优化 */
      .panel-title { 
        font-size: 1.5rem; 
        font-weight: bold;
      }
      
      /* 表格文字优化 */
      .dataTables_wrapper {
        font-size: 16px;
      }
      
      /* 按钮文字优化 */
      .btn {
        font-size: 16px;
        padding: 8px 16px;
      }
      
      /* 滑块文字优化 */
      .irs-min, .irs-max, .irs-single {
        font-size: 14px;
      }
      
      /* 选项卡优化 */
      .nav-tabs > li > a {
        font-size: 16px;
        padding: 10px 15px;
      }
      
      /* 表单控件文字 */
      .control-label {
        font-size: 16px;
        font-weight: 500;
      }
      
      /* 帮助图标放大 */
      .glyphicon-question-sign {
        font-size: 18px;
        margin-left: 5px;
        color: #337ab7;
      }
      
      /* 增加内容区域的内边距 */
      .main-content {
        padding: 20px 25px;
      }
      
      /* 欢迎页特定样式 */
      #welcome-content {
        font-size: 18px;
        line-height: 1.7;
      }
      
      #welcome-content li {
        margin-bottom: 10px;
        font-size: 17px;
      }
            /* 数据库计数样式 */
      #db-counter {
        text-align: center;
        font-weight: bold;
        margin-top: 10px;
        margin-bottom: 20px;
        color: #2c3e50;
        background-color: #f8f9fa;
        padding: 8px;
        border-radius: 4px;
        border: 1px solid #eaeaea;
      }
    "))
),
  titlePanel("Contamination Analysis and Tempering-An Automated Online Platform for Plasma Proteomics"),
  tabsetPanel(id = "Step",
              ## Welcome Tab ----
              tabPanel("Welcome",
                       # 外层容器：确保内容整体居中
                       div(style = "display: flex; justify-content: center; width: 100%;",
                           # 内容区域：固定最大宽度并左对齐
                           div(style = "max-width: 900px; width: 100%; text-align: left; padding: 20px;",
                               
                               # 标题（居中显示）
                               h2("Welcome to CAT-APP", 
                                  style = "text-align: center; font-size: calc(20px + 1vw); margin-bottom: 25px;"),
                               
                               # 图片（居中显示但内容左对齐）
                               div(style = "text-align: center; margin: 25px 0;",
                                   img(src = "Welcome_2.png", 
                                       style = "max-width: 100%; height: auto; border-radius: 0px;")
                               ),
                               
                               # 正文（强制左对齐）
                               div(style = "text-align: left;",  # 显式覆盖可能继承的居中样式
                                   p("This interactive tool allows you to analyze and correct for contamination in plasma proteomics data."),
                                   
                                   h3("Key features:", style = "margin-top: 25px;"),
                                   tags$ul(
                                     style = "padding-left: 20px;",
                                     tags$li("Multi-dimensional contamination assessment and adaptive contamination indexing"),
                                     tags$li("Mathematic model-based contamination correction"),
                                     tags$li("Data recovery evaluation with visualization")
                                   ),
                                   
                                   h3("How to use:", style = "margin-top: 25px;"),
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
                                 "CAT-APP is developed by R shiny (Version 1.11.1), and is free and open to all users with no login requirement. 
                                 For any questions or feedback, please contact us via email: ",
                                 tags$a(
                                   href = "mailto:zhangdong_0121@foxmail.com",
                                   target = "_blank",
                                   "zhangdong_0121@foxmail.com"
                                 )
                               ),
                               p("© 2025 CAT-APP - Contamination Analysis Tool for Plasma Proteomics")
                           )
                       )
              ),
              ## Step1 ----
              # 改为 data input
              tabPanel("Step 1: Data Input",
                       div(class = "main-content",
                       sidebarLayout(
                         sidebarPanel(h3("Step 1: Data Input",
                                         tags$span(
                                           id = 'span_step1_data_input',
                                           `data-toggle` = "tooltip",
                                           title = 'In this part, users can upload their own proteomics expression data and sample group data.The example data can be found when users click "Load example data" below.Detailed descriptions are provided in the "Help" part.',
                                           tags$span(class = "glyphicon glyphicon-question-sign")
                                         )),
                                      radioButtons("data_source", "Select Data Source",
                                                   choices = list("Load experimental data" = "experimental", 
                                                                  "Load example data" = "example"),
                                                   selected = "experimental"),
                                      downloadButton("download_example_data", "Download example data"),
                                      downloadButton("download_example_group_info", "Download example group info file"),
                                      conditionalPanel(
                                        condition = "input.data_source == 'experimental'",
                                        fileInput("data_file", "Upload Data File (CSV)", 
                                                  accept = ".csv"),
                                        fileInput("group_file", "Upload Group Info File (CSV)", 
                                                  accept = ".csv")
                                      ),
                                      h3("Grouping Settings",
                                         tags$span(
                                           id = 'span_step1_grouping_settings',
                                           `data-toggle` = "tooltip",
                                           title = 'In this step, users can set biological groups in the data to exclude markers that show significant differences between biological groups, preventing over-correction for variations arising from biological differences. Additionally, setting groups is fundamental for performing subsequent differential analysis. Here, group1 represents the experimental group, while group2 represents the control group.',
                                           tags$span(class = "glyphicon glyphicon-question-sign")
                                         )),
                                      # 新增去除生物学差异选项
                                      checkboxInput("remove_biological_diff", 
                                                    "Filter contamination marker proteins with inter-group bias", 
                                                    value = TRUE),
                                      # 条件面板控制分组设置显示
                                      conditionalPanel(
                                        condition = "input.remove_biological_diff == true",
                                        selectInput("group1", "Group 1", choices = NULL),
                                        selectInput("group2", "Group 2", choices = NULL)
                                      ),
                                      h3("Correlation Cutoff", style = "margin: 0; margin-right: 5px;",
                                      tags$span(
                                        id = 'span_step1_correlation_cutoff',
                                        `data-toggle` = "tooltip",
                                        title = 'This value is the cutoff for screening highly correlated markers. Any marker with a maximum correlation below this value will be removed. When markers are present, a higher value leads to stricter results.',
                                        tags$span(class = "glyphicon glyphicon-question-sign")
                                      )),
                                      sliderInput("cor_cutoff", 
                                                  # label = "Correlation Cutoff",
                                                  label = NULL,
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
                         ))))
              ),
              ## Step2 ----
              # 拆分为data evaluation和define markers
              tabPanel("Step 2: Check markers and contamination levels",
                       div(class = "main-content",
                       sidebarLayout(
                         sidebarPanel(
                           sliderInput("cor_cutoff_step2", "Correlation Cutoff",
                                       min = 0.5, max = 0.99, value = 0.9, step = 0.01),
                           actionButton("rerun_check_step2", "Re-run Data Check"),
                           h3("Step 2: Contamination Assessment"),
                           h3("Correction Type", style = "margin: 0; margin-right: 5px;",
                              tags$span(
                                id = 'span_step2_contamination_assessment',
                                `data-toggle` = "tooltip",
                                title = 'CAT-APP automatically screens markers based on the previous settings. If no markers meet the requirements, it indicates no significant contamination in this panel.CAT-APP will automatically select panels that require correction.',
                                tags$span(class = "glyphicon glyphicon-question-sign")
                              )),
                           # 修改后的单选按钮组
                           tags$div(class = "form-group",
                                    tags$label(class = "control-label", label = NULL),
                                    checkboxInput("type_all", "All", value = FALSE),
                                    uiOutput("erythrocyte_checkbox"),
                                    uiOutput("platelet_checkbox"),
                                    uiOutput("coagulation_checkbox")),
                           h3("constraint factor", style = "margin: 0; margin-right: 5px;",
                              tags$span(
                                id = 'span_step2_constraint_factor',
                                `data-toggle` = "tooltip",
                                title = 'If users consider the correction effect too strong or too weak, they can adjust the correction strength by setting the constraint factor.',
                                tags$span(class = "glyphicon glyphicon-question-sign")
                              )),
                           sliderInput("constraint_factor", label = NULL,
                                       # label = "constraint factor",
                                       min = 0.5, max = 1.5, value = 1.0, step = 0.01),
                           actionButton("run_correct", "Run Correction")),
                         mainPanel(h3("Data Quality Assessment"),
                                   h4("Contamination Summary",
                                      tags$span(
                                        id = 'span_step2_contamination_summary',
                                        `data-toggle` = "tooltip",
                                        title = "In this step, CAT-APP summarizes the status of markers in each panel:
- 'missing' indicates that the marker was not detected in the dataset;
- 'non-removable biological variation' indicates that the marker shows significant biological differences between the two groups and will be removed in subsequent analysis;
- 'low correlation' indicates that the marker does not have high correlation characteristics and will also be removed in subsequent analysis.",
                                        tags$span(class = "glyphicon glyphicon-question-sign")
                                      )),
                                   verbatimTextOutput("contamination_summary"),
                                   tabsetPanel(
                                     tabPanel("Erythrocyte", DTOutput("erythrocyte_marker_table")),
                                     tabPanel("Coagulation", DTOutput("coagulation_marker_table")),
                                     tabPanel("Platelet", DTOutput("platelet_marker_table"))
                                   ),
                                   h4("QC(Pre-correction)",
                                      tags$span(
                                        id = 'span_step2_pre-correction_qc',
                                        `data-toggle` = "tooltip",
                                        title = "CAT-APP performs basic quality control on the data, including:
- Correlation: Displays the correlation distribution between each protein and sample contamination levels. Normally, it should be a normal distribution centered around 0. If there is a peak near 1, it indicates the presence of highly contaminated proteins.
- PCA: Shows the distribution of samples on PCA plots, allowing users to check for potential contaminated outlier samples and verify whether samples cluster as expected.
- Heatmap: Reveals whether there are abnormally high-expressed proteins, helping to identify potential contaminated protein expression patterns.
- Boxplot: Displays the protein intensity distribution across samples, enabling users to detect potential contaminated samples.",
                                        tags$span(class = "glyphicon glyphicon-question-sign")
                                      )),
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
                                   h4("Contamination Marker Expression",
                                      tags$span(
                                        id = 'span_step2_contamination_marker_expression',
                                        `data-toggle` = "tooltip",
                                        title = "CAT-APP automatically calculates marker expression for each sample. 
If grouping information was provided earlier, it will filter out markers with significant differences between groups, which will not be included in subsequent analyses.",
                                        tags$span(class = "glyphicon glyphicon-question-sign")
                                      )),
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
                                   h4("Relevance of contamination markers",
                                      tags$span(
                                        id = 'span_step2_relevance_of_contamination_markers',
                                        `data-toggle` = "tooltip",
                                        title = "CAT-APP automatically calculates correlations between markers and screens for highly correlated markers based on the Correlation Cutoff setting. 
If contamination exists, the corresponding markers should exhibit highly correlated characteristics. Users can download the correlation results for further analysis.",
                                        tags$span(class = "glyphicon glyphicon-question-sign")
                                      )),
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
                                   h4("Contamination Levels(Pre-correction)",
                                      tags$span(
                                        id = 'span_step2_contamination_levels',
                                        `data-toggle` = "tooltip",
                                        title = "CAT-APP automatically calculates the CV value distribution of markers and other proteins in each panel:
- If there is a significant difference between the CV values of a panel and other normal proteins, it indicates significant contamination in that panel for this dataset;
- If there is no significant difference in CV values compared to other proteins, but highly correlated markers exist, it suggests potential mild contamination in that panel for this dataset;
- If there is no significant difference in CV values compared to other proteins, and no highly correlated markers exist, it indicates no contamination in that panel for this dataset, and no correction is needed.
Additionally, CAT-APP provides visualization of contamination indices for each sample across different panels, allowing for intuitive identification of samples with potential contamination issues.",
                                        tags$span(class = "glyphicon glyphicon-question-sign")
                                      )),
                                   tabsetPanel(
                                     tabPanel("CV Analysis", plotOutput("cv_pre_plot")),
                                     tabPanel("Erythrocyte", plotOutput("erythrocyte_marker_pre_plot")),
                                     tabPanel("Coagulation", plotOutput("coagulation_marker_pre_plot")),
                                     tabPanel("Platelet", plotOutput("platelet_marker_pre_plot"))
                                   )
                         )
                       ))
              ),
              ## Step3 ----
              # correction
              tabPanel("Step 3: Correction Results",
                       div(class = "main-content",
                       sidebarLayout(
                         sidebarPanel(h3("Step 3: Correction Analysis"),
                                      sliderInput("constraint_factor_step2", "constraint factor",
                                                  min = 0.5, max = 1.5, value = 1.0, step = 0.01),
                                      actionButton("run_correct_step2", "Re-run Correction"),
                                      h3("Step 4: DE Analysis",
                                         tags$span(
                                           id = 'span_step3_DE_Analysis',
                                           `data-toggle` = "tooltip",
                                           title = 'Running this differential analysis requires grouping information from Step 1. This analysis is based on limma and will automatically perform log2 transformation on the data.',
                                           tags$span(class = "glyphicon glyphicon-question-sign")
                                         )),
                                      actionButton("run_de", "Run Differential Expression Analysis")
                         ),
                         mainPanel(h3("Correction Outcomes"), 
                                   h4("QC(Post-correction)",
                                      tags$span(
                                        id = 'span_step3_post-correction_qc',
                                        `data-toggle` = "tooltip",
                                        title = "CAT-APP performs basic quality control on the data, including:
- Correlation: Displays the correlation distribution between each protein and sample contamination levels. Normally, it should be a normal distribution centered around 0. If there is a peak near 1, it indicates the presence of highly contaminated proteins.
- PCA: Shows the distribution of samples on PCA plots, allowing users to check for potential contaminated outlier samples and verify whether samples cluster as expected.
- Heatmap: Reveals whether there are abnormally high-expressed proteins, helping to identify potential contaminated protein expression patterns.
- Boxplot: Displays the protein intensity distribution across samples, enabling users to detect potential contaminated samples.",
                                        tags$span(class = "glyphicon glyphicon-question-sign")
                                      )),
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
                                   h4("Contamination Levels(Post-correction)",
                                      tags$span(
                                        id = 'span_step2_contamination_levels',
                                        `data-toggle` = "tooltip",
                                        title = "CAT-APP automatically calculates the CV value distribution of markers and other proteins in each panel:
- If there is a significant difference between the CV values of a panel and other normal proteins, it indicates significant contamination in that panel for this dataset;
- If there is no significant difference in CV values compared to other proteins, but highly correlated markers exist, it suggests potential mild contamination in that panel for this dataset;
- If there is no significant difference in CV values compared to other proteins, and no highly correlated markers exist, it indicates no contamination in that panel for this dataset, and no correction is needed.
Additionally, CAT-APP provides visualization of contamination indices for each sample across different panels, allowing for intuitive identification of samples with potential contamination issues.",
                                        tags$span(class = "glyphicon glyphicon-question-sign")
                                      )),
                                   tabsetPanel(
                                     tabPanel("CV Analysis", plotOutput("cv_post_plot")),
                                     tabPanel("Erythrocyte", plotOutput("erythrocyte_marker_post_plot")),
                                     tabPanel("Coagulation", plotOutput("coagulation_marker_post_plot")),
                                     tabPanel("Platelet", plotOutput("platelet_marker_post_plot"))
                                   ),
                                   h4("Post-correction Data Matrix"),
                                   downloadButton("download_data", "Download Post-correction Data"),
                                   DTOutput("data_correct_table")
                         )
                       ))),
              ## Step4 DE ----
              # DE & enrichment
              tabPanel("Step 4: Differential Expression",
                       div(class = "main-content",
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
                           tabPanel("Pre-correction Results",
                                    downloadButton("download_de_pre", "Download Pre-correction DE Results"),
                                    DTOutput("result_de_pre_table")),
                           tabPanel("Post-correction Results",
                                    downloadButton("download_de_post", "Download Post-correction DE Results"),
                                    DTOutput("result_de_post_table")
                           )),
                           tabsetPanel(
                             tabPanel("Venn Diagram", plotOutput("vn_plot")),
                             tabPanel("Volcano Plot (Pre-correction)", plotOutput("volc_de_pre")),
                             tabPanel("Volcano Plot (Post-correction)", plotOutput("volc_de_post"))
                           )
                         )
                       ))
              ),
              # User Manual ----
              tabPanel("User Manual",
                       div(class = "main-content",
                       div(style = "padding: 20px; max-width: 1000px; margin: 0 auto;",
                           h2("User Manual", style = "color: #2c3e50; border-bottom: 2px solid #2c3e50; padding-bottom: 10px;"),
                           
                           h3("1. Tool Overview", style = "color: #34495e;"),
                           p("This tool is designed for plasma proteomics data analysis and provides the following core functions:"),
                           tags$ul(
                             style = "padding-left: 20px;",
                             tags$li("Multi-dimensional contamination assessment and adaptive contamination indexing"),
                             tags$li("Mathematic model-based contamination correction"),
                             tags$li("Data recovery evaluation with visualization")
                           ),
                           
                           h3("2. User Guide", style = "color: #34495e;"),
                           h4("2.1 Data Input", style = "color: #7f8c8d;"),
                           tags$ol(
                             tags$li(strong("Data source selection:"), "Choose example data for reference or upload CSV files (gene expression matrix and group information)"),
                             tags$li(strong("File format requirements:"),
                                     tags$ul(
                                       tags$li("Expression matrix: First column contains protein names, columns represent samples. Requires missing value imputation. Do NOT perform log2 transformation (software will automatically apply log2 transformation)"),
                                       tags$li("Group information: Must contain id (matching expression matrix column names) and group columns")
                                     )),
                             tags$li(strong("Parameter settings:"), "Select comparison groups, set correlation coefficient threshold (default: 0.9)")
                           ),
                           
                           h4("2.2 Contamination Assessment", style = "color: #7f8c8d;"),
                           tags$ol(
                             tags$li(strong("Quality assessment:"), "View quality control plots including PCA, heatmap, correlation coefficient distribution"),
                             tags$li(strong("Marker selection:"),
                                     tags$ul(
                                       tags$li("Select contamination panels with high CV values from contamination type list"),
                                       tags$li("Filter effective markers through correlation analysis and differential expression")
                                     )),
                             tags$li(strong("Contamination level:"),
                                     tags$ul("Assess impact degree through CV distribution"),
                                     tags$ul("Evaluate sample-specific contamination through expression of contaminant markers"),
                                     tags$ul("If the CV value of a contamination panel is not significantly higher than other proteins, or if markers show no high correlation, the dataset has no significant contamination"),
                                     tags$ul("If contaminant markers show significant differential expression in both groups, correction cannot be performed as differences may originate from either contamination or biological variation"))
                           ),
                           
                           h4("2.3 Data Correction", style = "color: #7f8c8d;"),
                           tags$ol(
                             tags$li(strong("Correction type:"), "Select contamination types to correct (RBC, platelets, coagulation system).", strong("Do not select types without available markers")),
                             tags$li(strong("Constraint factor:"), "Adjust correction strength using slider (recommended range: 0.8-1.2, default: 1)"),
                             tags$li(strong("Quality control:"), "Compare quality metrics pre/post correction: PCA, contaminant marker CV changes")
                           ),
                           
                           h4("2.4 Differential Analysis", style = "color: #7f8c8d;"),
                           tags$ol(
                             tags$li(strong("Analysis method:"), "Differential expression analysis based on limma"),
                             tags$li(strong("Result interpretation:"),
                                     tags$ul(
                                       tags$li("Compare overlapping differential proteins pre/post correction using Venn diagrams"),
                                       tags$li("Visualize significant differential proteins via volcano plots")
                                     )),
                             tags$li(strong("Data export:"), "Download results in CSV format")
                           ),
                           
                           h3("3. Important Notes", style = "color: #34495e;"),
                           tags$ul(
                             tags$li(strong("Data preprocessing:"), "Perform missing value imputation before uploading"),
                             tags$li(strong("Marker validation:"), "Ensure selected contamination markers show stable expression in the dataset."),
                             tags$li(strong("Parameter optimization:"), "Adjust constraint factor using CV distribution, correlation plots and PCA results. Default values suffice for most cases"),
                             tags$li(strong("Result validation:"), "Post-correction should show: Significant reduction in CV values of contaminant markers and decreased high-correlation distribution"),
                             tags$li(strong("Technical support:"), "please contact us via email: zhangdong_0121@foxmail.com")
                           ),
                           
                           h3("4. Frequently Asked Questions", style = "color: #34495e;"),
                           tags$ul(
                             tags$li(strong("Q1: "), "Why do negative values appear after correction?",
                                     "A: This is normal and may occur with extremely small values due to automatic log2 transformation"),
                             tags$li(strong("Q2: "), "How to determine optimal correlation coefficient threshold?",
                                     "A: Default 0.9 works for most cases. Lower threshold if insufficient markers are identified"),
                             tags$li(strong("Q3: "), "What are Contamination Levels? How are they calculated?",
                                     "A: Contamination Levels are values calculated by CAT-APP for each sample, derived from the average expression of markers that are highly correlated within the dataset and show no significant differences between biological groups."),
                             tags$li(strong("Q4: "), "Is significant change in differential proteins post-correction normal?",
                                     "A: Yes. Removed proteins typically associate with contamination pathways, while new differential proteins often relate to biological pathways")
                           )
                       )
              ))
              
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
  run_count <- reactiveVal(get_run_count())
  #### 数据检查 ----
  na_check_status <- reactiveVal(NULL)
  negative_check_status <- reactiveVal(NULL)
  ##### 检查数据缺失值的函数 ----
  check_na_values <- function(data) {
    if (is.null(data)) return(NULL)
    
    total_values <- dim(data)[1] * dim(data)[2]
    na_count <- sum(is.na(data))
    na_percentage <- round((na_count / total_values) * 100, 2)
    
    return(list(
      has_na = na_count > 0,
      na_count = na_count,
      na_percentage = na_percentage,
      total_values = total_values
    ))
  }
  ##### 检查数据负值的函数 ----
  check_negative_values <- function(data) {
    if (is.null(data)) return(NULL)
    
    total_values <- dim(data)[1] * dim(data)[2]
    negative_count <- sum(data < 0, na.rm = TRUE)
    negative_percentage <- round((negative_count / total_values) * 100, 2)
    
    return(list(
      has_negative = negative_count > 0,
      negative_count = negative_count,
      negative_percentage = negative_percentage,
      total_values = total_values
    ))
  }
  ##### Function to display missing value warning ----
  show_na_warning <- function(na_info, data_type = "expression data") {
    if (!is.null(na_info) && na_info$has_na) {
      showModal(modalDialog(
        title = strong("Missing Data Warning", style = "color: #d9534f;"),
        tagList(
          p(icon("exclamation-triangle"), 
            sprintf("Missing values detected in your %s!", data_type)),
          hr(),
          p(sprintf("Total missing values: %d / %d (%.2f%%)", 
                    na_info$na_count, na_info$total_values, na_info$na_percentage)),
          hr(),
          p(strong("Recommended actions:")),
          tags$ul(
            tags$li("Please perform missing value imputation before uploading data"),
            tags$li("Common imputation methods: min value, median, KNN imputation, etc."),
            tags$li("Re-upload data after processing for analysis")
          )
        ),
        footer = tagList(
          # modalButton("Ignore and Continue"),
          actionButton("cancel_analysis", "Cancel Analysis", class = "btn-danger")
        ),
        size = "l"
      ))
    }
  }
  ##### Function to display negative value warning ----
  show_negative_warning <- function(negative_info, data_type = "expression data") {
    if (!is.null(negative_info) && negative_info$has_negative) {
      showModal(modalDialog(
        title = strong("Negative Values Warning", style = "color: #f0ad4e;"),
        tagList(
          p(icon("exclamation-triangle"), 
            sprintf("Negative values detected in your %s!", data_type)),
          hr(),
          p(sprintf("Total negative values: %d / %d (%.2f%%)", 
                    negative_info$negative_count, negative_info$total_values, 
                    negative_info$negative_percentage)),
          hr(),
          p(strong("Recommended actions:")),
          tags$ul(
            tags$li("Please check if your data has been properly processed"),
            tags$li("Ensure data has not been log-transformed before upload"),
            tags$li("Consider data normalization or transformation if needed"),
            tags$li("Verify data quality and processing pipeline")
          )
        ),
        footer = tagList(
          actionButton("cancel_analysis", "Cancel Analysis", class = "btn-warning")
        ),
        size = "l"
      ))
    }
  }
  
  ### 处理取消分析按钮 ----
  observeEvent(input$cancel_analysis, {
    removeModal()
    # updateTabsetPanel(session, "Step", selected = "Step 1: Data Input")
    session$reload()
  })
  ### cor_cutoff响应值 ----
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
  ### 反应式值存储用户选择 ----
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
  ## 下载示例数据 ----
  output$download_example_data <- downloadHandler(
    filename = function() {
      "example_data.csv"
    },
    content = function(file) {
      write.csv(example_data(), file)
    })
  output$download_example_group_info <- downloadHandler(
    filename = function() {
      "example_group_info.csv"
    },
    content = function(file) {
      write.csv(example_group(), file)
    })
  ## 读取数据 ----
  data <- reactive({
    if (input$data_source == "example") {
      data <- example_data()
    } else {
      req(input$data_file)
      df <- read.csv(input$data_file$datapath)
      rownames(df) <- df[, 1]
      data <- subset(df,select = -c(X))
    }
    na_info <- check_na_values(data)
    na_check_status(na_info)
    # 如果存在缺失值，显示警告
    if (!is.null(na_info) && na_info$has_na) {
      show_na_warning(na_info, "data file")
    }
    negative_info <- check_negative_values(data)
    negative_check_status(negative_info)
    # 如果存在负值，显示警告
    if (!is.null(negative_info) && negative_info$has_negative) {
      show_negative_warning(negative_info, "data file")
    }
    return(data)
  })

  ### data_group
  data_group <- reactive({
    if (input$data_source == "example") {
      data <- example_group()
    } else {
      req(input$group_file)
      df <- read.csv(input$group_file$datapath)
      rownames(df) <- df[, 1]
      data <- subset(df,select = -c(X))
      
      # 检查分组数据的缺失值
      na_info_group <- check_na_values(data)
      if (!is.null(na_info_group) && na_info_group$has_na) {
        show_na_warning(na_info_group, "data group info file")
      }
      negative_info <- check_negative_values(data)
      negative_check_status(negative_info)
      # 如果存在负值，显示警告
      if (!is.null(negative_info) && negative_info$has_negative) {
        show_negative_warning(negative_info, "data group info file")
      }
      return(data)
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
    # 记录日志并更新计数
    current_count <- log_run_event()
    run_count(current_count)
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
    
    # source("./R/data_correct_adaptive_model.R", local = TRUE)
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
    plot_stat_distribution(data = result_check()$rawdata,marker_list = result_check()$marker_list,
                           type = "pearson", statistic = "pvalue", alpha = 0.05)
  },height = 400,width = 500)
  output$correlation_r_pre_plot <- renderPlot({
    source("./R/plot_stat_distribution.R")
    req(result_check())
    plot_stat_distribution(data = result_check()$rawdata,marker_list = result_check()$marker_list,
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
    plot_stat_distribution(data = result_correct()$correct_data,marker_list = result_check()$marker_list,
                           type = "pearson", statistic = "pvalue", alpha = 0.05)
  },height = 400,width = 500)
  output$correlation_r_post_plot <- renderPlot({
    source("./R/plot_stat_distribution.R")
    req(result_check())  # 确保数据存在
    req(result_correct())
    plot_stat_distribution(data = result_correct()$correct_data,marker_list = result_check()$marker_list,
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
                                          displayNumbers = T,input_type = "correlation",
                                          title = "Expression Correlation Matrix(Erythrocyte markers)")
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
                                          displayNumbers = T,input_type = "correlation",
                                          title = "Expression Correlation Matrix(Coagulation markers)")
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
                                          displayNumbers = T,input_type = "correlation",
                                          title = "Expression Correlation Matrix(Platelet markers)")
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
      theme_classic(base_size = 14) +
      labs(y = "Coefficient of Variation", x = "Type",fill = "Type") +
      theme(
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 13),
        axis.title.y = element_text(size = 15,colour = "black",face = "bold"),
        axis.title.x = element_text(size = 15,colour = "black",face = "bold"),
        legend.title = element_text(size = 20,colour = "black",face = "bold"),
        legend.text = element_text(size = 13,colour = "black")
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
      theme_classic(base_size = 14) +
      labs(y = "Coefficient of Variation", x = "Type",fill = "Type") +
      theme(
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 13),
        axis.title.y = element_text(size = 15,colour = "black"),
        axis.title.x = element_text(size = 15,colour = "black"),
        legend.title = element_text(size = 15,colour = "black",face = "bold"),
        legend.text = element_text(size = 13,colour = "black")
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
    marker <- result_correct()$marker_list$erythrocyte
    source("./R/plot_protein_by_sample.R")
    if (length(marker) > 0) {
      plot_protein_by_sample(data = result_correct()$correct_data[rownames(result_correct()$correct_data)%in%marker,],
                             title = "Erythrocyte")
    }
  })
  
  #### platelet marker ----
  output$platelet_marker_post_plot <- renderPlot({
    req(result_check())  # 确保数据存在
    req(result_correct())
    marker <- result_correct()$marker_list$platelet
    source("./R/plot_protein_by_sample.R")
    if (length(marker) > 0) {
      plot_protein_by_sample(data = result_correct()$correct_data[rownames(result_correct()$correct_data)%in%marker,],
                             title = "Platelet")
    }
  })
  
  #### coagulation marker ----
  output$coagulation_marker_post_plot <- renderPlot({
    req(result_check())  # 确保数据存在
    req(result_correct())
    marker <- result_correct()$marker_list$coagulation
    source("./R/plot_protein_by_sample.R")
    if (length(marker) > 0) {
      plot_protein_by_sample(data = result_correct()$correct_data[rownames(result_correct()$correct_data)%in%marker,],
                             title = "Coagulation")
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
