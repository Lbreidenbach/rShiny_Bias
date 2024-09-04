library(shiny)
library(DiagrammeR)
library(HydeNet)
library(rjags)
library(MatchIt)
library(ggplot2)
library(plyr)

library(dplyr)
library(gtable)
library(grid)
library(GGally)
library(ggridges)
library("reshape2")
library(gridExtra)
source("plot_dag_ui2.R")
library(purrr)


set_p = function(p,model){
  p2 = log(p/(1-p))
  b0 = p2-model
  return(b0)
}

textInput3<-function (inputId, label, value = "",...) {
  div(style="display:inline-block",
      tags$label(label, `for` = inputId),
      tags$input(id = inputId, type = "text", value = value,...))
}

node_ui <- function(id) {
  ns <- NS(id)
  fluidRow(
    column(5,
           selectInput(ns("node_dist"),
                       label = paste0(id,"'s distribution"),
                       choices = c("binary", "continuous"))
    ),
    column(7,
           uiOutput(ns("ui_placeholder"))

    )
  )
}

node_server <- function(input, output, session) {
  return_value <- reactive({list(distribution = input$node_dist, prevalence = input$prev,
                                 mean = input$mean,
                                 std_dev = input$stdev
  )})
  ns <- session$ns
  node_name = sub("-", "", session$ns(""))
  output$ui_placeholder <- renderUI({
    type <- req(input$node_dist)
    if(type == "binary") {
      numericInput(ns("prev"), paste0(node_name,"'s prevalence (between 0-1):"), 0.5, min = 0, max = 1, step = 0.05)

    } else if (type == "continuous") {
      tagList(numericInput(ns("mean"), paste0("Mean of ", node_name,":"),0),
              numericInput(ns("stdev"), paste0("Std dev of ", node_name,":"), 1),
      )


    }
  })

  ## if we later want to do some more sophisticated logic
  ## we can add reactives to this list
  # list(return_value = reactive({input$node_dist}), try = reactive({input$inner_element}))
  return_value
}

ui <- fluidPage(
  titlePanel(h1("Bias Simulator", align ="center")),
  tabsetPanel(
    tabPanel("Set Bayesian Network", fluid = TRUE,
             fluidRow(column(12, mainPanel(HTML("<b>Input the directed acyclic graph formula in the text box below and follow this format:</b><br>
                                  1. Start formula with '~'<br>
                                  2. Write which variables cause others like this: 'effect|cause1*cause2*cause_n' or like this 'effect|cause1 + effect|cause2 + effect|cause_n'. <br>
                                          &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&#8226; if a causes b causes c, write 'b|a + c|b' or 'c|b + b|a'<br>
                                          &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&#8226; if a causes b, and c has no effect on either, write 'b|a + c'<br>
                                          &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&#8226; there must be at least one variable that affects another<br>
                                          &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&#8226; nodes can't have names with spaces<br>
                                          &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&#8226; formula can't be cyclical. cyclical formaulas include 'b|a + a|b' or 'c|b + b|a + a|c'<br>
                                  3. Once the formula is entered, fill out the additional node info on the side"),
                                           textInput("dag_text", label = h4("DAG formula box"),
                                                     value = "~exposure|confounder + outcome|confounder*exposure + collider|exposure*outcome",
                                                     width = "100%"),

                                           width = 12)
             )
             ),
             sidebarLayout(position = "right",

                           fluidRow(column(6, sidebarPanel(h4("---Set Distributions---", align = "center"),
                                                           uiOutput("node_box"),
                                                           div(id="placeholder"),
                                                           h4("---Set Beta Values---", align = "center"),
                                                           uiOutput("beta_box"),

                                                           width = 12
                           )),
                           ),
                           column(6, align = "center", mainPanel(h4("Dag plot"),
                                                                 grVizOutput("value"),


                                                                 width = 12
                           )),



             ),


    ),
    tabPanel("Set anaylsis", fluid = TRUE,
             sidebarLayout(
               sidebarPanel(h4("---Set Analysis---", align = "center"),
                            uiOutput("exp"),
                            uiOutput("out"),
                            uiOutput("sel"),
                            uiOutput("sel_op"),
                            uiOutput("conf"),
                            uiOutput("conf_op"),
                            #numericInput("n_data", "Number of samples per dataset \n(takes a minute to calculate when n>100,000 n>1,000,000 not recommended)", 10000, min = 5, step = 1),
                            #numericInput("iteration", "Number of datasets created and analyzed \n(i>300 not recommended)", 100, min = 2, step = 1),
                            actionButton("sim", "Simulate! (takes a few seconds)")
               ),
               mainPanel(fluidRow(column(12, mainPanel(h4("Plot appears here once analysis is set and the simulate button is pressed ", align ="center"),
                                                       plotOutput("bn_results"),
                                                       h5("Code for running set analysis", align ="center"),
                                                       verbatimTextOutput("analysis_info"),
                                                       width = 12)
               )
               )

               )
             ),
             fluidRow(column(12, mainPanel(h3("Code for your specified Bayesian Network", align ="center"),
                                           p("If the Hydenet library is in use, this code can be copy/pasted from this box to get a working dag object and results matrix"),
                                           verbatimTextOutput("bn_info"),
                                           width = 12)
             )
             )
    )
  ),


)

# Define server logic ----
server <- function(input, output, session) {

  handler = reactiveVal(list())
  node_ids = reactiveVal()
  beta_ids = reactiveVal()
  dag_code = reactiveVal()

  output$value = renderGrViz({ dag_ui(paste(input$dag_text)) })

  observeEvent(input$dag_text, {
    #beta interface
    beta_id = get_arrows(input$dag_text)
    output$beta_box = renderUI(
      map(beta_id, ~textInput3(.x, paste0( .x), value = 0))

    )

    #node distribution interface
    x = get_nodes(input$dag_text)
    ui_num = length(x)
    new_id = unlist(lapply(c(1:ui_num), function(y) paste(x[y])))

    output$node_box = renderUI({
      map(new_id, ~node_ui(.x))
    })

    node_ids(x)
    beta_ids(beta_id)

    handler_list <- isolate(handler())

    new_handler = lapply(1:ui_num, function(y)
      callModule(node_server, new_id[y])
    )
    handler_list <-unlist(new_handler)
    names(handler_list) <- new_id
    handler(handler_list)




  })

  output$exp = renderUI({
    x = node_ids()

    fluidRow(column(12,
                    selectInput("exp", "Exposure", choices = x)
    )
    )

  })

  output$out = renderUI({
    x = node_ids()
    fluidRow(column(12,
                    selectInput("out", "Outcome", choices = x[!x==input$exp])
    )
    )

    # fluidRow(column(12,
    #                 selectInput("out", "Outcome", choices = handler_df(handler)[!handler_df(handler)==input$exp])
    # )
    # )

  })

  output$sel = renderUI({
    if(length(handler_df(handler)[!handler_df(handler)==input$exp & !handler_df(handler)==input$out])==0){
      fluidRow(column(12,
                      selectInput("selection", "Choose a node that represents selection bias? (Node must be binary)", choices = c( "No"), selected = "No")
      )
      )

    }else{
      fluidRow(column(12,
                      selectInput("selection", "Choose a node that represents selection bias? (Node must be binary)", choices = c("Yes", "No"), selected = "No")
      )
      )

    }

    # if(length(handler_df(handler)[!handler_df(handler)==input$exp & !handler_df(handler)==input$out])==0){
    #   input$selection = "No"
    # }

  })

  output$sel_op = renderUI({
    if(input$selection == "Yes"){

      fluidRow(column(12,
                      radioButtons("sel", "Selection Node", choices = handler_df(handler)[!handler_df(handler)==input$exp & !handler_df(handler)==input$out])
      )
      )

    }else{

    }



  })

  output$conf = renderUI({
    fluidRow(column(12,
                    selectInput("adjust", "Adjust for confounders?", choices = c("Yes", "No"), selected = "No")
    )
    )

  })

  output$conf_op = renderUI({
    x = node_ids()
    if(input$adjust == "Yes"){
      if(input$selection == "Yes"){
        fluidRow(column(12,
                        checkboxGroupInput("conf", "Confounders", x[!x==input$exp & !x==input$out & !x==input$sel])
        )
        )
      }else{
        fluidRow(column(12,
                        checkboxGroupInput("conf", "Confounders", x[!x==input$exp & !x==input$out])
        )
        )

      }

    }else{

    }

  })

  output$bn_info = renderPrint({
    x = node_ids()
    beta_id = beta_ids()

    test_df = data.frame(unlist(lapply(handler(), function(handle) {
      handle()[["distribution"]]
    })), holder = 1)

    test_df$prevalence = rep(0.5, length(x))
    test_df$mean = rep(0, length(x))
    test_df$std_dev = rep(1, length(x))

    test_df$names = x

    if("binary" %in% test_df[[1]]){
      prev_vals = test_df[,3] = unlist(lapply(handler(), function(handle) {
        handle()[["prevalence"]]
      }))
      test_df[names(prev_vals), 3] = prev_vals

    }

    if("continuous" %in% test_df[[1]]){
      mean_vals = unlist(lapply(handler(), function(handle) {
        handle()[["mean"]]
      }))
      test_df[names(mean_vals), 4] = mean_vals

      std_vals = unlist(lapply(handler(), function(handle) {
        handle()[["std_dev"]]
      }))
      test_df[names(std_vals), 5] = std_vals
    }



    colnames(test_df) = c("distribution", "holder", "prevalence", "mean", "std_dev", "names")
    beta_df = data.frame(parent = unlist(lapply(beta_id, function(x){strsplit(x, " -> ")[[1]][1]})),
                         child = unlist(lapply(beta_id, function(x){strsplit(x, " -> ")[[1]][2]})),
                         beta_vals = unlist(lapply(beta_id, function(x) input[[paste(x)]]))
    )



    bi_nodes_df = test_df[test_df$distribution=="binary",]
    cont_nodes_df = test_df[test_df$distribution=="continuous",]

    parent_list = lapply(c(1:length(rownames(test_df))), function(y) beta_df[beta_df[2]==rownames(test_df)[[y]], c(1,3)]) #child == y from the apply functions
    names(parent_list) = rownames(test_df)

    node_info_list = lapply(c(1:length(rownames(test_df))), function(y) test_df[rownames(test_df)[[y]], c(1:6)])
    #node_info_list[name] if you enter name like this it'll return info without indexing issues
    names(node_info_list) = rownames(test_df)

    #parent_list[[name]][1], get parent names
    #parent_list[[name]][2], get beta values
    #node_info_list[[name]][1,3,4,5,6] get node 1.) distribution, 3.) prev, 4.) mean, 5.)std_dev 6.)name

    #make i node name, not index
    get_glm = function(i){
      if(nrow(parent_list[[i]]) == 0){
        return("")
      }
      out_glm = paste(unlist(lapply(c(1:nrow(parent_list[[i]])), function(y){
        paste0(
          parent_list[[i]][y,2],
          '," * ',
          parent_list[[i]][y,1],
          ' + ",')
      })), collapse = " " )

      return(out_glm)
    }
    get_linker = function(i){
      if(nrow(parent_list[[i]])==0){
        return(set_p(node_info_list[[i]][[3]], 0))
      }
      set_p(node_info_list[[i]][[3]],
            sum(unlist(lapply(c(1:length(parent_list[[i]][[1]])), function(y){
              as.numeric(parent_list[[i]][y,2]) *
                get_type(i,y)
            }))))
    }
    get_type =function(i,y){
      if(nrow(parent_list[[i]])==0){
        return(0)
      }
      if(node_info_list[[parent_list[[i]][y,1]]][["distribution"]]=="binary"){
        col_num = "prevalence"
      }else if(node_info_list[[parent_list[[i]][y,1]]][["distribution"]]=="continuous"){
        col_num = "mean"
      }
      node_info_list[[parent_list[[i]][y,1]]][[col_num]]

    }
    identity_link = function(i){
      if(nrow(parent_list[[i]])==0){
        return(set_p(node_info_list[[i]][[3]], 0))
      }
      sum(unlist(lapply(c(1:length(parent_list[[i]][[1]])), function(y){
        as.numeric(parent_list[[i]][y,2]) *
          get_type(i,y) *-1
      })))
    }

    if(nrow(bi_nodes_df)>0){
      bi_code = unlist(lapply(rownames(bi_nodes_df), function(i) {paste0(
        '\ndag = setNode(dag, ',
        node_info_list[[i]][6],
        ', nodeType = "dbern", prob = paste0("ilogit(",',
        paste(get_glm(i), get_linker(i), collapse = " "),',")"))')}

      ))
    }else{
      bi_code = NULL
    }

    if(nrow(cont_nodes_df)>0){
      cont_code = unlist(lapply(rownames(cont_nodes_df), function(i)
      {paste0('\ndag = setNode(dag, ',
              node_info_list[[i]][6], ', nodeType = "dnorm", mu = paste0(',
              node_info_list[[i]][4],
              " + ",
              paste(get_glm(i), identity_link(i), collapse = " "),
              "), tau = ",
              (1/sqrt(cont_nodes_df[i,5])),
              ")")
      }))
    }else{
      cont_code = NULL
    }

    code_dag = c(paste0("dag = HydeNetwork("
                        ,input$dag_text,")"),
                 bi_code,
                 cont_code, "\n","\n"
    )

    dag_code(code_dag)
    cat(code_dag)

    # lapply(c(1:length(x)), function(i) length(parent_list[[i]][[1]]))

  })

  user_code = eventReactive(input$sim, {
    x = node_ids()
    beta_id = beta_ids()
    beta_df = data.frame(parent = unlist(lapply(beta_id, function(x){strsplit(x, " -> ")[[1]][1]})),
                         child = unlist(lapply(beta_id, function(x){strsplit(x, " -> ")[[1]][2]})),
                         beta_vals = unlist(lapply(beta_id, function(x) input[[paste(x)]]))
    )

    set_ate = as.numeric(beta_df[beta_df$parent == input$exp & beta_df$child == input$out, "beta_vals"])
    if(length(set_ate)==0){
      set_ate=0
    }

    #Current Bug: if there's no valid selection node and the user has it set to Yes, then sb is ""
    if(input$selection == "Yes"){
      sb = paste0('"', input$sel, '"')
    }else{
      sb = format(NULL)
    }



    if(input$adjust == "Yes"){
      if(length(input$conf)==1){
        confounder = paste0('"', input$conf, '"')

      }else{
        confounder = list(format(input$conf))
      }

    }else{
      confounder = format(NULL)
    }

    outcome = paste0(
      '\nrun_1 = varied_runs(', 100,
      ', dag, exposure = "', input$exp,
      '" , outcome = "', input$out ,
      '" , covariates = ', confounder,
      ' , sb = ', sb,
      ', n = ', 20000,
      ')')
    user_code = capture.output(cat(c(dag_code(), outcome)))

  })
  output$bn_results = renderPlot({

    run_code(user_code())

  })
  output$analysis_info = renderPrint({
    x = node_ids()
    beta_id = beta_ids()
    beta_df = data.frame(parent = unlist(lapply(beta_id, function(x){strsplit(x, " -> ")[[1]][1]})),
                         child = unlist(lapply(beta_id, function(x){strsplit(x, " -> ")[[1]][2]})),
                         beta_vals = unlist(lapply(beta_id, function(x) input[[paste(x)]]))
    )

    set_ate = as.numeric(beta_df[beta_df$parent == input$exp & beta_df$child == input$out, "beta_vals"])
    if(length(set_ate)==0){
      set_ate=0
    }

    #Current Bug: if there's no valid selection node and the user has it set to Yes, then sb is ""
    if(input$selection == "Yes"){
      sb = paste0('"', input$sel, '"')
    }else{
      sb = format(NULL)
    }



    if(input$adjust == "Yes"){
      if(length(input$conf)==1){
        confounder = paste0('"', input$conf, '"')

      }else{
        confounder = list(format(input$conf))
      }

    }else{
      confounder = format(NULL)
    }

    outcome = paste0(
      'run_1 = varied_runs(', 100,
      ', dag, exposure = "', input$exp,
      '" , outcome = "', input$out ,
      '" , covariates = ', confounder,
      ' , sb = ', sb,
      ', n = ', 10000,
      ')')
    cat(outcome)

  })



}

# Run the app ----
shinyApp(ui = ui, server = server)
