library(DiagrammeR)

# nodes and edges files

label1<-c("mvrnorm\n_rml","generator\n_graph","generator_param","generator_data","generator\n_ZpZhat","graph",
         "diagnostics","save\n_params","save_scores","na.omit.\nlist","compare\n_methods","record",
         "bootstrap\n_summary","simu")
label2<-c("F_\nNegLikelihood","F_NegGradient","SetLambda","FitBetaStatic","Laplacian","SumTree","Kirshner",
          "Edge\nProba","TreeGGM")
label<-c(label1,label2)
nodes<-tibble(labels=label,id=1:length(label))
nodes<-cbind(nodes,type=c(rep("Simu",length(label1)),rep("Tree",length(label2))))

edges<-tibble(from=1,to=4) %>%
  rbind(c(2,8)) %>%
  rbind(c(2,14)) %>%
  rbind(c(3,8)) %>%
  rbind(c(3,14)) %>%
  rbind(c(4,11)) %>%
  rbind(c(5,13)) %>%
  rbind(c(7,6)) %>%
  rbind(c(8,14)) %>%
  rbind(c(9,13)) %>%
  rbind(c(10,9)) %>%
  rbind(c(11,13)) %>%
  rbind(c(12,13)) %>%
  rbind(c(12,14)) %>%
  rbind(c(12,8)) %>%
  rbind(c(13,14)) %>%
  rbind(c(14,7)) %>%
  rbind(c(15,18)) %>%
  rbind(c(17,16)) %>%
  rbind(c(19,20)) %>%
  rbind(c(19,21)) %>%
  rbind(c(20,22))%>%
  rbind(c(21,16))%>%
  rbind(c(20,15))%>%
  rbind(c(22,18)) %>%
  rbind(c(18,23))%>%
  rbind(c(21,23))%>%
  rbind(c(23,11)) %>%
  rbind(c(16,18))


a<-create_graph() %>%
  add_nodes_from_table(
    table = nodes,
    label_col = labels,type_col=type)

b<- a %>%
  add_edges_from_table(
    table = edges,
    from_col = from,
    to_col = to,
    from_to_map = id_external)%>%
  set_edge_attrs(edge_attr = color,
                 value = "gray27") %>%
  set_edge_attrs(edge_attr = penwidth,value = 1.3) %>%
  set_edge_attrs(edge_attr = arrowsize,value = 1.2) %>%
  set_node_attrs(node_attr = fontcolor,
                 value = "black") %>%
  set_node_attrs(node_attr = fontsize,
                 value = 22) %>%
  set_node_attrs(node_attr = width,
                 value = 0.6) %>%
  select_nodes(
    conditions = type == "Simu") %>%
  set_node_attrs_ws(
    node_attr = fillcolor,
    value = "steelblue1") %>%
  invert_selection() %>%
  set_node_attrs_ws(
    node_attr = fillcolor,
    value = "olivedrab2") %>%
  delete_node(19) %>%
  delete_node(10)

render_graph(b,layout="neato")

tkid <- tkplot(b) #tkid is the id of the tkplot that will open
l <- tkplot.getcoords(tkid) # grab the coordinates from tkplot
plot(net_adj, layout=l)
