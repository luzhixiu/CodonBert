import pydot

graph = pydot.Dot("my_graph", graph_type="digraph", bgcolor="yellow")

# Add nodes
my_node = pydot.Node("a", label="AGT")
graph.add_node(my_node)
# Or, without using an intermediate variable:
graph.add_node(pydot.Node("b", shape="circle",color='gray'))

# Add edges
help(pydot.Edge)
my_edge = pydot.Edge("a", "b", color="red",weight=3,label=3,penwidth = 5)

graph.add_edge(my_edge)
# Or, without using an intermediate variable:
graph.add_edge(pydot.Edge("b", "c", color="blue"))
output_graphviz_svg = graph.create_svg()
graph.write_png("output.png")
