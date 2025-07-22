from graphviz import Digraph
from IPython.display import display, Image


def Visualize(tree, filename = 'binary_tree'):
    ''' Visualizes a binary tree using Graphviz.
        Input tree format: [node_list, [left_child, right_child]] where node_list is a list and left_child, right_child are similar subtrees or None.
    '''

    dot = Digraph(comment = 'Binary Tree')
    dot.attr(rankdir = 'TB')  # Top to bottom layout
    dot.attr('node', shape = 'plaintext')  # Nodes as rectangles for list display
    dot.attr('edge', arrowhead='none') # no arrowhead

    def add_node(current_tree, parent_id=None, edge_label=None):
        if current_tree is None or not current_tree:
            return None
        
        # Extract node list and children
        node_list = current_tree[0]
        children = current_tree[1] if len(current_tree) > 1 else [None, None]
        
        # Create a unique node ID
        node_id = str(id(current_tree))
        
        # Convert node list to string for display
        node_label = str(node_list)
        dot.node(node_id, node_label)
        
        # Connect to parent if exists
        if parent_id:
            dot.edge(parent_id, node_id, label = edge_label)
        
        try:
            # Recursively process left and right children
            add_node(children[0], node_id) if children[0] else None
            add_node(children[1], node_id) if children[1] else None
        except IndexError:
            pass
        
        return node_id

    # Start processing from root
    add_node(tree)
    
    # # Render and save the graph
    dot.render(filename, format = 'png', cleanup = True)
    print(f"Tree visualization saved as {filename}.png\n")

    # # display the graph in notebook
    # dot.format = 'png'
    # display(Image(dot.pipe()))

