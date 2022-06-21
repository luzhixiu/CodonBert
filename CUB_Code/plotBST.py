from graphviz import Digraph
dot = Digraph()

#Part of plotting Code is from this github io: https://h1ros.github.io/posts/introduction-to-graphviz-in-jupyter-notebook/

class TreeNode(object):
    def __init__(self, x):
        self.val = x
        self.left = None
        self.right = None

#This funtion is used to merge two trees, not sure if I will need it, just keeping it for the future
class Solution(object):
    def mergeTrees(self, t1, t2):
        """
        :type t1: TreeNode
        :type t2: TreeNode
        :rtype: TreeNode
        """

        if t1 and t2:
            # Merge the values from two trees
            node = TreeNode(t1.val + t2.val)

            # Merge the left trees using the function recursively
            node.left = self.mergeTrees(t1.left, t2.left)

            # Merge the right trees using the function recursively
            node.right = self.mergeTrees(t1.right, t2.right)

            return node

        else:
            return t1 or t2

#scrath input rscuList[1.4,0.5,1.7]
def createBSTFromRSCUList(rscuList):
    tree=TreeNode(0)
    root=tree
    for rscu in rscuList:
        if rscu>=1:
            tree.left=TreeNode(1)
            tree=tree.left
        else:
            tree.right=TreeNode(1)
            tree=tree.right
    return root
# rscuList=[1.4,0.5,1.7,0.8]
# tree=createBSTFromRSCUList(rscuList)
# print(tree.left.right.val)

def createBSTfromMultipleRSCUList(rscuLists):
    tarzan=None
    s=Solution()
    for rscuList in rscuLists:
        tree=createBSTFromRSCUList(rscuList)
        tarzan=s.mergeTrees(tarzan,tree)
    return tarzan
# rscuLists=[[1.4,0.5,1.7,0.8],[1.4,0.5,1.7]]
# tarzan=createBSTfromMultipleRSCUList(rscuLists)
# print(tarzan.left.right.val)

def visualize_tree(tree):
    def add_nodes_edges(tree, dot=None):
        # Create Digraph object
        if dot is None:
            dot = Digraph()
            dot.node(name=str(tree), label=str(tree.val))

        # Add nodes
        if tree.left:
            dot.node(name=str(tree.left), label=str(tree.left.val),color="red")
            dot.edge(str(tree), str(tree.left))
            dot = add_nodes_edges(tree.left, dot=dot)

        if tree.right:
            dot.node(name=str(tree.right), label=str(tree.right.val),color="green")
            dot.edge(str(tree), str(tree.right))
            dot = add_nodes_edges(tree.right, dot=dot)
        return dot

    # Add nodes recursively and create a list of edges
    dot = add_nodes_edges(tree)

    # Visualize the graph
    # display(dot)
    dot.render(directory='doctest-output', view=True)
    return dot
