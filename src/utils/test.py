import dendropy
import os
def test(tree,node):
	l = node.child_nodes()
	s = str()	
	for a in l:
		s += str(a.label)
		if a.label=="node5":
			node.clear_child_nodes()
			return "done"
	return s
		
