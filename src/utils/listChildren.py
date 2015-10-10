#!/usr/bin/env python

import dendropy
import sys
import os
import copy

def listChildren(Tree conTree):
	tree_set=set()
	for e in conTree.postorder_node_iter():
		if e.is_leaf():
		tree_set.add(e)

	return(tree_set)	
