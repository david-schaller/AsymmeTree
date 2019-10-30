# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


__author__ = "David Schaller"
__copyright__ = "Copyright (C) 2019, David Schaller"


class GeneTreeVis:
    
    def __init__(self, tree):
        self.tree = tree
        
        self.symbolsize = 0.03
        self.symbollw = 0.04
        self.leafs_per_vertical_unit = 15
        
        
        print(tree.to_newick())
        
        self.distance_dict = {}
        #self.divtime_dict {}
        self.leaf_counter = 0
        self.node_positions = {}
        
        self.draw()
        
    
    def draw(self):
        self.fig, self.ax = plt.subplots()
        self.ax.set_aspect('equal')
        self.ax.invert_yaxis()
        
        self.initial_traversal()
        self.assign_positions()
        self.draw_edges()
        self.draw_nodes()
        
        plt.tight_layout()
        plt.show()
        
    
    def initial_traversal(self):
        
        xmax = 0.0
        
        for v in self.tree.preorder():
            if not v.parent:
                self.distance_dict[v] = 0.0
            else:
                self.distance_dict[v] = self.distance_dict[v.parent] + v.dist
                if self.distance_dict[v] > xmax:
                    xmax = self.distance_dict[v]
            if not v.children:
                self.leaf_counter += 1
                
        self.ax.set_xlim(-0.1,xmax+0.5)
                
    
    def assign_positions(self):
        
        ymax = (self.leaf_counter-1)/self.leafs_per_vertical_unit
        self.ax.set_ylim(ymax+0.1, -0.1)
        
        yposition = 0
        for v in self.tree.postorder():
            if not v.children:
                self.node_positions[v] = (self.distance_dict[v],
                                          yposition)
                yposition += 1/self.leafs_per_vertical_unit
            else:
                ymean = (self.node_positions[v.children[0]][1] +
                         self.node_positions[v.children[-1]][1])/2
                self.node_positions[v] = (self.distance_dict[v],
                                          ymean)
    
    def draw_edges(self):
        
        for v in self.tree.preorder():
            if v.parent:
                self.ax.plot([self.node_positions[v.parent][0],
                              self.node_positions[v][0]],
                             [self.node_positions[v][1],
                              self.node_positions[v][1]],
                     color='black',
                     linestyle='-', linewidth=1)
            if v.children:
                self.ax.plot([self.node_positions[v][0],
                              self.node_positions[v][0]],
                             [self.node_positions[v.children[0]][1],
                              self.node_positions[v.children[-1]][1]],
                     color='black',
                     linestyle='-', linewidth=1)
    
    
    def draw_nodes(self):
        for v in self.tree.preorder():
            if not v.parent:
                self.draw_root(*self.node_positions[v])
            elif v.children:
                if v.label == 'D':
                    self.draw_dupl(*self.node_positions[v])
                elif v.label == 'S':
                    self.draw_spec(*self.node_positions[v])
                if v.label == 'H':
                    self.draw_hgt(*self.node_positions[v])
            else:
                if v.label == '*':
                    self.draw_loss(*self.node_positions[v])
                    print(self.node_positions[v])
                else:
                    self.draw_leaf(*self.node_positions[v])
                
    
    def draw_leaf(self, x, y, color='white', leftalign=True):
        if leftalign:
            x += self.symbolsize/2
        fill = mpatches.Circle((x, y), self.symbolsize/2,
                               color=color, fill=True)
        self.ax.add_patch(fill)
        outer = mpatches.Circle((x, y), self.symbolsize/2,
                                color='black', fill=False,
                                lw=self.symbolsize/self.symbollw)
        self.ax.add_patch(outer)
        little = mpatches.Circle((x, y), self.symbolsize/8,
                                 color='black', fill=True)
        self.ax.add_patch(little)
        
    
    def draw_loss(self, x, y):
        print(x,y)
        self.ax.plot([x,x],
                     [y-self.symbolsize/2, y+self.symbolsize/2],
                     color='black',
                     linestyle='-', linewidth=1)
    
    
    def draw_root(self, x, y):
        fill = mpatches.Circle((x, y), self.symbolsize/2,
                               color='yellow', fill=True)
        self.ax.add_patch(fill)
        outer = mpatches.Circle((x, y), self.symbolsize/2,
                                color='black', fill=False,
                                lw=self.symbolsize/self.symbollw)
        self.ax.add_patch(outer)
        little = mpatches.Circle((x, y), self.symbolsize/5,
                                 color='black', fill=False,
                                 lw=self.symbolsize/self.symbollw)
        self.ax.add_patch(little)
    
    
    def draw_spec(self, x, y):
        fill = mpatches.Circle((x, y), self.symbolsize/2,
                               color='black', fill=True)
        self.ax.add_patch(fill)
        
    
    def draw_dupl(self, x, y):
        square = mpatches.Rectangle((x-self.symbolsize/2,
                                     y-self.symbolsize/2),
                                    width=self.symbolsize,
                                    height=self.symbolsize,
                                    color='white', fill=True)
        self.ax.add_patch(square)
        border = mpatches.Rectangle((x-self.symbolsize/2,
                                     y-self.symbolsize/2),
                                    width=self.symbolsize,
                                    height=self.symbolsize,
                                    color='black', fill=False,
                                    lw=self.symbolsize/self.symbollw)
        self.ax.add_patch(border)
    
    
    def draw_hgt(self, x, y):
        coord = np.asarray([[x,y-self.symbolsize/2],
                            [x+self.symbolsize/2,y+self.symbolsize/2],
                            [x-self.symbolsize/2,y+self.symbolsize/2]])

        inner = mpatches.Polygon(coord, closed=True,
                                 color='white', fill=True)
        self.ax.add_patch(inner)
        outer = mpatches.Polygon(coord, closed=True,
                                 color='black', fill=False,
                                 lw=self.symbolsize/self.symbollw)
        self.ax.add_patch(outer)


if __name__ == "__main__":
    
    import simulator.TreeSimulator as ts
    import simulator.TreeImbalancer as tm
    
    S = ts.build_species_tree(3, planted=True)
    
    TGT = ts.build_gene_tree(S, (1.0,1.0,1.0))
    TGT = tm.imbalance_tree(TGT, S, baseline_rate=1,
                            lognormal_v=0.2,
                            gamma_param=(0.5, 1.0, 2.2),
                            weights=(1/3, 1/3, 1/3),
                            copy_tree=False)
    
    gtv = GeneTreeVis(TGT)