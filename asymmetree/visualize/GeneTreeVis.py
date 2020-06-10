# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


__author__ = 'David Schaller'


class GeneTreeVis:
    
    def __init__(self, tree):
        
        self.tree = tree
        
        self.symbolsize = 0.03
        self.symbollw = 0.04
        self.leafs_per_vertical_unit = 15
        self.symbol_zorder = 3        
        
        print(tree.to_newick())
        
        self.distance_dict = {}
        self.colors = {}
        #self.divtime_dict {}
        self.leaf_counter = 0
        self.node_positions = {}
        
        self.draw()
        
    
    def draw(self):
        
        self.fig, self.ax = plt.subplots()
        self.ax.set_aspect('equal')
        self.ax.invert_yaxis()
        
        self.initial_traversal()
        self.assign_colors()
        self.assign_positions()
        self.draw_edges()
        self.draw_nodes()
        
        self.ax.axvline(x=0, linewidth=1, color='grey', linestyle='--')
        self.ax.set_yticks([])
        self.ax.spines["top"].set_visible(False)
        self.ax.spines["left"].set_visible(False)
        self.ax.spines["right"].set_visible(False)
        
        xmin, xmax = self.ax.get_xlim()
        ymin, ymax = self.ax.get_ylim()
        self.fig.set_size_inches(5*abs(xmax-xmin), 5*abs(ymax-ymin)+0.4)
        plt.tight_layout()
        plt.show()
        
    
    def initial_traversal(self, distance="evolutionary"):
        
        xmax = 0.0
        
        for v in self.tree.preorder():
            if distance == 'evolutionary':
                if not v.parent:
                    self.distance_dict[v] = 0.0
                else:
                    self.distance_dict[v] = self.distance_dict[v.parent] + v.dist
                    if self.distance_dict[v] > xmax:
                        xmax = self.distance_dict[v]
            elif distance == 'divtime':
                raise ValueError('divergence time not yet implemented')
            else:
                raise ValueError("distance mode '{}' not supported".format(distance))
            if not v.children:
                self.leaf_counter += 1
                if not v.is_loss() and v.color not in self.colors:
                    self.colors[v.color] = len(self.colors)
                
        self.ax.set_xlim(-0.1,xmax+0.5)
        
    
    def assign_colors(self):
        
        if len(self.colors) <= 10:
            cmap = plt.get_cmap('tab10')(np.arange(len(self.colors), dtype=int))
        else:
            cmap = plt.get_cmap('jet')(np.linspace(0, 1.0, len(self.colors)))
        for color_label, color in zip(self.colors.keys(), cmap):
            self.colors[color_label] = color
        print(self.colors)
    
    
    def assign_positions(self):
        
        ymax = (self.leaf_counter-1)/self.leafs_per_vertical_unit
        self.ax.set_ylim(ymax+0.1, -self.symbolsize*0.6)
        
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
                if v.is_loss():
                    self.draw_loss(*self.node_positions[v])
                else:
                    x, y = self.node_positions[v]
                    self.draw_leaf(x, y,
                                   color=self.colors[v.color])
                    self.write_label(x+self.symbolsize+0.02, y,
                                     str(v.label))
                
    
    def draw_leaf(self, x, y, color='white', leftalign=True):
        
        if leftalign:
            x += self.symbolsize/2
        fill = mpatches.Circle((x, y), self.symbolsize/2,
                               color=color, fill=True,
                               zorder=self.symbol_zorder)
        self.ax.add_patch(fill)
        outer = mpatches.Circle((x, y), self.symbolsize/2,
                                color='black', fill=False,
                                lw=self.symbolsize/self.symbollw,
                                zorder=self.symbol_zorder)
        self.ax.add_patch(outer)
        little = mpatches.Circle((x, y), self.symbolsize/8,
                                 color='black', fill=True,
                                 zorder=self.symbol_zorder)
        self.ax.add_patch(little)
        
    
    def draw_loss(self, x, y):
        
        self.ax.plot([x,x],
                     [y-self.symbolsize/2, y+self.symbolsize/2],
                     color='black',
                     linestyle='-', linewidth=1)
    
    
    def draw_root(self, x, y, rightalign=True):
        
        if rightalign:
            x -= self.symbolsize/2
        fill = mpatches.Circle((x, y), self.symbolsize/2,
                               color='white', fill=True,
                               zorder=self.symbol_zorder)
        self.ax.add_patch(fill)
        outer = mpatches.Circle((x, y), self.symbolsize/2,
                                color='black', fill=False,
                                lw=self.symbolsize/self.symbollw,
                                zorder=self.symbol_zorder)
        self.ax.add_patch(outer)
        little = mpatches.Circle((x, y), self.symbolsize/5,
                                 color='black', fill=False,
                                 lw=self.symbolsize/self.symbollw,
                                 zorder=self.symbol_zorder)
        self.ax.add_patch(little)
    
    
    def draw_spec(self, x, y):
        
        fill = mpatches.Circle((x, y), self.symbolsize/2,
                               color='black', fill=True,
                               zorder=self.symbol_zorder)
        self.ax.add_patch(fill)
        
    
    def draw_dupl(self, x, y):
        
        square = mpatches.Rectangle((x-self.symbolsize/2,
                                     y-self.symbolsize/2),
                                    width=self.symbolsize,
                                    height=self.symbolsize,
                                    color='white', fill=True,
                                    zorder=self.symbol_zorder)
        self.ax.add_patch(square)
        border = mpatches.Rectangle((x-self.symbolsize/2,
                                     y-self.symbolsize/2),
                                    width=self.symbolsize,
                                    height=self.symbolsize,
                                    color='black', fill=False,
                                    lw=self.symbolsize/self.symbollw,
                                    zorder=self.symbol_zorder)
        self.ax.add_patch(border)
    
    
    def draw_hgt(self, x, y):
        
        coord = np.asarray([[x,y-self.symbolsize/2],
                            [x+self.symbolsize/2,y+self.symbolsize/2],
                            [x-self.symbolsize/2,y+self.symbolsize/2]])

        inner = mpatches.Polygon(coord, closed=True,
                                 color='white', fill=True,
                                 zorder=self.symbol_zorder)
        self.ax.add_patch(inner)
        outer = mpatches.Polygon(coord, closed=True,
                                 color='black', fill=False,
                                 lw=self.symbolsize/self.symbollw,
                                 zorder=self.symbol_zorder)
        self.ax.add_patch(outer)
    
    
    def write_label(self, x, y, text):
        
        self.ax.text(x, y, text,
                     horizontalalignment='left',
                     verticalalignment='center')


if __name__ == '__main__':
    
    import asymmetree.treeevolve as te
    
    S = te.simulate_species_tree(6, planted=True, non_binary_prob=0.2)
    
    TGT_simulator = te.GeneTreeSimulator(S)
    TGT = TGT_simulator.simulate(DLH_rates=(1.0,1.0,1.0),
                                 dupl_polytomy=0.5)
    TGT = te.assign_rates(TGT, S)
    
    gtv = GeneTreeVis(TGT)