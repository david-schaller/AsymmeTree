# -*- coding: utf-8 -*-

"""Species and gene tree visualization."""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


__author__ = 'David Schaller'


def assign_species_colors(tree):
    """Assign a unique color to each (extant) species in a species tree.
    
    Parameters
    ----------
    tree : Tree
        The species tree with uniquely labeled (non-loss) leaves.
    
    Returns
    -------
    dict
        A dictonary containing the labels of the extant species as keys and
        the assigned colors as values.
    """
    
    return _assign_cmap({v.label for v in tree.leaves() if v.event != 'L'})


def assign_gene_colors(tree, species_colors=None):
    """Assign a color to each (extant) gene in a gene tree according to its
    species.
    
    Parameters
    ----------
    tree : Tree
        The gene tree whose (non-loss) leaves have the 'reconc' attribute, i.e,
        the information to which species they belong.
    species_colors : dict, optional
        A dictonary containing the values of the 'color' attribute appearing
        in the tree as keys and assigned colors as values.
    
    Returns
    -------
    dict
        A dictonary containing the labels of the extant genes as keys and
        the assigned colors as values.
    """
    
    if species_colors is None:
        
        species_colors = _assign_cmap({v.reconc for v in tree.leaves()
                                       if v.event != 'L'})
    
    return {v.label: species_colors[v.reconc] for v in tree.leaves()
            if v.event != 'L'}


def assign_colors(species_tree, gene_tree):
    """Assign a color to each (extant) species and gene.
    
    Parameters
    ----------
    species_tree : Tree
        The species tree with uniquely labeled (non-loss) leaves.
    gene_tree : Tree
        The gene tree whose (non-loss) leaves have the 'reconc' attribute, i.e,
        the information to which species they belong, that furthermore appear
        as labels of the (non-loss) leaves in the species tree.
    
    Returns
    -------
    tuple of two dicts
        Dictonaries containing the labels of the extant species/genes as keys
        and the assigned colors as values.
    """
    
    species_colors = assign_species_colors(species_tree)
    
    return species_colors, assign_gene_colors(gene_tree,
                                              species_colors=species_colors)


def _assign_cmap(colors):
    
    color_dict = {}
    
    if len(colors) <= 10:
        cmap = plt.get_cmap('tab10')(np.arange(len(colors), dtype=int))
    else:
        cmap = plt.get_cmap('jet')(np.linspace(0, 1.0, len(colors)))
        
    for label, color in zip(colors, cmap):
        color_dict[label] = color
    
    return color_dict


def visualize(tree, color_dict=None, save_as=False, scale_symbols=1.0,
              fontsize='medium'):
    """Visualize a tree.
    
    Parameters
    ----------
    tree : Tree
        A tree whose nodes have the 'dist' attribute.
    color_dict : dict, optional
        A dictonary containing as values the assigned color for each label
        of the (non-loss) leaves.
    save_as : str, optional
        Save the image to a file. The default is False, in which case the
        image is not saved.
    scale_symbols : float, optional
        Adjust the size of the event type symbols. The default is 1.0.
    fontsize : str or float or int, optional
        Adjust the size of the text. The default is 'medium'.
    """
    
    vis = Visualizer(tree, color_dict=color_dict, scale_symbols=scale_symbols,
                     fontsize=fontsize)
    vis.draw(save_as=save_as)


class Visualizer:
    """Species and gene tree visualization."""
    
    def __init__(self, tree, color_dict=None,
                 scale_symbols=1.0, fontsize='medium',
                 species_info=False):
        """
        Parameters
        ----------
        tree : Tree
            A tree whose nodes have the 'dist' attribute.
        color_dict : dict, optional
            A dictonary containing as values the assigned color for each label
            of the (non-loss) leaves.
        scale_symbols : float, optional
            Adjust the size of the event type symbols. The default is 1.0.
        fontsize : str or float or int, optional
            Adjust the size of the text. The default is 'medium'.
        species_info : bool, False
            If True, write the reconciliation information behind the leaf
            label if available. The default is False.
        """
        
        self.tree = tree
        
        if color_dict:
            self.colors = color_dict
        else:
            self.colors = {v.label: 'white' for v in tree.preorder()}
        
        self.species_info = species_info
        
        self.symbolsize = 0.03 * scale_symbols
        self.fontsize = fontsize
        self.symbollw = 0.04
        self.leafs_per_vertical_unit = 15
        self.symbol_zorder = 3
        
        self.distance_dict = {}
        #self.divtime_dict {}
        self.leaf_counter = sum(1 for _ in tree.leaves())
        self.node_positions = {}
        
    
    def draw(self, distance_mode='evolutionary', save_as=False):
        """Draw the tree and optionally save it to file.
        
        Parameters
        ----------
        distance_mode : str, optional
            Information that is used for the length of the species edges.
            Currently only 'evolutionary' is supported, in which case the
            'dist' attribute is used.
        save_as : str, optional
            Save the image to a file. The default is False, in which case the
            image is not saved.
        """
        
        self.fig, self.ax = plt.subplots()
        self.ax.set_aspect('equal')
        self.ax.invert_yaxis()
        
        self.initial_traversal(distance_mode)
        # self.assign_colors()
        self.assign_positions()
        self.draw_edges()
        self.draw_nodes()
        
        self.ax.axvline(x=0, linewidth=1, color='grey', linestyle='--')
        self.ax.set_yticks([])
        self.ax.spines['top'].set_visible(False)
        self.ax.spines['left'].set_visible(False)
        self.ax.spines['right'].set_visible(False)
        self.ax.tick_params(axis='x', labelsize=self.fontsize)
        
        xmin, xmax = self.ax.get_xlim()
        ymin, ymax = self.ax.get_ylim()
        self.fig.set_size_inches(5*abs(xmax-xmin), 5*abs(ymax-ymin)+0.4)
        plt.tight_layout()
        if save_as:
            plt.savefig(save_as, dpi=300)
        else:
            plt.show()
        
    
    def initial_traversal(self, distance_mode):
        
        xmax = 0.0
        
        for v in self.tree.preorder():
            if distance_mode == 'evolutionary':
                if not v.parent:
                    self.distance_dict[v] = 0.0
                else:
                    self.distance_dict[v] = self.distance_dict[v.parent] + v.dist
                    if self.distance_dict[v] > xmax:
                        xmax = self.distance_dict[v]
            elif distance_mode == 'divergence time':
                raise ValueError('divergence time not yet implemented')
            else:
                raise ValueError(f"distance mode '{distance_mode}' not supported")
                
        self.ax.set_xlim(-0.1,xmax+0.5)
    
    
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
                if hasattr(v, 'transferred') and v.transferred:
                    color = 'red'
                else:
                    color = 'black'
                self.ax.plot([self.node_positions[v.parent][0],
                              self.node_positions[v][0]],
                             [self.node_positions[v][1],
                              self.node_positions[v][1]],
                     color=color,
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
            if not v.parent and not v.event:
                self.draw_root(*self.node_positions[v])
            elif v.children:
                if v.event == 'D':
                    self.draw_dupl(*self.node_positions[v])
                elif v.event == 'S':
                    self.draw_spec(*self.node_positions[v])
                elif v.event == 'H':
                    self.draw_hgt(*self.node_positions[v])
                elif v.event == 'GC':
                    self.draw_gc(*self.node_positions[v])
            else:
                if v.event == 'L':
                    self.draw_loss(*self.node_positions[v])
                else:
                    x, y = self.node_positions[v]
                    self.draw_leaf(x, y, color=self.colors[v.label])
                    if not self.species_info:
                        text = str(v.label)
                    elif hasattr(v, 'reconc'):
                        text = '{} <{}>'.format(v.label, v.reconc)
                    self.write_label(x+self.symbolsize+0.02, y, text)
                
    
    def draw_leaf(self, x, y, color='white', leftalign=False):
        
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
    
    
    def draw_root(self, x, y, rightalign=False):
        
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
        
        coord = np.asarray([[x-self.symbolsize/2,y],
                            [x+self.symbolsize/2,y-self.symbolsize/2],
                            [x+self.symbolsize/2,y+self.symbolsize/2]])

        inner = mpatches.Polygon(coord, closed=True,
                                 color='white', fill=True,
                                 zorder=self.symbol_zorder)
        self.ax.add_patch(inner)
        outer = mpatches.Polygon(coord, closed=True,
                                 color='black', fill=False,
                                 lw=self.symbolsize/self.symbollw,
                                 zorder=self.symbol_zorder)
        self.ax.add_patch(outer)
    
    
    def draw_gc(self, x, y):
        
        coord = []
        R = self.symbolsize/2
        for i in range(5):
            coord.append([x + 0.5 * R * np.cos(np.radians(0+i*72)),
                          y - 0.5 * R * np.sin(np.radians(0+i*72))])
            coord.append([x + R * np.cos(np.radians(36+i*72)),
                          y - R * np.sin(np.radians(36+i*72))])

        inner = mpatches.Polygon(coord, closed=True,
                                 color='black', fill=True,
                                 zorder=self.symbol_zorder)
        self.ax.add_patch(inner)
    
    
    def write_label(self, x, y, text):
        
        self.ax.text(x, y, text,
                     fontsize=self.fontsize,
                     horizontalalignment='left',
                     verticalalignment='center')
