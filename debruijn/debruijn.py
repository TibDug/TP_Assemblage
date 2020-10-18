#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Perform assembly based on debruijn graph."""

from operator import itemgetter
from random import randint
import argparse
import os
import sys
import statistics
import random
import networkx as nx
#import matplotlib.pyplot as plt
random.seed(9001)

__author__ = "Thibault Dugauquier"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Thibault Dugauquier"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Thibault Dugauquier"
__email__ = "thibault.dug@gmail.com"
__status__ = "Developpement"

#==============================================================
# Main program
#==============================================================

def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=21, help="K-mer size (default 21)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file")
    return parser.parse_args()


def read_fastq(fastq_file):
    """prend un seul argument correspondant au fichier fastq et retourne un \
    générateur de séquences"""
    with open(fastq_file) as file_in:
        for ligne in enumerate(file_in):
            yield next(file_in)[:-1]
            next(file_in)
            next(file_in)


def cut_kmer(seq, kmer_size):
    """prend une séquence, une taille de k-mer et retourne un générateur de k-mer"""
    for nuc in range(len(seq[:-kmer_size+1])):
        yield seq[nuc:nuc+kmer_size]


def build_kmer_dict(fastq_file, kmer_size):
    """prend un fichier fastq, une taille k- mer et retourne un dictionnaire \
    ayant pour clé le k-mer et pour valeur le nombre d’occurrence de ce k-mer"""
    kmer_dict = {}
    for seq in read_fastq(fastq_file):
        for kmer in cut_kmer(seq, kmer_size):
            if kmer in kmer_dict:
                kmer_dict[kmer] += 1
            else:
                kmer_dict[kmer] = 1
    return kmer_dict


def build_graph(kmer_dict):
    """prend en entrée un dictionnaire de k-mer et crée l’arbre de k-mers \
    préfixes et suffixes"""
    graph = nx.DiGraph()
    for kmer, poids in kmer_dict.items():
        graph.add_node(kmer[:-1])
        graph.add_node(kmer[1:])
        graph.add_edge(kmer[:-1], kmer[1:], weight = poids)
    return graph


def get_starting_nodes(graph):
    """prend en entrée un graphe et retourne une liste de noeuds d’entrée"""
    starting_nodes = []
    for node in graph.nodes:
        if list(graph.predecessors(node)) == []:
            starting_nodes.append(node)
    return starting_nodes


def get_sink_nodes(graph):
    """prend en entrée un graphe et retourne une liste de noeuds de sortie"""
    sink_nodes = []
    for node in graph.nodes:
        if list(graph.successors(node)) == []:
            sink_nodes.append(node)
    return sink_nodes


def get_contigs(graph, starting_nodes, sink_nodes):
    """prend un graphe, une liste de noeuds d’entrée et une liste de \
    sortie et retourne une liste de tuple(contig, taille du contig)"""
    contigs_list = []
    for starting_node in starting_nodes:
        for sink_node in sink_nodes:
            for paths in nx.all_simple_paths(graph, starting_node, sink_node):
                contig = paths[0]
                for path in paths[1:]:
                    contig = contig + path[-1]
                contigs_list.append((contig, len(contig)))
    return contigs_list


def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))


def save_contigs(contigs_list, output_file):
    """prend une liste de tuple (contig, taille du contig) et un nom \
    de fichier de sortie et écrit un fichier de sortie contenant les \
    contigs selon le format fasta"""
    with open(output_file, "w") as file_out:
        for indice_contig, contig in enumerate(contigs_list):
            file_out.write(">contig_" + str(indice_contig) + " len=" + str(contig[1]) + \
                "\n" + fill(contig[0]) +"\n")


def std(val_list):
    """prend une liste de valeur, qui retourne l’écart type"""
    return statistics.stdev(val_list)


def path_average_weight(graph, path):
    """prend un graphe et un chemin et qui retourne un poids moyen"""
    weight_list = []
    for indice_node in range(len(path)-1):
        weight_list.append(graph.get_edge_data(path[indice_node], path[indice_node+1])['weight'])
    return statistics.mean(weight_list)


def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    """prend un graphe et une liste de chemin, la variable booléenne \
    delete_entry_node pour indiquer si les noeuds d’entrée seront supprimés \
    et la variable booléenne delete_sink_node pour indiquer si les noeuds de \
    sortie seront supprimés et retourne un graphe nettoyé des chemins indésirables"""
    for path in path_list:
        for indice_node in range(len(path)):
            if (indice_node == 0 and not delete_entry_node) or \
                (indice_node == len(path)-1 and not delete_sink_node):
                continue
            graph.remove_node(path[indice_node])
    return graph


def select_best_path(graph, path_list, path_size_list, average_weight_list, \
    delete_entry_node = False, delete_sink_node = False):
    """prend un graphe, une liste de chemin, une liste donnant la longueur \
    de chaque chemin, une liste donnant le poids moyen de chaque chemin, \
    delete_entry_node pour indiquer si les noeuds d’entrée seront supprimés \
    et delete_sink_node pour indiquer si les noeuds de sortie seront supprimés \
    et retourne un graphe nettoyé des chemins indésirables"""
    
    """
    if average_weight_list[0] > average_weight_list[1]:
        remove_paths(graph, [path_list[1]], delete_entry_node, delete_sink_node)
    elif average_weight_list[0] < average_weight_list[1]:
        remove_paths(graph, [path_list[0]], delete_entry_node, delete_sink_node)
    else:
        if path_size_list[0] > path_size_list[1]:
            remove_paths(graph, [path_list[1]], delete_entry_node, delete_sink_node)
        elif path_size_list[0] < path_size_list[1]:
            remove_paths(graph, [path_list[0]], delete_entry_node, delete_sink_node)
        else:
            remove_paths(graph, [path_list[random.randint(0,1)]], delete_entry_node, \
                delete_sink_node)
    return graph
    """
    del_path_list = []
    while len(path_list) != 1:
        if average_weight_list[0] > average_weight_list[1]:
            del_path_list.append(path_list[1])
            del path_list[1]
            del average_weight_list[1]
            del path_size_list[1]
        elif average_weight_list[0] < average_weight_list[1]:
            del_path_list.append(path_list[0])
            del path_list[0]
            del average_weight_list[0]
            del path_size_list[0]
        else:
            if path_size_list[0] > path_size_list[1]:
                del_path_list.append(path_list[1])
                del path_list[1]
                del average_weight_list[1]
                del path_size_list[1]
            elif path_size_list[0] < path_size_list[1]:
                del_path_list.append(path_list[0])
                del path_list[0]
                del average_weight_list[0]
                del path_size_list[0]
            else:
                del_value = random.randint(0,1)
                del_path_list.append(path_list[del_value])
                del path_list[del_value]
                del average_weight_list[del_value]
                del path_size_list[del_value]
    graph = remove_paths(graph, del_path_list, delete_entry_node, delete_sink_node)
    return graph
    
    '''
    ##average_weight_list[:average_weight_list.index(max_weight_indice)] + average_weight_list[average_weight_list.index(max_weight_indice)+1:]
    max_weight_indice = average_weight_list.index(max(average_weight_list))
    #max_size_indice = max(path_size_list)
    best_indice = max_weight_indice
    equal_weight_indice_list = [max_weight_indice]
    for weight_indice, weight in enumerate(average_weight_list):
        if weight_indice != max_weight_indice:
            if std([average_weight_list[max_weight_indice], weight]) > 0:
                equal_weight_indice_list.append(weight_indice)
    if len(equal_weight_indice_list) s> 1:
        sizes_to_compare = []
        for size_indice in equal_weight_indice_list:
            sizes_to_compare.append(path_size_list[size_indice])
        max_size_indice = path_size_list.index(max(sizes_to_compare))
        best_indice = max_size_indice
        equal_size_indice_list = [max_size_indice]
        for size_indice, size in enumerate(path_size_list):
            if size_indice != max_size_indice and size_indice in equal_weight_indice_list:
                if std([path_size_list[max_size_indice], size]) > 0:
                    equal_size_indice_list.append(size_indice)
        if len(equal_size_indice_list) > 1:
            best_indice = path_size_list.index(equal_size_indice_list[random.randint(0,len(equal_size_indice_list)-1)])
    remove_paths(graph, path_list[:best_indice] + path_list[best_indice+1:], delete_entry_node, delete_sink_node)
    return graph
    '''
            
    
    
            


def solve_bubble(graph, start_node, stop_node):
    """prend un graphe, un noeud ancêtre, un noeud descendant et retourne un graph \
    nettoyé de la bulle se trouvant entre ces deux noeuds"""
    path_list = []
    path_size_list = []
    average_weight_list = []
    for path in nx.all_simple_paths(graph, start_node, stop_node):
        path_list.append(path)
        path_size_list.append(len(path))
        average_weight_list.append(path_average_weight(graph, path))
    graph = select_best_path(graph, path_list, path_size_list, average_weight_list)
    return graph
    


def simplify_bubbles(graph) :
    """prend un graphe et retourne un graphe sans bulle"""
    for starting_node in get_starting_nodes(graph):
        for sink_node in get_sink_nodes(graph):
            paths = list(nx.all_simple_paths(graph, starting_node, sink_node))
            while len(paths) > 1:
                for indice_node, node in enumerate(paths[0]):
                    if not node in paths[1]:
                        node_begin_bubble = paths[0][indice_node-1]
                        break
                for node in paths[0][paths[0].index(node_begin_bubble)+1:]:
                    if node in paths[1]:
                        node_end_bubble = node
                        break
                graph = solve_bubble(graph, node_begin_bubble, node_end_bubble)
                paths = list(nx.all_simple_paths(graph, starting_node, sink_node))
            '''
            if len(paths) > 1:
                for path in paths[1:]:
                    differences_list = [x for x in paths[0] if x not in path]
                    begin_bubble_list = []
                    end_bubble_list = []
                    while True:
                        #begin_bubble_list.append(path[0][path[0].index(differences_list[0])-1])
                        begin_tmp = paths[0][paths[0].index(differences_list[0])-1]
                        del differences_list[0]
                        cpt = 2
                        while True:
                            if begin_tmp + cpt >= len(paths[0]):
                                break
                            if differences_list[0] != paths[0][paths[0].index(begin_tmp + cpt)]:
                                end_tmp = paths[0][paths[0].index(differences_list[0]+1)]
                                break
                            cpt += 1
                            del differences_list[0]
                        begin_bubble_list.append(begin_tmp)
                        end_bubble_list.append(end_tmp)
                        if len(differences_list) >= 1:
                            break
                    for bubble_indice in enumerate(begin_bubble_list):
                        graph = solve_bubble(graph, begin_bubble_list[bubble_indice], end_bubble_list[bubble_indice])
           '''
    return graph


def solve_entry_tips(graph, starting_nodes) :
    """prend un graphe et une liste de noeuds d’entrée et retourne graphe sans \
    chemin d’entrée indésirable"""
    node = starting_nodes[0]
    while len(list(graph.successors(node))) != 0:
        if len(list(graph.predecessors(node))) > 1:
            last_node = node
        node = list(graph.successors(node))[0]
    path_list = []
    path_size_list = []
    average_weight_list = []
    for starting_node in starting_nodes:
        for path in nx.all_simple_paths(graph, starting_node, last_node):
            path_list.append(path)
            path_size_list.append(len(path))
            average_weight_list.append(path_average_weight(graph, path))
    select_best_path(graph, path_list, path_size_list, average_weight_list, delete_entry_node = True, delete_sink_node = False)
    return graph
    
        


def solve_out_tips(graph, sink_nodes) :
    """prend un graphe et une liste de noeuds de sortie et retourne graphe sans \
    chemin de sortie indésirable"""
    node = sink_nodes[0]
    while len(list(graph.predecessors(node))) != 0:
        if len(list(graph.successors(node))) > 1:
            first_node = node
        node = list(graph.predecessors(node))[0]
    path_list = []
    path_size_list = []
    average_weight_list = []
    for sink_node in sink_nodes:
        for path in nx.all_simple_paths(graph, first_node, sink_node):
            path_list.append(path)
            path_size_list.append(len(path))
            average_weight_list.append(path_average_weight(graph, path))
    select_best_path(graph, path_list, path_size_list, average_weight_list, delete_entry_node = False, delete_sink_node = True)
    return graph


def main():
    """Main program function"""
    # Get arguments
    args = get_arguments()

    # Création du graphe de de Bruijn
    kmer_dict = build_kmer_dict(args.fastq_file, args.kmer_size)

    graph = build_graph(kmer_dict)

    #nx.draw(G, with_labels=True, font_weight='bold')
    #matplotlib.pyplot.show()

    # Parcours du graphe de de Bruijn
    starting_nodes = get_starting_nodes(graph)
    sink_nodes = get_sink_nodes(graph)
    contigs_list = get_contigs(graph, starting_nodes, sink_nodes)
    save_contigs(contigs_list, args.output_file)

    # Simplification du graphe de De Bruijn
    path_average_weight(graph, path)
    remove_paths(graph, path_list, delete_entry_node, delete_sink_node)


if __name__ == '__main__':
    main()
