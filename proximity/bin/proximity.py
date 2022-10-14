# -*- conding: utf-8 -*-
# Function: Compute network proximity
# Usage: python proximity.py parm1 parm2 parm3 parm4 parm5
# parm1: The shortest distance file for any two nodes in a PPI network.
# parm2: Drugs and their corresponding targets.
# parm3: Gene pairs of disease states.
# parm4: At least 100 proteins were included in each bin.
# parm5: The location of the output file.
# Example: python proximity.py ../examples/input/PPI.pkl ../examples/input/drug_data.txt ../examples/input/HMI_pair.txt ../examples/input/degree-bin100.txt ../examples/output/HMI

import numpy as np
import sys
import networkx as nx
import random
import pickle
import os

def read_network(file):
    with open(file, 'rb') as f:
        f = pickle.load(f)
    return f

def read_gene_pair(file):
    f= open (file)
    gene_pair=[]
    for i in f:
        gene_pair.append(i)
    gene_pair_num = len(gene_pair)
    return gene_pair,gene_pair_num

def read_drug(file):
    f= open(file)
    drug_data = dict()
    drug_drug_data_num = dict()
    for line in f:
        line=line.rstrip('\n')
        line=line.rstrip()
        drug=line.split('\t')
        drug_data[drug[0]] = drug[1:]
        drug_drug_data_num[drug[0]] = len(drug[1:])
    return drug_data,drug_drug_data_num

def read_bin(file):
    f= open(file)
    bin_node = dict()
    node_bin = dict()
    for line in f:
        line=line.rstrip('\n')
        line=line.rstrip('\t')
        line=line.split('\t')
        line = [i for i in line if i != '']
        bin_node[line[0]] = line[1:]
        for i in line[1:]:
            node_bin[i] = line[0]
    return bin_node,node_bin

def choice(node,fil1,fil2):
    node_degree = fil1[node]
    random_node = random.choice(fil2[node_degree])
    return random_node



def proximity(G,drug_data,drug_drug_data_num,gene_pair,gene_pair_num,bin_node, node_bin,output):
    with open(output+'/proximity.txt', 'w+') as file_object:
        for drug in drug_data:
            dis_drug_data = []
            tar_pair = []
            for t in list(drug_data[drug]):
                dis1 = []
                tar_pair1 = {}
                for n in gene_pair:
                    n = n.rstrip('\n')
                    pair = n.split("\t")
                    if int(pair[0]) in G[int(t)] and int(pair[1]) in G[int(t)]:
                        p1 = G[int(t)][int(pair[0])]
                        p2 = G[int(t)][int(pair[1])]
                        dis = p1 + p2
                        dis1.append(dis)
                        pa = pair[0] + '-' + pair[1]
                        tar_pair1.setdefault(dis, []).append(t)  
                        tar_pair1.setdefault(dis, []).append(pa)
                    else:
                        dis1.append(100)
                dis2 = min(dis1)
                tar_pair.append(tar_pair1[min(dis1)])
                dis_drug_data.append(dis2)

            dis_pair = []
            for n in gene_pair:
                n = n.rstrip('\n')
                pair = n.split("\t")
                dis3 = []
                for t in list(drug_data[drug]):
                    if int(t) in G[int(pair[0])] and int(t) in G[int(pair[1])]:
                        p1 = G[int(pair[0])][int(t)]
                        p2 = G[int(pair[1])][int(t)]
                        dis = p1 + p2
                        dis3.append(dis)
                    else:
                        dis3.append(100)
                dis4 = min(dis3)
                dis_pair.append(dis4)
            distance1 = (1 / drug_drug_data_num[drug]) * sum(dis_drug_data)
            distance2 = (1 / gene_pair_num) * sum(dis_pair)
            distance = distance1 + distance2
            file_object.writelines(str(drug) + '\t' + str(distance) + '\t' + str(distance1) + '\t' + str(distance2) + '\n')

def permutation(G,drug_data,drug_drug_data_num,gene_pair,gene_pair_num,bin_node, node_bin,output):
    isExists = os.path.exists(output+'/permutation')
    if not isExists:
        os.makedirs(output+'/permutation')
    for i in range(1,1001):
        print(str(i)+'\t')
        with open(output+'/permutation/out'+str(i)+'.txt','w+') as file_object:
            for drug in drug_data:
                dis_drug_data=[]
                ran_drug_data=[]
                for t in list(drug_data[drug]):
                    for k in range(1,1000):
                        r=choice(t,node_bin,bin_node)
                        if r not in ran_drug_data:
                           t=r
                           ran_drug_data.append(t)
                           break

                    dis1 = []
                    for n in gene_pair:
                        n = n.rstrip('\n')
                        pair = n.split("\t")
                        if int(pair[0]) in G[int(t)] and int(pair[1]) in G[int(t)]:
                            p1 = G[int(t)][int(pair[0])]
                            p2 = G[int(t)][int(pair[1])]
                            dis = p1 + p2
                            dis1.append(dis)
                        else:
                            dis1.append(100)
                    dis2 = min(dis1)
                    dis_drug_data.append(dis2)

                dis_pair = []
                for n in gene_pair:
                    n = n.rstrip('\n')
                    pair = n.split("\t")
                    dis3 = []
                    for t in ran_drug_data:
                        if int(t) in G[int(pair[0])] and int(t) in G[int(pair[1])]:
                            p1 = G[int(pair[0])][int(t)]
                            p2 = G[int(pair[1])][int(t)]
                            dis = p1 + p2
                            dis3.append(dis)
                        else:
                            dis3.append(100)
                    dis4 = min(dis3)
                    dis_pair.append(dis4)

                distance1 = (1 / drug_drug_data_num[drug]) * sum(dis_drug_data)
                distance2 = (1 / gene_pair_num) * sum(dis_pair)
                distance = distance1 + distance2
                file_object.write(
                    str(drug) + '\t' + str(distance) + '\t' + str(distance1) + '\t' + str(distance2) + '\n')


def combine_1000(output):
    distance1 = []
    distance2 = []
    distance = []

    with open(output+'/permutation/out1.txt') as file_object:
        lines = file_object.readlines()
    for line in lines:
        distance.append([])
        distance1.append([])
        distance2.append([])

    for i in range(1, 1001):
        print(str(i) + '\t')
        with open(output+'/permutation/out' + str(i) + '.txt') as file_object:
            lines = file_object.readlines()
        j = 0
        for line in lines:
            line = line.rstrip('\n')
            line = line.rstrip()
            pair = line.split('\t')
            distance[j].append(float(pair[1]))
            distance1[j].append(float(pair[2]))
            distance2[j].append(float(pair[3]))
            j += 1

    ave_distance = []
    for i in range(len(distance)):
        t1 = [np.average(distance[i]), np.std(distance[i])]
        t2 = [np.average(distance1[i]), np.std(distance1[i])]
        t3 = [np.average(distance2[i]), np.std(distance2[i])]
        ave_distance.append([t1, t2, t3])
    with open(output+'/permutation.txt', 'w+') as save:
        for i in range(len(ave_distance)):
            save.writelines(str(ave_distance[i][0][0]) + '\t' + str(ave_distance[i][0][1]) + '\t' + str(
                ave_distance[i][1][0]) + '\t' \
                            + str(ave_distance[i][1][1]) + '\t' + str(ave_distance[i][2][0]) + '\t' + str(
                ave_distance[i][2][1]) + '\n')



def z_score(output):
    with open(output+'/proximity.txt', encoding='utf-8') as file_obj:
        lines_pro = file_obj.readlines()
    with open(output+'/permutation.txt', encoding='utf-8') as file_obj:
        lines_per = file_obj.readlines()
    with open(output+'/z_score.txt', 'w+') as zfile:
        z_col = []
        i=0
        for line_pro in lines_pro:
            line_pro = line_pro.rstrip('\n').rstrip()
            pair_pro = line_pro.split('\t')

            line_per =lines_per[i].rstrip('\n').rstrip()
            pair_per = line_per.split('\t')

            t1=  (float(pair_pro[1]) - float(pair_per[0])) / float(pair_per[1])
            t2 = (float(pair_pro[2]) - float(pair_per[2])) / float(pair_per[3])
            t3 = (float(pair_pro[3]) - float(pair_per[4])) / float(pair_per[5])
            i = i + 1
            z_col.append([pair_pro[0],t1,t2,t3])

        for i in range(len(z_col)):
            s =z_col[i][0]+ '\t'+str(z_col[i][1]) + '\t' + str(z_col[i][2]) + '\t' + str(z_col[i][3]) + '\n'
            zfile.writelines(s)

if __name__ == '__main__':
    input1 = sys.argv[1]
    input2 = sys.argv[2]
    input3 = sys.argv[3]
    input4 = sys.argv[4]
    output = sys.argv[5]
    isExists = os.path.exists(output)
    if not isExists:
        os.makedirs(output)


    G = read_network(input1)
    drug_data,drug_drug_data_num = read_drug(input2)
    gene_pair,gene_pair_num = read_gene_pair(input3)
    bin_node, node_bin = read_bin(input4)

    proximity(G,drug_data,drug_drug_data_num,gene_pair,gene_pair_num,bin_node, node_bin,output)
    permutation(G,drug_data,drug_drug_data_num,gene_pair,gene_pair_num,bin_node, node_bin,output)
    combine_1000(output)
    z_score(output)


