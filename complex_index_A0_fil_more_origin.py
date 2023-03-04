# -*- coding: utf-8 -*-
"""
Created on Wed May 25 09:35:27 2022

@author: bijialin
"""

import gudhi 
from scipy.sparse import dok_matrix
from scipy.io import loadmat
import numpy as np
import networkx as nx
import itertools
import math
import scipy.linalg
import numpy.matlib 
import pandas as pd
import csv


#Sorting simplices by dimension k and sorting simplices by value in the same dimension
def get_diff_len_simplex(Simplex_list):
    #key=k value:k-simplex 
    Simplex_dict={}
    for simplex in Simplex_list:
        dimen=len(simplex)-1
        if dimen not in Simplex_dict:
            Simplex_dict[dimen]=[simplex]
        else:
            Simplex_dict[dimen].append(simplex)
    for key, value in Simplex_dict.items():
        Simplex_dict[key]=sorted(Simplex_dict[key])
    return Simplex_dict  



#i -simplex i-boundary lower adjacency
def boundary_operator_fromdict(Simplex_dict, i):
    if i in Simplex_dict:
        source_simplices = list(Simplex_dict[i])
    else:
        source_simplices=[]
    if i-1 in Simplex_dict:    
        target_simplices = list(Simplex_dict[i-1])
    else:
        target_simplices=[]

    if len(target_simplices)==0 or len(source_simplices)==0:
        S = dok_matrix((1, len(source_simplices)), dtype=np.float64)
        S[0, 0:len(source_simplices)] = 0
    else:
        source_simplices_dict = {source_simplices[j]: j for j in range(len(source_simplices))}
        target_simplices_dict = {target_simplices[i]: i for i in range(len(target_simplices))}

        S = dok_matrix((len(target_simplices), len(source_simplices)), dtype=np.float64)
        for source_simplex in source_simplices:
            for a in range(len(source_simplex)):
                target_simplex = source_simplex[:a]+source_simplex[(a+1):]
                i = target_simplices_dict[target_simplex]
                j = source_simplices_dict[source_simplex]
                S[i, j] = -1 if a % 2==1 else 1
    
    return S

#incidence_matrix：
# Step1
#k1-simplex  k2_Simplex k2_incidence matrix 
def get_simplex_incidence_matrix(k1_Simplex,k2_Simplex):
    N=len(k1_Simplex)
    M=len(k2_Simplex)
    simplex_incidence_matrix=np.zeros((N,M))
    for i in range(N):
        for j in range(M):
            if set(k1_Simplex[i]).issubset(set(k2_Simplex[j])):
                simplex_incidence_matrix[i][j]=1  
    return simplex_incidence_matrix

# Step2
#upper_adjacency_matrix,degree
def get_adjacency_matrix(incidence_matrix):
    adjacency_matrix=np.dot(incidence_matrix,incidence_matrix.T)
    diagonal=np.diag_indices(len(adjacency_matrix[0]))
    degrees_matrix=adjacency_matrix[diagonal]
    adjacency_matrix[diagonal]=0
    return degrees_matrix,adjacency_matrix

# Step2
#lower_adjacency_matrix,degree 
def get_edge_adjacency_matrix(incidence_matrix): 
    edge_adjacency_matrix=np.dot(incidence_matrix.T,incidence_matrix)
    diagonal=np.diag_indices(len(edge_adjacency_matrix[0]))
    degrees_matrix=edge_adjacency_matrix[diagonal]
    edge_adjacency_matrix[diagonal]=0
    return degrees_matrix,edge_adjacency_matrix

         
# Step3
#upper_adjacency_matrix,lower_adjacency_matrix,general adjacency matrix
def lower_uper_k_adjacency_matrix(Simplex_dict,dimen):
    upper_adjacency_matrix=None
    if Simplex_dict.__contains__(dimen+1):
        simplex_incidence_matrix=get_simplex_incidence_matrix(Simplex_dict[dimen],Simplex_dict[dimen+1])#L_k+1    k_num=m  k+_num=n m*n    
        upper_degrees_matrix,upper_adjacency_matrix=get_adjacency_matrix(simplex_incidence_matrix)#uper L_k+1*L_k+1^T uper m*m
    simplex_incidence_matrix2=get_simplex_incidence_matrix(Simplex_dict[dimen-1],Simplex_dict[dimen])#L_k    k-_num=l  k_num=m   l*m  
    lower_degrees_matrix,lower_adjacency_matrix=get_edge_adjacency_matrix(simplex_incidence_matrix2)# lower L_k^T*L_k  m*m
    simplex_adjacency_matrix=lower_adjacency_matrix.copy()#m*m
    if Simplex_dict.__contains__(dimen+1):
        simplex_adjacency_matrix[upper_adjacency_matrix==1]=0
    return upper_adjacency_matrix,lower_adjacency_matrix,simplex_adjacency_matrix

#shortest path
#and shortest distance
def shortest_path(simplex_adjacency_matrix):
    mat = simplex_adjacency_matrix
    Graph = nx.from_numpy_matrix(mat)
    p = nx.shortest_path(Graph)
    shortest_path_mat=np.matrix(np.ones(mat.shape) * np.inf)
    count=0
    nodes_dict=dict()
    for a in Graph.nodes():
        nodes_dict[a]=count
        count+=1   
    for a in Graph.nodes():
        pi=p[a]#node a 可以到达的点的路径
        for b in pi:
            i=nodes_dict[a]
            j=nodes_dict[b]
            shortest_path_mat[i, j]=len(pi[b])-1
    return p,shortest_path_mat


#shortest path 
#and shortest distance
def shortest_path_from_graph(Graph):
    p = nx.shortest_path(Graph)  
    n=nx.number_of_nodes(Graph)
    shortest_path_mat=np.matrix(np.ones((n,n)) * np.inf)
    count=0
    nodes_dict=dict()
    for a in Graph.nodes():
        nodes_dict[a]=count
        count+=1   
    for a in Graph.nodes():
        pi=p[a]#node a 可以到达的点的路径
        for b in pi:
            i=nodes_dict[a]
            j=nodes_dict[b]
            shortest_path_mat[i, j]=len(pi[b])-1
    return p,shortest_path_mat


# shortest path number      
def shortest_path_num(P):
    P_key=P.keys()
    shortest_num=0
    for key in P_key:
        shortest_num+=len(P[key].values())-1
    return shortest_num/2  


#matrixPow
def matrixPow(Matrix,n):
    if(type(Matrix)==list):
        Matrix=np.array(Matrix)
    if(n==1):
        return Matrix
    else:
        return np.matmul(Matrix,matrixPow(Matrix,n-1))
    
    
def generate_simplex_from_clique(lis):
    simple_set=[]
    simple_list=[]
    for li in lis:
        for i in range(len(li)-1):
            simple_set.extend(list(itertools.combinations(list(li),i+1)))
    simple_set=set(simple_set)
    for i in set(simple_set):
        simple_list.append(list(i))
    simple_list.extend(lis)
    simple_list=sorted(simple_list, key= lambda x: len(x))
    return simple_list

def degree_heterogeneity(G):
    degree_h=0
    G_degree=G.degree
    G_edges=G.edges
    for e in G_edges:
        degree_h+=((G_degree[e[0]])**(-1/2)-(G_degree[e[1]])**(-1/2))**2
    degree_h=degree_h/(len(G.nodes)-2*np.sqrt(len(G.nodes)-1))#sd
    return degree_h    


if __name__ == "__main__":
    
    Note2=open('name.txt',mode='r')
    file_name= eval(Note2.read())
    Note2.close()
    
    
    resid1='C'
    resid2='C'
    area='_mutation'#binding
    filtration_range=[i/10.0 for i in range(10,101)]
    for ename in range(344,345):#333 527
        print(ename)
        '''
        #generate simplices from distance matrix--------
        file=pd.read_excel('./coodinates/origin_mutation_distance_'+resid1+'_'+resid2+'/'+'6M0J_origin_'+str(ename)+'_'+resid1+'_'+resid2+area+'_distance.xlsx',index_col=0)
        file_col=file.columns
        file_index=file.index
        file_value=file.values
        '''
        

        #generate simplices from-coodinates
        file=pd.read_excel('./coodinates/origin_mutation_coodinates_'+resid1+'_'+resid2+'/'+'6M0J_'+str(ename)+'_'+resid1+'_'+resid2+'_origin'+area+'_coodinates.xlsx',index_col=0)
        file_col=file.columns
        file_index=file.index
        file_value=file.values
        #'''
        result_name='6M0J_origin_'+str(ename)+'_'+resid1+'_'+resid2+area+'_A0_index.csv'
        File1=open('./origin_mutation_result_A0_'+resid1+'_'+resid2+'/'+result_name,"w+",newline='')
        writing_File1=csv.writer(File1)
        
        
        flag=1
        for filtration_value in filtration_range:
            #generate simplices from distance matrix--------
            # rips_complex = gudhi.RipsComplex(distance_matrix=file_value, max_edge_length=filtration_value)
            #generate simplices from-coodinates------------------------------------------------------------------ 
            rips_complex = gudhi.RipsComplex( points=file_value[:,-3:], max_edge_length=filtration_value)
            
            simplex_tree = rips_complex.create_simplex_tree(max_dimension=2)#统计simplex的最高维数
            val = simplex_tree.get_filtration()     
            simplices_list=[]
            simplices = set()                
            for filtered_value in val:
                simplices.add(tuple(filtered_value[0]))
                simplices_list.append(filtered_value[0])
            simplex_list = sorted(simplices, key = lambda x: len(x))
            # simplex_list1 = sorted(simplices, key = lambda x: len(x))
            
            #generate simplices from clique--------------------------------------------------------------------
            #lapalance
            #simplex_list=generate_simplex_from_clique([[1,2],[1,3],[3,5],[4,5],[2,3,4],[7,8],[6,7],[9]])
        
            #shorted by the dimension        
            simplex_dict=get_diff_len_simplex(simplex_list)
            if len(simplex_dict.keys()) !=1:
                #I_1
                incidence_1=get_simplex_incidence_matrix(simplex_dict[0],simplex_dict[1])
                # #I_2
                # incidence_2=get_simplex_incidence_matrix(simplex_dict[1],simplex_dict[2])
                #A0--------------------------------------------------------------------
                a0_degree,a0=get_adjacency_matrix(incidence_1)
                # print("a0")
                # #A1--------------------------------------------------------------------
                # upper_adjacency_matrix,lower_adjacency_matrix,simplex_adjacency_matrix=lower_uper_k_adjacency_matrix(simplex_dict,1)
                # print("a1")
                # #A2--------------------------------------------------------------------
                # upper_adjacency_matrix2,lower_adjacency_matrix2,simplex_adjacency_matrix2=lower_uper_k_adjacency_matrix(simplex_dict,2)
                # print("a2")
                
               
                #generate graph from laplacian
                mat=a0
                graph = nx.from_numpy_matrix(mat) 
                nodes=graph.nodes
                edges=graph.edges
                
                #是基于adjacency的index---------------------------------------------------------------------   
                
                n=len(nodes)
                m=len(edges)
                # print("n",n)
                # print("m",m)
                # print ("n", nx.number_of_nodes(graph))
                # print ("m", nx.number_of_edges(graph))
                
                #edge density 2m/(n*(n-1))
                delta=2*m/(n*(n-1))
                
                #average degree
                degree_mean=2*m/n
                
                #degree_heterogeneity
                rho=degree_heterogeneity(graph)  
                
                #clustering cofficient
                c_node = nx.clustering(graph)
                C_mean=nx.average_clustering(graph)
                C_sum=np.sum(list(c_node.values()))
                C_max=np.max(list(c_node.values()))
                C_min=np.min(list(c_node.values()))
                C_sd=np.std(list(c_node.values()), ddof=1)
                
                #average_shortest_path_length "Graph is connected."
                # if nx.is_connected(graph):
                #     L=nx.average_shortest_path_length(graph)
                #     #average diameter
                #     ad=nx.diameter(graph)
                
                #betweenness_centrality
                #使用了归一化2/((N-1)*(N-2)) c_B(v) =\sum_{s,t \in V} \frac{\sigma(s, t|v)}{\sigma(s, t)}
                bc_node =nx.betweenness_centrality(graph,normalized=True)
                bc_mean=np.mean(list(bc_node.values()))
                bc_sum=np.sum(list(bc_node.values()))
                bc_max=np.max(list(bc_node.values()))
                bc_min=np.min(list(bc_node.values()))
                bc_sd=np.std(list(bc_node.values()), ddof=1)
                
                #degree_assortativity
                r=nx.degree_assortativity_coefficient(graph)
                
                # #eigenvector_centrality  
                # e_node = nx.eigenvector_centrality(graph)
                # e_mean=np.mean(list(e_node.values()))  
                
                #subgraph_centrality
                s_node = nx.subgraph_centrality(graph)
                s_mean=np.mean(list(s_node.values()))
                s_sum=np.sum(list(s_node.values()))
                s_max=np.max(list(s_node.values()))
                s_min=np.min(list(s_node.values()))
                s_sd=np.std(list(s_node.values()), ddof=1)
                # sc_node = nx.subgraph_centrality_exp(graph)
                # sc_mean=np.mean(list(sc_node.values()))
                
                pr = nx.pagerank(graph)
                pr_mean=np.mean(list(pr.values()))
                pr_sum=np.sum(list(pr.values()))
                pr_max=np.max(list(pr.values()))
                pr_min=np.min(list(pr.values()))
                pr_sd=np.std(list(pr.values()), ddof=1)
                #k1/(N-1)
                d = nx.degree_centrality(graph)
                d_mean=np.mean(list(d.values()))
                d_sum=np.sum(list(d.values()))
                d_max=np.max(list(d.values()))
                d_min=np.min(list(d.values()))
                d_sd=np.std(list(d.values()), ddof=1)
                
                cl = nx.closeness_centrality(graph)
                cl_mean=np.mean(list(cl.values()))
                cl_sum=np.sum(list(cl.values()))
                cl_max=np.max(list(cl.values()))
                cl_min=np.min(list(cl.values()))
                cl_sd=np.std(list(cl.values()), ddof=1)
            
                #average communicability  & average communicabolity angle
                Gpq = nx.communicability_exp(graph)#字典中包含字典
                Gpq_val_list=[]
                Gpq_val=0
                Cpq_val=0
                count=0
                for i in Gpq:
                    Gpq_i=Gpq[i]
                    for key, value in Gpq_i.items():
                        if key>i:
                            Gpq_val+=value
                            Gpq_val_list.append(value)
                            Cpq_val+=math.degrees(np.arccos(Gpq[i][key]/(np.sqrt(Gpq[i][i]*Gpq[key][key]))))
                Gpq_mean=2*Gpq_val/(n*(n-1))    
                Cpq_mean=2*Cpq_val/(n*(n-1)) 
                Gpq_sum=np.sum(Gpq_val_list)
                Gpq_max=np.max(Gpq_val_list)
                Gpq_min=np.min(Gpq_val_list)
                Gpq_sd=np.std(Gpq_val_list, ddof=1)   
             
                
                # nodelist = list(graph)  # ordering of nodes in matrix
                # A = nx.to_numpy_array(graph, nodelist)
                # # convert to 0-1 matrix
                # A[A != 0.0] = 1
                # A2=matrixPow(A,2)
                # expA = scipy.linalg.expm(A)
                # expA22=scipy.linalg.expm(A2/2)
                # I=np.matlib.identity(len(A))
                       
                # k=np.sqrt(math.pi)*math.log(2)#√πln(2). 
                # E=scipy.linalg.tanhm(k*A/np.sqrt(2))
                
                  
                # Z=0.5*np.matmul((np.sqrt(2*math.pi)*E+2*I),expA22)   
                # Z=np.array(Z)
                
                # Z_mean=0   
                # for i in range(len(Z)):
                #         Z_mean+=Z[i][i]           
                # Z_mean=Z_mean/len(Z)    
                
                # Z_pq=0
                # for p in range(len(Z)):
                #     Z_pq+=np.sum(Z[p,p+1:n])           
                # Z_pq=2*Z_pq/(n*(n-1)) 
                
                
                #gravity model
                p,shortest_path_mat=shortest_path_from_graph(graph)
                degree = graph.degree()
                deg={}
                for (node, val) in degree:
                    deg[node]=val     
                gm={}
                for node in graph.nodes():
                    gm[node]=0
                    for nei in list(p[node].keys()):
                        if nei != node:
                            gm[node]+= (deg[node]*deg[nei])/(len(p[node][nei])-1)
                gm_mean=np.mean(list(gm.values()))
                gm_sum=np.sum(list(gm.values()))
                gm_max=np.max(list(gm.values()))
                gm_min=np.min(list(gm.values()))
                gm_sd=np.std(list(gm.values()), ddof=1)
                            
                            
                result={}
                result['n']=n
                result['m']=m
                result['delta']=delta
                result['degree_mean']=degree_mean
                result['rho']=rho
                result['C_mean']=C_mean
                result['C_sum']=C_sum
                result['C_max']=C_max
                result['C_min']=C_min
                result['C_sd']=C_sd
                # result['L']=[L]
                result['r']=r
                result['bc_mean']=bc_mean
                result['bc_sum']=bc_sum
                result['bc_max']=bc_max
                result['bc_min']=bc_min
                result['bc_sd']=bc_sd
                # result['e_mean']=e_mean
                result['s_mean']=s_mean
                result['s_sum']=s_sum
                result['s_max']=s_max
                result['s_min']=s_min
                result['s_sd']=s_sd
                result['pr_mean']=pr_mean
                result['pr_sum']=pr_sum
                result['pr_max']=pr_max
                result['pr_min']=pr_min
                result['pr_sd']=pr_sd
                result['d_mean']=d_mean
                result['d_sum']=d_sum
                result['d_max']=d_max
                result['d_min']=d_min
                result['d_sd']=d_sd
                result['cl_mean']=cl_mean
                result['cl_sum']=cl_sum
                result['cl_max']=cl_max
                result['cl_min']=cl_min
                result['cl_sd']=cl_sd
                result['Gpq_mean']=Gpq_mean
                result['Gpq_sum']=Gpq_sum
                result['Gpq_max']=Gpq_max
                result['Gpq_min']=Gpq_min
                result['Gpq_sd']=Gpq_sd
                #result['Cpq_mean']=Cpq_mean
                # result['Z_mean']=[Z_mean]
                # result['Z_pq']=[Z_pq]
                result['gm_mean']=gm_mean
                result['gm_sum']=gm_sum
                result['gm_max']=gm_max
                result['gm_min']=gm_min
                result['gm_sd']=gm_sd
                if flag==1:
                    writing_File1.writerow(['filtration_value']+list(result.keys()))
                    flag=0
                writing_File1.writerow([filtration_value]+list(result.values()))
            if len(simplex_dict.keys()) ==1:
                result={}
                result['n']=len(simplex_dict[0])
                result['m']=0
                result['delta']=0
                result['degree_mean']=0
                result['rho']=0
                result['C_mean']=0
                result['C_sum']=0
                result['C_max']=0
                result['C_min']=0
                result['C_sd']=0
                # result['L']=[L]
                result['r']=0
                result['bc_mean']=0
                result['bc_sum']=0
                result['bc_max']=0
                result['bc_min']=0
                result['bc_sd']=0
                # result['e_mean']=e_mean
                result['s_mean']=0
                result['s_sum']=0
                result['s_max']=0
                result['s_min']=0
                result['s_sd']=0
                result['pr_mean']=0
                result['pr_sum']=0
                result['pr_max']=0
                result['pr_min']=0
                result['pr_sd']=0
                result['d_mean']=0
                result['d_sum']=0
                result['d_max']=0
                result['d_min']=0
                result['d_sd']=0
                result['cl_mean']=0
                result['cl_sum']=0
                result['cl_max']=0
                result['cl_min']=0
                result['cl_sd']=0
                result['Gpq_mean']=0
                result['Gpq_sum']=0
                result['Gpq_max']=0
                result['Gpq_min']=0
                result['Gpq_sd']=0
                #result['Cpq_mean']=Cpq_mean
                # result['Z_mean']=[Z_mean]
                # result['Z_pq']=[Z_pq]
                result['gm_mean']=0
                result['gm_sum']=0
                result['gm_max']=0
                result['gm_min']=0
                result['gm_sd']=0
                if flag==1:
                    writing_File1.writerow(['filtration_value']+list(result.keys()))
                    flag=0
                writing_File1.writerow([filtration_value]+list(result.values()))
        File1.close()
    
    
    
=======
# -*- coding: utf-8 -*-
"""
Created on Wed May 25 09:35:27 2022

@author: bijialin
"""

import gudhi 
from scipy.sparse import dok_matrix
from scipy.io import loadmat
import numpy as np
import networkx as nx
import itertools
import math
import scipy.linalg
import numpy.matlib 
import pandas as pd
import csv


#Sorting simplices by dimension k and sorting simplices by value in the same dimension
def get_diff_len_simplex(Simplex_list):
    #key=k value:k-simplex 
    Simplex_dict={}
    for simplex in Simplex_list:
        dimen=len(simplex)-1
        if dimen not in Simplex_dict:
            Simplex_dict[dimen]=[simplex]
        else:
            Simplex_dict[dimen].append(simplex)
    for key, value in Simplex_dict.items():
        Simplex_dict[key]=sorted(Simplex_dict[key])
    return Simplex_dict  



#i -simplex i-boundary lower adjacency
def boundary_operator_fromdict(Simplex_dict, i):
    if i in Simplex_dict:
        source_simplices = list(Simplex_dict[i])
    else:
        source_simplices=[]
    if i-1 in Simplex_dict:    
        target_simplices = list(Simplex_dict[i-1])
    else:
        target_simplices=[]

    if len(target_simplices)==0 or len(source_simplices)==0:
        S = dok_matrix((1, len(source_simplices)), dtype=np.float64)
        S[0, 0:len(source_simplices)] = 0
    else:
        source_simplices_dict = {source_simplices[j]: j for j in range(len(source_simplices))}
        target_simplices_dict = {target_simplices[i]: i for i in range(len(target_simplices))}

        S = dok_matrix((len(target_simplices), len(source_simplices)), dtype=np.float64)
        for source_simplex in source_simplices:
            for a in range(len(source_simplex)):
                target_simplex = source_simplex[:a]+source_simplex[(a+1):]
                i = target_simplices_dict[target_simplex]
                j = source_simplices_dict[source_simplex]
                S[i, j] = -1 if a % 2==1 else 1
    
    return S

#incidence_matrix：
# Step1
#k1-simplex  k2_Simplex k2_incidence matrix 
def get_simplex_incidence_matrix(k1_Simplex,k2_Simplex):
    N=len(k1_Simplex)
    M=len(k2_Simplex)
    simplex_incidence_matrix=np.zeros((N,M))
    for i in range(N):
        for j in range(M):
            if set(k1_Simplex[i]).issubset(set(k2_Simplex[j])):
                simplex_incidence_matrix[i][j]=1  
    return simplex_incidence_matrix

# Step2
#upper_adjacency_matrix,degree
def get_adjacency_matrix(incidence_matrix):
    adjacency_matrix=np.dot(incidence_matrix,incidence_matrix.T)
    diagonal=np.diag_indices(len(adjacency_matrix[0]))
    degrees_matrix=adjacency_matrix[diagonal]
    adjacency_matrix[diagonal]=0
    return degrees_matrix,adjacency_matrix

# Step2
#lower_adjacency_matrix,degree 
def get_edge_adjacency_matrix(incidence_matrix): 
    edge_adjacency_matrix=np.dot(incidence_matrix.T,incidence_matrix)
    diagonal=np.diag_indices(len(edge_adjacency_matrix[0]))
    degrees_matrix=edge_adjacency_matrix[diagonal]
    edge_adjacency_matrix[diagonal]=0
    return degrees_matrix,edge_adjacency_matrix

         
# Step3
#upper_adjacency_matrix,lower_adjacency_matrix,general adjacency matrix
def lower_uper_k_adjacency_matrix(Simplex_dict,dimen):
    upper_adjacency_matrix=None
    if Simplex_dict.__contains__(dimen+1):
        simplex_incidence_matrix=get_simplex_incidence_matrix(Simplex_dict[dimen],Simplex_dict[dimen+1])#L_k+1    k_num=m  k+_num=n m*n    
        upper_degrees_matrix,upper_adjacency_matrix=get_adjacency_matrix(simplex_incidence_matrix)#uper L_k+1*L_k+1^T uper m*m
    simplex_incidence_matrix2=get_simplex_incidence_matrix(Simplex_dict[dimen-1],Simplex_dict[dimen])#L_k    k-_num=l  k_num=m   l*m  
    lower_degrees_matrix,lower_adjacency_matrix=get_edge_adjacency_matrix(simplex_incidence_matrix2)# lower L_k^T*L_k  m*m
    simplex_adjacency_matrix=lower_adjacency_matrix.copy()#m*m
    if Simplex_dict.__contains__(dimen+1):
        simplex_adjacency_matrix[upper_adjacency_matrix==1]=0
    return upper_adjacency_matrix,lower_adjacency_matrix,simplex_adjacency_matrix

#shortest path
#and shortest distance
def shortest_path(simplex_adjacency_matrix):
    mat = simplex_adjacency_matrix
    Graph = nx.from_numpy_matrix(mat)
    p = nx.shortest_path(Graph)
    shortest_path_mat=np.matrix(np.ones(mat.shape) * np.inf)
    count=0
    nodes_dict=dict()
    for a in Graph.nodes():
        nodes_dict[a]=count
        count+=1   
    for a in Graph.nodes():
        pi=p[a]#node a 可以到达的点的路径
        for b in pi:
            i=nodes_dict[a]
            j=nodes_dict[b]
            shortest_path_mat[i, j]=len(pi[b])-1
    return p,shortest_path_mat


#shortest path 
#and shortest distance
def shortest_path_from_graph(Graph):
    p = nx.shortest_path(Graph)  
    n=nx.number_of_nodes(Graph)
    shortest_path_mat=np.matrix(np.ones((n,n)) * np.inf)
    count=0
    nodes_dict=dict()
    for a in Graph.nodes():
        nodes_dict[a]=count
        count+=1   
    for a in Graph.nodes():
        pi=p[a]#node a 可以到达的点的路径
        for b in pi:
            i=nodes_dict[a]
            j=nodes_dict[b]
            shortest_path_mat[i, j]=len(pi[b])-1
    return p,shortest_path_mat


# shortest path number      
def shortest_path_num(P):
    P_key=P.keys()
    shortest_num=0
    for key in P_key:
        shortest_num+=len(P[key].values())-1
    return shortest_num/2  


#matrixPow
def matrixPow(Matrix,n):
    if(type(Matrix)==list):
        Matrix=np.array(Matrix)
    if(n==1):
        return Matrix
    else:
        return np.matmul(Matrix,matrixPow(Matrix,n-1))
    
    
def generate_simplex_from_clique(lis):
    simple_set=[]
    simple_list=[]
    for li in lis:
        for i in range(len(li)-1):
            simple_set.extend(list(itertools.combinations(list(li),i+1)))
    simple_set=set(simple_set)
    for i in set(simple_set):
        simple_list.append(list(i))
    simple_list.extend(lis)
    simple_list=sorted(simple_list, key= lambda x: len(x))
    return simple_list

def degree_heterogeneity(G):
    degree_h=0
    G_degree=G.degree
    G_edges=G.edges
    for e in G_edges:
        degree_h+=((G_degree[e[0]])**(-1/2)-(G_degree[e[1]])**(-1/2))**2
    degree_h=degree_h/(len(G.nodes)-2*np.sqrt(len(G.nodes)-1))#sd
    return degree_h    


if __name__ == "__main__":
    
    Note2=open('name.txt',mode='r')
    file_name= eval(Note2.read())
    Note2.close()
    
    
    resid1='C'
    resid2='C'
    area='_mutation'#binding
    filtration_range=[i/10.0 for i in range(10,101)]
    for ename in range(344,345):#333 527
        print(ename)
        '''
        #generate simplices from distance matrix--------
        file=pd.read_excel('./coodinates/origin_mutation_distance_'+resid1+'_'+resid2+'/'+'6M0J_origin_'+str(ename)+'_'+resid1+'_'+resid2+area+'_distance.xlsx',index_col=0)
        file_col=file.columns
        file_index=file.index
        file_value=file.values
        '''
        

        #generate simplices from-coodinates
        file=pd.read_excel('./coodinates/origin_mutation_coodinates_'+resid1+'_'+resid2+'/'+'6M0J_'+str(ename)+'_'+resid1+'_'+resid2+'_origin'+area+'_coodinates.xlsx',index_col=0)
        file_col=file.columns
        file_index=file.index
        file_value=file.values
        #'''
        result_name='6M0J_origin_'+str(ename)+'_'+resid1+'_'+resid2+area+'_A0_index.csv'
        File1=open('./origin_mutation_result_A0_'+resid1+'_'+resid2+'/'+result_name,"w+",newline='')
        writing_File1=csv.writer(File1)
        
        
        flag=1
        for filtration_value in filtration_range:
            #generate simplices from distance matrix--------
            # rips_complex = gudhi.RipsComplex(distance_matrix=file_value, max_edge_length=filtration_value)
            #generate simplices from-coodinates------------------------------------------------------------------ 
            rips_complex = gudhi.RipsComplex( points=file_value[:,-3:], max_edge_length=filtration_value)
            
            simplex_tree = rips_complex.create_simplex_tree(max_dimension=2)#统计simplex的最高维数
            val = simplex_tree.get_filtration()     
            simplices_list=[]
            simplices = set()                
            for filtered_value in val:
                simplices.add(tuple(filtered_value[0]))
                simplices_list.append(filtered_value[0])
            simplex_list = sorted(simplices, key = lambda x: len(x))
            # simplex_list1 = sorted(simplices, key = lambda x: len(x))
            
            #generate simplices from clique--------------------------------------------------------------------
            #lapalance
            #simplex_list=generate_simplex_from_clique([[1,2],[1,3],[3,5],[4,5],[2,3,4],[7,8],[6,7],[9]])
        
            #shorted by the dimension        
            simplex_dict=get_diff_len_simplex(simplex_list)
            if len(simplex_dict.keys()) !=1:
                #I_1
                incidence_1=get_simplex_incidence_matrix(simplex_dict[0],simplex_dict[1])
                # #I_2
                # incidence_2=get_simplex_incidence_matrix(simplex_dict[1],simplex_dict[2])
                #A0--------------------------------------------------------------------
                a0_degree,a0=get_adjacency_matrix(incidence_1)
                # print("a0")
                # #A1--------------------------------------------------------------------
                # upper_adjacency_matrix,lower_adjacency_matrix,simplex_adjacency_matrix=lower_uper_k_adjacency_matrix(simplex_dict,1)
                # print("a1")
                # #A2--------------------------------------------------------------------
                # upper_adjacency_matrix2,lower_adjacency_matrix2,simplex_adjacency_matrix2=lower_uper_k_adjacency_matrix(simplex_dict,2)
                # print("a2")
                
               
                #generate graph from laplacian
                mat=a0
                graph = nx.from_numpy_matrix(mat) 
                nodes=graph.nodes
                edges=graph.edges
                
                #是基于adjacency的index---------------------------------------------------------------------   
                
                n=len(nodes)
                m=len(edges)
                # print("n",n)
                # print("m",m)
                # print ("n", nx.number_of_nodes(graph))
                # print ("m", nx.number_of_edges(graph))
                
                #edge density 2m/(n*(n-1))
                delta=2*m/(n*(n-1))
                
                #average degree
                degree_mean=2*m/n
                
                #degree_heterogeneity
                rho=degree_heterogeneity(graph)  
                
                #clustering cofficient
                c_node = nx.clustering(graph)
                C_mean=nx.average_clustering(graph)
                C_sum=np.sum(list(c_node.values()))
                C_max=np.max(list(c_node.values()))
                C_min=np.min(list(c_node.values()))
                C_sd=np.std(list(c_node.values()), ddof=1)
                
                #average_shortest_path_length "Graph is connected."
                # if nx.is_connected(graph):
                #     L=nx.average_shortest_path_length(graph)
                #     #average diameter
                #     ad=nx.diameter(graph)
                
                #betweenness_centrality
                #使用了归一化2/((N-1)*(N-2)) c_B(v) =\sum_{s,t \in V} \frac{\sigma(s, t|v)}{\sigma(s, t)}
                bc_node =nx.betweenness_centrality(graph,normalized=True)
                bc_mean=np.mean(list(bc_node.values()))
                bc_sum=np.sum(list(bc_node.values()))
                bc_max=np.max(list(bc_node.values()))
                bc_min=np.min(list(bc_node.values()))
                bc_sd=np.std(list(bc_node.values()), ddof=1)
                
                #degree_assortativity
                r=nx.degree_assortativity_coefficient(graph)
                
                # #eigenvector_centrality  
                # e_node = nx.eigenvector_centrality(graph)
                # e_mean=np.mean(list(e_node.values()))  
                
                #subgraph_centrality
                s_node = nx.subgraph_centrality(graph)
                s_mean=np.mean(list(s_node.values()))
                s_sum=np.sum(list(s_node.values()))
                s_max=np.max(list(s_node.values()))
                s_min=np.min(list(s_node.values()))
                s_sd=np.std(list(s_node.values()), ddof=1)
                # sc_node = nx.subgraph_centrality_exp(graph)
                # sc_mean=np.mean(list(sc_node.values()))
                
                pr = nx.pagerank(graph)
                pr_mean=np.mean(list(pr.values()))
                pr_sum=np.sum(list(pr.values()))
                pr_max=np.max(list(pr.values()))
                pr_min=np.min(list(pr.values()))
                pr_sd=np.std(list(pr.values()), ddof=1)
                #k1/(N-1)
                d = nx.degree_centrality(graph)
                d_mean=np.mean(list(d.values()))
                d_sum=np.sum(list(d.values()))
                d_max=np.max(list(d.values()))
                d_min=np.min(list(d.values()))
                d_sd=np.std(list(d.values()), ddof=1)
                
                cl = nx.closeness_centrality(graph)
                cl_mean=np.mean(list(cl.values()))
                cl_sum=np.sum(list(cl.values()))
                cl_max=np.max(list(cl.values()))
                cl_min=np.min(list(cl.values()))
                cl_sd=np.std(list(cl.values()), ddof=1)
            
                #average communicability  & average communicabolity angle
                Gpq = nx.communicability_exp(graph)#字典中包含字典
                Gpq_val_list=[]
                Gpq_val=0
                Cpq_val=0
                count=0
                for i in Gpq:
                    Gpq_i=Gpq[i]
                    for key, value in Gpq_i.items():
                        if key>i:
                            Gpq_val+=value
                            Gpq_val_list.append(value)
                            Cpq_val+=math.degrees(np.arccos(Gpq[i][key]/(np.sqrt(Gpq[i][i]*Gpq[key][key]))))
                Gpq_mean=2*Gpq_val/(n*(n-1))    
                Cpq_mean=2*Cpq_val/(n*(n-1)) 
                Gpq_sum=np.sum(Gpq_val_list)
                Gpq_max=np.max(Gpq_val_list)
                Gpq_min=np.min(Gpq_val_list)
                Gpq_sd=np.std(Gpq_val_list, ddof=1)   
             
                
                # nodelist = list(graph)  # ordering of nodes in matrix
                # A = nx.to_numpy_array(graph, nodelist)
                # # convert to 0-1 matrix
                # A[A != 0.0] = 1
                # A2=matrixPow(A,2)
                # expA = scipy.linalg.expm(A)
                # expA22=scipy.linalg.expm(A2/2)
                # I=np.matlib.identity(len(A))
                       
                # k=np.sqrt(math.pi)*math.log(2)#√πln(2). 
                # E=scipy.linalg.tanhm(k*A/np.sqrt(2))
                
                  
                # Z=0.5*np.matmul((np.sqrt(2*math.pi)*E+2*I),expA22)   
                # Z=np.array(Z)
                
                # Z_mean=0   
                # for i in range(len(Z)):
                #         Z_mean+=Z[i][i]           
                # Z_mean=Z_mean/len(Z)    
                
                # Z_pq=0
                # for p in range(len(Z)):
                #     Z_pq+=np.sum(Z[p,p+1:n])           
                # Z_pq=2*Z_pq/(n*(n-1)) 
                
                
                #gravity model
                p,shortest_path_mat=shortest_path_from_graph(graph)
                degree = graph.degree()
                deg={}
                for (node, val) in degree:
                    deg[node]=val     
                gm={}
                for node in graph.nodes():
                    gm[node]=0
                    for nei in list(p[node].keys()):
                        if nei != node:
                            gm[node]+= (deg[node]*deg[nei])/(len(p[node][nei])-1)
                gm_mean=np.mean(list(gm.values()))
                gm_sum=np.sum(list(gm.values()))
                gm_max=np.max(list(gm.values()))
                gm_min=np.min(list(gm.values()))
                gm_sd=np.std(list(gm.values()), ddof=1)
                            
                            
                result={}
                result['n']=n
                result['m']=m
                result['delta']=delta
                result['degree_mean']=degree_mean
                result['rho']=rho
                result['C_mean']=C_mean
                result['C_sum']=C_sum
                result['C_max']=C_max
                result['C_min']=C_min
                result['C_sd']=C_sd
                # result['L']=[L]
                result['r']=r
                result['bc_mean']=bc_mean
                result['bc_sum']=bc_sum
                result['bc_max']=bc_max
                result['bc_min']=bc_min
                result['bc_sd']=bc_sd
                # result['e_mean']=e_mean
                result['s_mean']=s_mean
                result['s_sum']=s_sum
                result['s_max']=s_max
                result['s_min']=s_min
                result['s_sd']=s_sd
                result['pr_mean']=pr_mean
                result['pr_sum']=pr_sum
                result['pr_max']=pr_max
                result['pr_min']=pr_min
                result['pr_sd']=pr_sd
                result['d_mean']=d_mean
                result['d_sum']=d_sum
                result['d_max']=d_max
                result['d_min']=d_min
                result['d_sd']=d_sd
                result['cl_mean']=cl_mean
                result['cl_sum']=cl_sum
                result['cl_max']=cl_max
                result['cl_min']=cl_min
                result['cl_sd']=cl_sd
                result['Gpq_mean']=Gpq_mean
                result['Gpq_sum']=Gpq_sum
                result['Gpq_max']=Gpq_max
                result['Gpq_min']=Gpq_min
                result['Gpq_sd']=Gpq_sd
                #result['Cpq_mean']=Cpq_mean
                # result['Z_mean']=[Z_mean]
                # result['Z_pq']=[Z_pq]
                result['gm_mean']=gm_mean
                result['gm_sum']=gm_sum
                result['gm_max']=gm_max
                result['gm_min']=gm_min
                result['gm_sd']=gm_sd
                if flag==1:
                    writing_File1.writerow(['filtration_value']+list(result.keys()))
                    flag=0
                writing_File1.writerow([filtration_value]+list(result.values()))
            if len(simplex_dict.keys()) ==1:
                result={}
                result['n']=len(simplex_dict[0])
                result['m']=0
                result['delta']=0
                result['degree_mean']=0
                result['rho']=0
                result['C_mean']=0
                result['C_sum']=0
                result['C_max']=0
                result['C_min']=0
                result['C_sd']=0
                # result['L']=[L]
                result['r']=0
                result['bc_mean']=0
                result['bc_sum']=0
                result['bc_max']=0
                result['bc_min']=0
                result['bc_sd']=0
                # result['e_mean']=e_mean
                result['s_mean']=0
                result['s_sum']=0
                result['s_max']=0
                result['s_min']=0
                result['s_sd']=0
                result['pr_mean']=0
                result['pr_sum']=0
                result['pr_max']=0
                result['pr_min']=0
                result['pr_sd']=0
                result['d_mean']=0
                result['d_sum']=0
                result['d_max']=0
                result['d_min']=0
                result['d_sd']=0
                result['cl_mean']=0
                result['cl_sum']=0
                result['cl_max']=0
                result['cl_min']=0
                result['cl_sd']=0
                result['Gpq_mean']=0
                result['Gpq_sum']=0
                result['Gpq_max']=0
                result['Gpq_min']=0
                result['Gpq_sd']=0
                #result['Cpq_mean']=Cpq_mean
                # result['Z_mean']=[Z_mean]
                # result['Z_pq']=[Z_pq]
                result['gm_mean']=0
                result['gm_sum']=0
                result['gm_max']=0
                result['gm_min']=0
                result['gm_sd']=0
                if flag==1:
                    writing_File1.writerow(['filtration_value']+list(result.keys()))
                    flag=0
                writing_File1.writerow([filtration_value]+list(result.values()))
        File1.close()
    
    
    
>>>>>>> 75a591fa05177761e37a2f428e0173a1ec357a48
    
