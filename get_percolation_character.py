## calculate cluster, loop, break_bond, percolation new version
#!/usr/bin/python3
import hoomd
import hoomd.md
import gsd
import gsd.hoomd
import time
import numpy as np
import random 
import math
import sys
import freud
import matplotlib.pyplot as plt
import networkx as nx


sample=sys.argv[1]
n_polymer = int(sys.argv[2])
N_length = int(sys.argv[3])
f_sticker=int(sys.argv[4])
L_final = float(sys.argv[5])
eps1= float(sys.argv[6])  ## Interaction strength of sticky potential
sigmas= 1 

real_l=int((N_length-f_sticker)/(f_sticker+1))

filename1="trajectory_cnf_"+str(sample)+"_epsilon1_"+str(eps1)+"_polymer_"+str(int(n_polymer))+"_length_"+str(N_length)+"_sticker_"+str(f_sticker)+"_sigmas_"+str(sigmas)+"_boxsize_"+str(L_final)+".gsd"

traj = gsd.hoomd.open(filename1, 'rb')

bondarray  = traj[0].bonds.group
types = traj[0].particles.typeid

sticker_index=types==1


N_tot=n_polymer*N_length
re_size = n_polymer*(N_length-1)
sticker_index=sticker_index[0:N_tot]


Lx=traj[0].configuration.box[0]
Ly=traj[0].configuration.box[1]
Lz=traj[0].configuration.box[2]
box=freud.box.Box(Lx=Lx,Ly=Ly,Lz=Lz)
box1=np.asarray([Lx,Ly,Lz])
dt=0.005*(traj[1].configuration.step-traj[0].configuration.step)

def find_spanning_path(points, bonds_array, box, tol=1.5):
    
    N = len(points)
    Lx, Ly, Lz = box

    # shifts along xyz
    shifts = [N, 2*N, 3*N]
    lengths = [Lx, Ly, Lz]

    for dim in range(3):  # 0=x,1=y,2=z
        shift = shifts[dim]
        L = lengths[dim]

        G = nx.Graph()
        # 1. original bond
        for i, j in bonds_array:
            G.add_edge(i, j)
        # 2. image bond
        bonds_img = bonds_array + shift
        for i, j in bonds_img:
            G.add_edge(i, j)
        # 3. cross boundary bond (+ axis)
        for i, j in bonds_array:
            if points[i, dim] >= L/2 - tol and points[j, dim] <= -L/2 + tol:
                G.add_edge(i, j+shift)
                G.add_edge(i+shift, j)
                G.remove_edge(i, j)
                G.remove_edge(i+shift, j+shift)
            if points[j, dim] >= L/2 - tol and points[i, dim] <= -L/2 + tol:
                G.add_edge(i, j+shift)
                G.add_edge(i+shift, j)
                G.remove_edge(i, j)
                G.remove_edge(i+shift, j+shift)

        start_node = np.random.choice(list(range(N)))
        try:
        	path = nx.shortest_path(G, source=start_node, target=start_node+shift)
        	# print("find spanning path")
        	return path, dim, 1
        except nx.NetworkXNoPath:
        	# print("no path")
        	continue
    return None, None, 0


for i in range(len(traj)):
	if i<50:
		continue
	points =traj[i].particles.position[0:N_tot]
	bondarray  = traj[0].bonds.group[0:re_size]

	type1_indices = np.where(types == 1)[0]
	
	cluster = freud.cluster.Cluster()
	cluster.compute((box,points[sticker_index]),neighbors={"r_max": 0.5})
	r=cluster.cluster_idx
	b=cluster.cluster_keys
	
	keyswith2=[]
	iii=0
	inddd=0
	for j in range(len(b)):
		if len(b[j])>2:
			# print(b[i])
			inddd=1
		if len(b[j])==2:
			iii+=1
	
	cl_props = freud.cluster.ClusterProperties()
	cl_props.compute((box,points[sticker_index]), cluster.cluster_idx)
	n=cl_props.sizes  ###n[i] will return the value of the size of each cluster
	cl2index=n==2
	keys=np.asarray(b,dtype=object)
	keyswith2=np.asarray(keys[cl2index].tolist())
	if len(keyswith2)==0:
		filename2="./binding_inf_new/"+str(sample)+"_2_epsilon1_"+str(eps1)+"_1_binding_information_polymer_"+str(int(n_polymer))+"_length_"+str(N_length)+"_sticker_"+str(f_sticker)+"_boxsize_"+str(L_final)+"_.txt"
		a=np.append(i*dt,len(keyswith2)*2/(n_polymer*f_sticker))
		a=np.append(a,len(keyswith2))
		a=np.append(a,0)
		a=np.append(a,0)
		a=np.append(a,0)
		a=np.append(a,inddd)
		a=np.append(a,0)
		a=np.append(a,0)
		a=np.append(a,0)
		a=np.append(a,0)
		a=a.reshape(1,a.shape[0])
		with open(filename2,'a+') as f3:
			np.savetxt(f3,a)
		continue


	if len(keyswith2)>0:
		keyswith2=type1_indices[keyswith2]
		# print(keyswith2)

		newlis=keyswith2-keyswith2%N_length
		# print(newlis)

		n_intra=0
		for k in range(len(newlis)):
			if newlis[k][0]==newlis[k][1]:
				n_intra+=1
		n_inter=len(newlis)-n_intra
		#print(i,'_binding_information')
		#print(iii*2/(f_sticker*n_polymer),n_intra,n_inter)

	#continue

	bonds_array = np.copy(bondarray)
	bonds_array  = np.append(bondarray,keyswith2,axis=0)
	bonds_array=np.sort(bonds_array,axis=1)
	bonds_array=np.unique(bonds_array,axis=0)

	system = freud.AABBQuery.from_system((box,points))
	distances = np.linalg.norm(box.wrap(points[bonds_array[:, 1]] - points[bonds_array[:, 0]]),axis=1)
	neighbors = freud.locality.NeighborList.from_arrays(len(points),len(points),
	                bonds_array[:, 0],
	                bonds_array[:, 1],
	                distances,
	            )

	cl = freud.cluster.Cluster()
	cl1=cl.compute(system=system, neighbors=neighbors)

	cl_props = freud.cluster.ClusterProperties()
	cl_props.compute(system, cl.cluster_idx)
	n=cl_props.sizes
	b=cl1.cluster_keys
	cluster_idx = cl.cluster_idx

	bonds_array = np.column_stack((neighbors.query_point_indices, neighbors.point_indices))
	cluster_keys1 = list(range(len(b)))  # Replace lists with unique integer IDs
	
	maxcl_size=len(b[0])/N_length
	
	bonds_by_cluster = {key: [] for key in cluster_keys1}

	# Filter bonds based on cluster indices
	for bond in bonds_array:
	    p1, p2 = bond  # Get particle indices in the bond
	    cluster1 = cluster_idx[p1]
	    cluster2 = cluster_idx[p2]

	    # Add the bond to the corresponding cluster if both particles belong to the same cluster
	    if cluster1 == cluster2:
	        bonds_by_cluster[cluster1].append((p1, p2))
	# Convert bond lists to arrays for easier handling
	for key in bonds_by_cluster:
	    bonds_by_cluster[key] = np.array(bonds_by_cluster[key])

	spanning_clusters=[]
	for kkk in range(len(b)):
		if len(b[kkk])<=N_length:
			continue
		cluster_indices = b[kkk]
		cl_point=points[b[kkk]]
		
		global_to_local = {global_idx: local_idx for local_idx, global_idx in enumerate(cluster_indices)}
		local_to_global = {local_idx: global_idx for local_idx, global_idx in enumerate(cluster_indices)}
		cl_bonds_array_global = bonds_by_cluster[kkk]
		cl_bonds_array = np.array([[global_to_local[i], global_to_local[j]] 
                               for i, j in cl_bonds_array_global])
		path, dim, spanclust=find_spanning_path(cl_point, cl_bonds_array, box1, tol=1.5)
		
		spanning_clusters.append(spanclust)
	
	spanclust=sum(spanning_clusters)
	
	interbind_bond=[]
	
	nnnn=0
	for j in range(len(bonds_by_cluster)):
		ppp=0
		newlis=[]
		for k in range(len(bonds_by_cluster[j])):
			a1=bonds_by_cluster[j][k][0]-(bonds_by_cluster[j][k][0]%N_length)
			a2=bonds_by_cluster[j][k][1]-(bonds_by_cluster[j][k][1]%N_length)
			if a1!=a2:
				if ppp==0:
					newlis=np.append(a1,a2)
					newlis=newlis.reshape(1,newlis.shape[0])
					ppp=1
				else:
					c=np.append(a1,a2)
					c=c.reshape(1,c.shape[0])
					newlis=np.append(newlis,c,axis=0)
		nnnn+=len(newlis)
		interbind_bond.append(newlis)
	singlenum=[]
	cycle_number=[]
	edges_to_remove_break_cycles=[]
	for ii in range(len(interbind_bond)):
		if len(interbind_bond[ii])!=0:
			vertices = np.unique(np.array(interbind_bond[ii]).flatten())  # Numeric vertices
			edges = interbind_bond[ii]  # Numeric edges
			G = nx.MultiGraph()
			G.add_nodes_from(vertices)  # Add vertices
			G.add_edges_from(edges)     # Add edges
			
			## get singlenum
			single_nodes = [node for node, degree in G.degree() if degree == 1]
			singlenum=np.append(singlenum,len(single_nodes))

			simple_G = nx.Graph(G)
			E = simple_G.number_of_edges()
			V = simple_G.number_of_nodes()
			C = nx.number_connected_components(simple_G)
			cycle_number = np.append(cycle_number,E - V + C+(len(edges)-E))
			
			mst = nx.minimum_spanning_tree(simple_G)

			edges_to_remove = list(set(simple_G.edges) - set(mst.edges))

			edges_to_remove_break_cycles=np.append(edges_to_remove_break_cycles,len(edges_to_remove)+(len(edges)-E))

	if spanclust>=1:
		spanclust=1
	else:
		spanclust=0

	filename2="./binding_inf_new/"+str(sample)+"_2_epsilon1_"+str(eps1)+"_1_binding_information_polymer_"+str(int(n_polymer))+"_length_"+str(N_length)+"_sticker_"+str(f_sticker)+"_boxsize_"+str(L_final)+"_.txt"
	a=np.append(i*dt,len(keyswith2)*2/(n_polymer*f_sticker))
	a=np.append(a,len(keyswith2))
	a=np.append(a,np.sum(singlenum))
	a=np.append(a,np.sum(cycle_number)+n_intra)
	a=np.append(a,np.sum(edges_to_remove_break_cycles)+n_intra)
	a=np.append(a,inddd)
	a=np.append(a,n_intra)
	a=np.append(a,n_inter)
	a=np.append(a,spanclust)
	a=np.append(a,maxcl_size/n_polymer)
	a=a.reshape(1,a.shape[0])
	with open(filename2,'a+') as f3:
		np.savetxt(f3,a)

exit()

