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
# import pandas
import freud
# import matplotlib as plt
import matplotlib.pyplot as plt
import networkx as nx


sample=sys.argv[1]
n_polymer = int(sys.argv[2])
N_length = int(sys.argv[3])
f_sticker=int(sys.argv[4])
L_final = float(sys.argv[5])
eps1= float(sys.argv[6])  ## Interaction strength of sticky potential
sigmas= 1 # float(sys.argv[6]) # samll radius>=(sqrt(2)-1)sigma=0.414sigma

real_l=int((N_length-f_sticker)/(f_sticker+1))


# filename1="trajectory_cnf_"+str(sample)+"_epsilon1_"+str(eps1)+"_polymer_"+str(int(n_polymer))+"_length_"+str(N_length)+"_sticker_"+str(f_sticker)+"_boxsize_"+str(L_final)+".gsd"
filename1="trajectory_cnf_"+str(sample)+"_epsilon1_"+str(eps1)+"_polymer_"+str(int(n_polymer))+"_length_"+str(N_length)+"_sticker_"+str(f_sticker)+"_sigmas_"+str(sigmas)+"_boxsize_"+str(L_final)+".gsd"


traj = gsd.hoomd.open(filename1, 'rb')

# with gsd.hoomd.open("output_frame.gsd", 'wb') as outfile:
#     outfile.append(traj[5000])

# particleid0 = traj[0].particles.typeid==0
# particleid1 = traj[0].particles.typeid==1
# particleid2 = traj[0].particles.typeid==2
# particleid3 = traj[0].particles.typeid==3
bondarray  = traj[0].bonds.group
types = traj[0].particles.typeid

# print(types)

sticker_index=types==1
# list11=np.where(types[0:65] == 1)[0]

# print(np.where(types[0:65] == 1)[0])

# print(sum(sticker_index))


# exit()

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
	# print('------')
	# if i<len(traj)-1:
	# 	# print(i)
	# 	continue
	# if i==0:
	# 	filename2=str(sample)+"_2_epsilon1_"+str(eps1)+"_1_binding_information_polymer_"+str(int(n_polymer))+"_length_"+str(N_length)+"_sticker_"+str(f_sticker)+"_boxsize_"+str(L_final)+"_.txt"
	# 	a=np.append(i,0)
	# 	a=np.append(a,0)
	# 	a=np.append(a,0)
	# 	a=a.reshape(1,a.shape[0])
	# 	with open(filename2,'a+') as f3:
	# 		np.savetxt(f3,a)
	# 	continue
	if i<50:
		continue
	# if i%50!=0:
	# 	continue
	# print(i)
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
			# print(i)
			# print(type1_indices[b[j]])
			# p10=points[type1_indices[b[j]][0]]
			# p20=points[type1_indices[b[j]][1]]
			# p30=points[type1_indices[b[j]][2]]
			# print(box.compute_distances(p10,p20),box.compute_distances(p10,p30),box.compute_distances(p20,p30))
			# for oo in range(len(type1_indices[b[j]])):
			# bb1=bondarray[np.where(bondarray==type1_indices[b[j]][0])[0][0]][0]
			# bb2=bondarray[np.where(bondarray==type1_indices[b[j]][1])[0][0]][0]
			# bb3=bondarray[np.where(bondarray==type1_indices[b[j]][2])[0][0]][0]
			# print(bb1,bb2,bb3)
			# p1=points[bb1]
			# p2=points[bb2]
			# p3=points[bb3]
			# print(box.compute_distances(p1,p2),box.compute_distances(p1,p3),box.compute_distances(p2,p3))
			# print(box.compute_distances(p1,p10),box.compute_distances(p2,p20),box.compute_distances(p3,p30))
		if len(b[j])==2:
			iii+=1
		# 	newb=type1_indices[b[i]]%N_length
		# 	print(newb)
			# keyswith2.append(type1_indices[b[i]])
	
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

	## old spanning cluster judgement
	# points_x=np.append(points,points-[Lx,0,0],axis=0)
	# bonds_array_x = np.append(bonds_array,np.copy(bonds_array+[len(points),len(points)]),axis=0)
	
	# points_y=np.append(points,points-[0,Ly,0],axis=0)
	# bonds_array_y = np.append(bonds_array,np.copy(bonds_array+[len(points),len(points)]),axis=0)
	
	# points_z=np.append(points,points-[0,0,Lz],axis=0)
	# bonds_array_z = np.append(bonds_array,np.copy(bonds_array+[len(points),len(points)]),axis=0)
	
	# for j in range(len(points_x)):
	# 	if points_x[j,0]>=Lx/2.0-1.5:
	# 		index1=np.where(bonds_array_x==j)
	# 		n1=np.shape(index1)[1]
	# 		for l in range(n1):
	# 			if index1[1][l]==0:
	# 				pair_p=bonds_array_x[index1[0][l]][1]
	# 				if points_x[pair_p][0]<=(-Lx/2.0+1.5):
	# 					bonds_array_x[index1[0][l]][1]=pair_p+len(points)
	# 					index2=np.where(bonds_array_x==bonds_array_x[index1[0][l]][0]+len(points))
	# 					n2=np.shape(index2)[1]
	# 					for mm in range(n2):
	# 						if index2[1][mm]==0 and bonds_array_x[index2[0][mm]][1]==pair_p+len(points):
	# 							bonds_array_x[index2[0][mm]][1]=pair_p
	# 			else:
	# 				pair_p=bonds_array_x[index1[0][l]][0]
	# 				if points_x[pair_p][0]<=(-Lx/2.0+1.5):
	# 					bonds_array_x[index1[0][l]][0]=pair_p+len(points)
	# 					index2=np.where(bonds_array_x==bonds_array_x[index1[0][l]][1]+len(points))
	# 					n2=np.shape(index2)[1]
	# 					for mm in range(n2):
	# 						if index2[1][mm]==1 and bonds_array_x[index2[0][mm]][0]==pair_p+len(points):
	# 							bonds_array_x[index2[0][mm]][0]=pair_p
		
	# 	if points_y[j,1]>=Ly/2.0-1.5:
	# 		index1=np.where(bonds_array_y==j)
	# 		n1=np.shape(index1)[1]
	# 		for l in range(n1):
	# 			if index1[1][l]==0:
	# 				pair_p=bonds_array_y[index1[0][l]][1]
	# 				if points_y[pair_p][1]<=(-Ly/2.0+1.5):
	# 					bonds_array_y[index1[0][l]][1]=pair_p+len(points)
	# 					index2=np.where(bonds_array_y==bonds_array_y[index1[0][l]][0]+len(points))
	# 					n2=np.shape(index2)[1]
	# 					for mm in range(n2):
	# 						if index2[1][mm]==0 and bonds_array_y[index2[0][mm]][1]==pair_p+len(points):
	# 							bonds_array_y[index2[0][mm]][1]=pair_p
	# 			else:
	# 				pair_p=bonds_array_y[index1[0][l]][0]
	# 				if points_y[pair_p][1]<=(-Ly/2.0+1.5):
	# 					bonds_array_y[index1[0][l]][0]=pair_p+len(points)
	# 					index2=np.where(bonds_array_y==bonds_array_y[index1[0][l]][1]+len(points))
	# 					n2=np.shape(index2)[1]
	# 					for mm in range(n2):
	# 						if index2[1][mm]==1 and bonds_array_y[index2[0][mm]][0]==pair_p+len(points):
	# 							bonds_array_y[index2[0][mm]][0]=pair_p
		
	# 	if points_z[j,2]>=Lz/2.0-1.5:
	# 		index1=np.where(bonds_array_z==j)
	# 		n1=np.shape(index1)[1]
	# 		for l in range(n1):
	# 			if index1[1][l]==0:
	# 				pair_p=bonds_array_z[index1[0][l]][1]
	# 				if points_z[pair_p][2]<=(-Lz/2.0+1.5):
	# 					bonds_array_z[index1[0][l]][1]=pair_p+len(points)
	# 					index2=np.where(bonds_array_z==bonds_array_z[index1[0][l]][0]+len(points))
	# 					n2=np.shape(index2)[1]
	# 					for mm in range(n2):
	# 						if index2[1][mm]==0 and bonds_array_z[index2[0][mm]][1]==pair_p+len(points):
	# 							bonds_array_z[index2[0][mm]][1]=pair_p
	# 			else:
	# 				pair_p=bonds_array_z[index1[0][l]][0]
	# 				if points_z[pair_p][2]<=(-Lz/2.0+1.5):
	# 					bonds_array_z[index1[0][l]][0]=pair_p+len(points)
	# 					index2=np.where(bonds_array_z==bonds_array_z[index1[0][l]][1]+len(points))
	# 					n2=np.shape(index2)[1]
	# 					for mm in range(n2):
	# 						if index2[1][mm]==1 and bonds_array_z[index2[0][mm]][0]==pair_p+len(points):
	# 							bonds_array_z[index2[0][mm]][0]=pair_p
	
	# spanclust=0
	# spanning_index=[]

	# for k in range(len(b)):
	# 	# print(k)
	# 	partids=np.asarray(b[k])
	# 	if len(b[k])>65:
	# 		# print(len(b[k]))
	# 		# n_cl+=1
	# 		pos = points[partids]
	# 		check=0
	# 		for j in range(3):
	# 			fstcheck=pos[:,j]<=-box1[j]/2+2.0
	# 			sndcheck=pos[:,j]>=box1[j]/2-2.0
	# 			trdcheck=(pos[:,j]>=-2.0)&(pos[:,j]<=2.0)
	# 			forcheck=(pos[:,j]>=-box1[j]/4-2.0)&(pos[:,j]<=-box1[j]/4+2.0)
	# 			fifcheck=(pos[:,j]>=box1[j]/4-2.0)&(pos[:,j]<=box1[j]/4+2.0)
	# 			sum1=np.sum(fstcheck)
	# 			sum2=np.sum(sndcheck)
	# 			sum3=np.sum(trdcheck)
	# 			sum4=np.sum(forcheck)
	# 			sum5=np.sum(fifcheck)
	# 			check=np.logical_and(np.logical_and(np.logical_and(np.logical_and(sum1,sum2),sum3),sum4),sum5)
				
	# 			if check==1:
	# 				b_b=np.asarray(b[k])
	# 				b_2=np.append(b_b,b_b+len(points),axis=0)
	# 				b_2=np.sort(b_2)
	# 				# print(b_2)
	# 				# print(len(b_2))

	# 				# exit()
	# 				if j==0:
	# 					# print('------candidate--------')
	# 					pos_x=np.append(pos,pos-[Lx,0,0],axis=0)
	# 					# print(len(pos_x))
	# 					front=b_2[0]
	# 					queue=np.asarray(front)
	# 					father=np.asarray([[front,front]])
	# 					while (b_2[0]+len(points)) not in queue:
	# 						index1=np.where(bonds_array_x==front)
	# 						n1=np.shape(index1)[1]
	# 						for l in range(n1):
	# 							if index1[1][l]==0:
	# 								pair_p=bonds_array_x[index1[0][l]][1]
	# 								if pair_p not in father[:,0]:
	# 									queue=np.append(queue,pair_p)
	# 									father=np.append(father,[[pair_p,front]],axis=0)
	# 							else:
	# 								pair_p=bonds_array_x[index1[0][l]][0]
	# 								if pair_p not in father[:,0]:
	# 									queue=np.append(queue,pair_p)
	# 									father=np.append(father,[[pair_p,front]],axis=0)
	# 						queue=np.delete(queue,0)
	# 						# print(queue)
	# 						if len(father[:,0])==len(b_2) or len(queue)==0:
	# 							check=0
	# 							break
	# 						front=queue[0]
	# 				elif j==1:
	# 					# print('------candidate--------')
	# 					pos_y=np.append(pos,pos-[0,Ly,0],axis=0)
	# 					# print(len(pos_y))
	# 					front=b_2[0]
	# 					queue=np.asarray(front)
	# 					father=np.asarray([[front,front]])
	# 					while (b_2[0]+len(points)) not in queue:
	# 						index1=np.where(bonds_array_y==front)
	# 						n1=np.shape(index1)[1]
	# 						for l in range(n1):
	# 							if index1[1][l]==0:
	# 								pair_p=bonds_array_y[index1[0][l]][1]
	# 								if pair_p not in father[:,0]:
	# 									queue=np.append(queue,pair_p)
	# 									father=np.append(father,[[pair_p,front]],axis=0)
	# 							else:
	# 								pair_p=bonds_array_y[index1[0][l]][0]
	# 								if pair_p not in father[:,0]:
	# 									queue=np.append(queue,pair_p)
	# 									father=np.append(father,[[pair_p,front]],axis=0)
	# 						queue=np.delete(queue,0)
	# 						# print(queue)
	# 						if len(queue)==0:
	# 							check=0
	# 							break
	# 						front=queue[0]
							
	# 				else:
	# 					# print('------candidate--------')
	# 					pos_z=np.append(pos,pos-[0,0,Lz],axis=0)
	# 					# print(len(pos_z))
	# 					front=b_2[0]
	# 					queue=np.asarray(front)
	# 					father=np.asarray([[front,front]])
	# 					while (b_2[0]+len(points)) not in queue:
	# 						index1=np.where(bonds_array_z==front)
	# 						n1=np.shape(index1)[1]
	# 						for l in range(n1):
	# 							if index1[1][l]==0:
	# 								pair_p=bonds_array_z[index1[0][l]][1]
	# 								if pair_p not in father[:,0]:
	# 									queue=np.append(queue,pair_p)
	# 									father=np.append(father,[[pair_p,front]],axis=0)
	# 							else:
	# 								pair_p=bonds_array_z[index1[0][l]][0]
	# 								if pair_p not in father[:,0]:
	# 									queue=np.append(queue,pair_p)
	# 									father=np.append(father,[[pair_p,front]],axis=0)
	# 						queue=np.delete(queue,0)
	# 						# print(queue)
	# 						if len(father[:,0])==len(b_2) or len(queue)==0:
	# 							check=0
	# 							break
	# 						front=queue[0]
							
	# 				# print(b_2)
	# 				# print(len(b[k]))
	# 				# print(np.sort(father[:,0]))
	# 				# print(len(b_2))
	# 				# print(len(np.sort(father[:,0])))
					
	# 				# print(np.sort(father[:,0])==b_2)
	# 				# exit()
	# 				if check==0:
	# 					# print('-----------no path------------')
	# 					continue
	# 				else:
	# 					spanning_index=np.append(spanning_index,k)
	# 					# if (n[k]/N_tot)>=mmm:
	# 					# 	mmm=(n[k]/N_tot)
	# 					spanclust+=1
	# 					break
	# #print(spanclust,spanning_index)

	bonds_array = np.column_stack((neighbors.query_point_indices, neighbors.point_indices))
	cluster_keys1 = list(range(len(b)))  # Replace lists with unique integer IDs
	# print(n[int(spanning_index)],len(b[0]))

	# print(cluster_keys1)
	maxcl_size=len(b[0])/N_length
	
	# Create a dictionary to store bonds for each cluster
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
		# print(len(b[kkk]))
		
		global_to_local = {global_idx: local_idx for local_idx, global_idx in enumerate(cluster_indices)}
		local_to_global = {local_idx: global_idx for local_idx, global_idx in enumerate(cluster_indices)}
		cl_bonds_array_global = bonds_by_cluster[kkk]
		cl_bonds_array = np.array([[global_to_local[i], global_to_local[j]] 
                               for i, j in cl_bonds_array_global])
		# spanning_clusters=find_spanning_clusters(cl_point, cl_bonds_array, box1, tol=1.5)
		path, dim, spanclust=find_spanning_path(cl_point, cl_bonds_array, box1, tol=1.5)
		# if spanclust==1:
		# 	print(kkk,len(b[kkk]),dim)
		# 	print(np.array(path))
		# 	patheee=np.array(path)%len(b[kkk])
		# 	cl_bonds_array = np.array([[local_to_global[i]] 
        #                        for i in patheee])
		# 	print((np.array(path)%len(b[kkk]))%N_length)
		# 	exit(0)
		spanning_clusters.append(spanclust)
	# print(spanning_clusters)
	# exit()
	spanclust=sum(spanning_clusters)
	# print(spanning_clusters)
	# exit()

	interbind_bond=[]
	# trimer_list_in_spannning=[]
	# RNA_list_in_spannning=[]
	# totalbondnum=[]

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
			# if ppp==0:
			# 	newlis=np.append(a1,a2)
			# 	newlis=newlis.reshape(1,newlis.shape[0])
			# 	ppp=1
			# else:
			# 	c=np.append(a1,a2)
			# 	c=c.reshape(1,c.shape[0])
			# 	newlis=np.append(newlis,c,axis=0)
		# newlis=bonds_by_cluster[j]-bonds_by_cluster[j]%N_length
		nnnn+=len(newlis)
		interbind_bond.append(newlis)
	singlenum=[]
	# loopnum=0
	# loopbond=0
	# cyclebond=0
	cycle_number=[]
	edges_to_remove_break_cycles=[]
	# bridgenum=0
	# bondnum=0
	for ii in range(len(interbind_bond)):
		if len(interbind_bond[ii])!=0:
			# bondnum+=len(interbind_bond[ii])
			vertices = np.unique(np.array(interbind_bond[ii]).flatten())  # Numeric vertices
			edges = interbind_bond[ii]  # Numeric edges
			# print(len(edges)-len(vertices)+1)
			G = nx.MultiGraph()
			G.add_nodes_from(vertices)  # Add vertices
			G.add_edges_from(edges)     # Add edges
			# print(edges)
			# print('circle number:')
			# print(G.number_of_edges()-G.number_of_nodes()+1)

			## get singlenum
			single_nodes = [node for node, degree in G.degree() if degree == 1]
			singlenum=np.append(singlenum,len(single_nodes))

			# ## get loopnum
			# def get_edge_multiplicity(graph):
			#     # Use a dictionary to store edge multiplicity
			#     edge_multiplicity = {}
			#     for u, v in graph.edges(keys=False):
			#         # Edges in MultiGraph are treated as undirected (u, v) and (v, u) are the same
			#         edge = tuple(sorted((u, v)))  # Ensure consistent order
			#         edge_multiplicity[edge] = edge_multiplicity.get(edge, 0) + 1
			#     return edge_multiplicity

			# # Get edge multiplicity
			# edge_multiplicities = get_edge_multiplicity(G)

			# count_multiplicity_2 = sum(1 for edge, multiplicity in edge_multiplicities.items() if multiplicity == 2)
			# count_multiplicity_3 = sum(1 for edge, multiplicity in edge_multiplicities.items() if multiplicity == 3)
			# loopnum+=count_multiplicity_2+2*count_multiplicity_3
			# loopbond+=2*count_multiplicity_2+3*count_multiplicity_3
			# print('loop number:')
			# print(count_multiplicity_2+2*count_multiplicity_3)
			# print('real circle number:')
			# print(G.number_of_edges()-G.number_of_nodes()+1-(count_multiplicity_2+2*count_multiplicity_3))
			
			simple_G = nx.Graph(G)
			E = simple_G.number_of_edges()
			V = simple_G.number_of_nodes()
			C = nx.number_connected_components(simple_G)
			cycle_number = np.append(cycle_number,E - V + C+(len(edges)-E))
			# print(cycle_number)

			mst = nx.minimum_spanning_tree(simple_G)

			edges_to_remove = list(set(simple_G.edges) - set(mst.edges))

			edges_to_remove_break_cycles=np.append(edges_to_remove_break_cycles,len(edges_to_remove)+(len(edges)-E))

			# cycles = nx.cycle_basis(simple_G)
			# cycle_edges = set()  # Use a set to avoid duplicates
			# for cycle in cycles:
			#     edges_in_cycle = [(cycle[i], cycle[(i + 1) % len(cycle)]) for i in range(len(cycle))]
			#     cycle_edges.update(tuple(sorted(edge)) for edge in edges_in_cycle)  # Sort to avoid direction ambiguity
			# cycle_edges_matrix = np.array(list(cycle_edges))
			# cyclebond+=len(cycle_edges_matrix)
			
			# print('loop number:')
			# print(count_multiplicity_2+2*count_multiplicity_3)
			# print('real circle number:')
			# print(G.number_of_edges()-G.number_of_nodes()+1-(count_multiplicity_2+2*count_multiplicity_3))
			# print(edge_multiplicities)
			# print(G.degree())

			# ## plot graph
			# single_nodes = [node for node, degree in simple_G.degree() if degree == 1 and node<trimer_tot]
			# simple_G.remove_nodes_from(single_nodes)
			# node_colors = ["red" if v >= trimer_tot else "blue" for v in simple_G.nodes]
			# pos = nx.spring_layout(simple_G)  # Use spring layout for visualization
			# # colors = ["red", "green", "blue", "orange", "purple"]
			# for i, cycle in enumerate(cycles):
			# 	# Create edges for the cycle
			# 	plt.figure(figsize=(12, 10))  # Adjust figure size if needed
			# 	nx.draw(simple_G, pos, with_labels=True, node_color=node_colors, edge_color="gray",node_size=2,font_size=1)
			# 	cycle_edges = [(cycle[j], cycle[(j + 1) % len(cycle)]) for j in range(len(cycle))]
			# 	print(cycle_edges)
			# 	nx.draw_networkx_edges(simple_G, pos, edgelist=cycle_edges, edge_color="red", width=2)
			# 	plt.title("Graph with Two Types of Particles")
			# 	plt.show()
	#print(i*dt)
	#print(singlenum)
	#print(edges_to_remove_break_cycles)

	#exit()
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

