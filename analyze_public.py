import os.path
import numpy as np
from operator import itemgetter
from ete2 import Tree
import math


# path and filenames

#combined maf file
path_to_maf_file="WRITE THE PATH HERE"
maf_filename="combined_maf.txt"

#annotation file
path_annotation_file="WRITE THE PATH HERE"
annotation_filename="tcga.pancanAtlas_sample_type_annotation.txt"

#2018 seg file
path_to_2018_seg_file="WRITE THE PATH HERE"
seg_filename_2018="pan_can_atlas_mc3_caller_2p_data_cna_hg19.seg"

# PhyloWGS result
path_to_phylowgs_results="WRITE THE PATH HERE"

#output file
path_output="WRITE THE PATH HERE"
output_filename="Table.txt"


# some parameters
nb_best_trees=50
seed=[1,2,3,4,5,6,7,8,9,10] # list of trials


cancer=["acc_2015","brca_2015","chol","crc_2015","esca","hnsc","kirc","laml","lihc_2015","lusc","ov","pcpg","skcm_2015","tgct","thym","ucs","blca_2015","cesc_2015","dlbc","gbm","kich_2015","kirp_2015","lgg","luad_2015","meso","paad","prad_2015","sarc","stad_2015","thca_2015","ucec_2015","uvm_2015"]


all_sample=[]
all_sample_short=[]
all_cancer_type=[]
all_subtype=[]



############################################################
# all list of samples with cancer type and subtype

input_fn=path_annotation_file+annotation_filename
input_f=open(input_fn,"r")
a=input_f.readline()

#Assuming the annotation file has the following header
#sample_id       patient_id      type    subtype class

for line in input_f:
    lline=line.split("\t")
    s_id=lline[0]
    p_id=lline[1]

    cancer_type=lline[2].lower()
    cancer_subtype=lline[3]
    
        
    if cancer_type=="coad" or cancer_type=="read": #aggreate COAD and READ into a single CRC cancer_type
        cancer_type="crc"
            
    all_sample.append(s_id)
    all_sample_short.append(p_id)
    all_cancer_type.append(cancer_type)
    all_subtype.append(cancer_subtype)

input_f.close()

############################################################
# subset of samples analyzed by phylowgs
subset_sample=[]
subset_sample_short=[]
subset_cancer_type=[]
subset_subtype=[]

subset_fn="WRITE THE PATH HERE/list_sample_done.txt" #list of samples done
subset_f=open(subset_fn,"r")
for line in subset_f:
    lline=line.split()
    s=lline[0]
    subset_sample.append(s)
    idx_s=all_sample.index(s)
    subset_sample_short.append(all_sample_short[idx_s])
    subset_cancer_type.append(all_cancer_type[idx_s])
    sub=all_subtype[idx_s]
    subset_subtype.append(sub)
    
subset_f.close()



############################################################
# calculate number of segments
subset_nb_seg=np.zeros(len(subset_sample))

input_fn=path_to_2018_seg_file+seg_filename_2018
input_f=open(input_fn,"r")
a=input_f.readline()

for line in input_f:
    lline=line.split()
    s=lline[0]
    ss=s[:15] #IDs in the seg file have long name
    seg_mean=float(lline[5])
    if ss in subset_sample:
        idx_ss=subset_sample.index(ss)
        if abs(seg_mean) >0.3:
            subset_nb_seg[idx_ss]=subset_nb_seg[idx_ss]+1
input_f.close()


############################################################
# calculate phylogenetic properties
# requires ete2 --> from ete2 import Tree

def Create_tree_from_file(tree_fn):
    tree_f=open(tree_fn,"r")
    Tumor=Tree()
    N0=Tumor.get_tree_root()
    N0.name="N0"
    count=0
    for line in tree_f:
        lline=line.split()
        if count >0:
            parent=lline[0]
            n_p="N"+parent
            child=lline[1]
            n_c="N"+child

            node_p=Tumor.search_nodes(name=n_p)[0]
            node_c=node_p.add_child(name=n_c)
        count=1
    tree_f.close()
    return Tumor

subset_nb_nodes=np.zeros(len(subset_sample))
subset_tree_score=np.zeros(len(subset_sample))

#to store samples with nan
subset_blacklist_sample_phylo=[]

# 10 trials for each sample (variable seed)
# 50 best trees are ketp for each trial: top_k_trees1, top_k_trees2, ....

#with the following format, output from PloWGS
#-405.573595465
#0 1
#1 2 
#1 3 

#meaning:
#likelyhood of the tree
#clone 1 connected to clone 0
#clone 2 connected to clone 1
#clone 3 connected to clone 1


for i_s,sample in enumerate(subset_sample):
    cancer=subset_cancer_type[i_s]
    print cancer,sample
    list_of_llh=[]
    
    for i in seed:
        count_tree=0
        for ii in range(0,nb_best_trees,1):
            llh_fn="%s/%s/%s_%s/top_k_trees%d"%(path_to_phylowgs_results,cancer,sample,i,ii)
            if os.path.isfile(llh_fn):
                llh_f=open(llh_fn,"r")
                for line in llh_f:
                    lline=line.split()
                    if float(lline[0]) < 0 :# we get the llh
                        temp=[float(lline[0]),count_tree,i]
                        count_tree=count_tree+1
                        list_of_llh.append(temp)
                llh_f.close()
    list_of_llh.sort(key=itemgetter(0),reverse=True)

    likelyhood=np.zeros(nb_best_trees)

    if len(list_of_llh)>=nb_best_trees:
        for i in range (0,nb_best_trees,1):
            likelyhood[i]=float(list_of_llh[i][0])
    else:
        # in case all phylogenies get a nan score
        subset_blacklist_sample_phylo.append(sample)
        for i in range (0,nb_best_trees,1):
            likelyhood[i]=-1.0
            
    minllh=min(likelyhood)
    maxllh=max(likelyhood)
    
    if minllh==maxllh:
        weight=np.ones(nb_best_trees)
    else:
        weight=np.zeros(nb_best_trees)
        for i in range (0,nb_best_trees,1):
            weight[i]=(likelyhood[i]-minllh)/(maxllh-minllh)

    
    score=[]
    trunk_length=[]
    nnnode=[]
    lleaves=[]

    if len(list_of_llh)>=nb_best_trees:
        for i0 in range (0,nb_best_trees,1):
            llh=list_of_llh[i0][0]
            tree=list_of_llh[i0][1]
            seeed=list_of_llh[i0][2]
            top_trees_fn="%s/%s/%s_%s/top_k_trees%d"%(path_to_phylowgs_results,cancer,sample,seeed,tree)

            Tumor=Create_tree_from_file(top_trees_fn)

            N0=Tumor.get_tree_root()
            N1=Tumor.search_nodes(name="N1")[0]
            avdistance=0
            nbleaves=0
            nbnode=0

            flag_trunk=0
            trunk=N0

            for node in Tumor.traverse(strategy='levelorder'):
                nbnode=nbnode+1
                if flag_trunk==0 and len(node.children)>1:
                    trunk=node
                    flag_trunk=1

                if node.is_leaf():
                    nbleaves=nbleaves+1
                    avdistance=avdistance+(Tumor.get_distance(N1,node))

            if nbnode==2:
                scor=1
            else:
                scor=1.0*avdistance/(1.0*nbleaves*(nbnode-2))

            if trunk==N0:
                trunk_length.append(nbnode-1)
            else:
                trunk_length.append(Tumor.get_distance(N1,trunk))

            score.append(1.-scor)
            nnnode.append(nbnode-1)
            lleaves.append(nbleaves)

        final_score=np.average(score,weights=weight)
        final_node=np.average(nnnode,weights=weight)
    else:
        final_score=0
        final_node=0
        
    subset_nb_nodes[i_s]=final_node
    subset_tree_score[i_s]=final_score


############################################################
# calculate nb of mut and reads
subset_nb_mut=np.zeros(len(subset_sample))
subset_nb_reads=np.zeros(len(subset_sample))

maf_fn=path_to_maf_file+combined_maf.txt
maf_f=open(maf_fn,"r")

#with the following header
#Tymor_type	Tumor_Sample	Hugo_Symbol	Chromosome	Start_Position	Variant_Classification	Protein_Change	t_ref	t_alt

a=maf_f.readline()

for line in maf_f:
    lline=line.split()
    s=lline[1]
    if s in subset_sample:
        idx_s=subset_sample.index(s)
        
        subset_nb_mut[idx_s]+=1.
        
        t_ref=int(lline[7])
        t_alt=int(lline[8])
        r=1.*(t_ref+t_alt)
        subset_nb_reads[idx_s]+=r
        
maf_f.close()

for i in range(len(subset_nb_reads)):
    subset_nb_reads[i]/=subset_nbmut[i] #mean number of reads per mutation




    
############################################################
# calculate number of mutations in first clone

subset_nb_mut_first=np.zeros(len(subset_sample))
subset_nb_mut_rest=np.zeros(len(subset_sample))

for i_s,sample in enumerate(subset_sample):
    cancer=subset_cancer_type[i_s]
    list_of_llh=[]
    for i in seed:
        for ii in range(0,nb_best_trees,1):
            llh_fn="%s/%s/%s_%s/top_k_trees%d"%(path_to_phylowgs_results,cancer,sample,i,ii)
            if os.path.isfile(llh_fn):
                llh_f=open(llh_fn,"r")
                for line in llh_f:
                    lline=line.split()
                    llh=float(lline[0])
                    if math.isnan(llh):
                        llh=-1000000000
                    if llh < 0 :# we get the llh
                        list_of_llh.append([llh])

                llh_f.close()


# here we use the output top_k_trees from PhyloWGS with muational content of each clone for all the phylogenies

#typical format
#id      phi     nChildren       nGenes  genes
#0       1.000000        1       0       
#1       0.738000        1       7     s1     s9     s3     s4     s5     s6     s13     
#2       0.224000        0       6      s7   s8     s2     s10      s11     s12    
    
    count_tree=0
    for i in seed:
        toptree_fn="%s/%s/%s_%s/top_k_trees"%(path_to_phylowgs_results,cancer,sample,i)
        if os.path.isfile(toptree_fn):
            toptree_f=open(toptree_fn,"r")
            for line in toptree_f:
                lline=line.split()
                if lline[0]=="1":
                    nba=int(lline[3])
                    my_lst_str = ''.join(map(str, lline))
                    nbm=my_lst_str.count("s")
                    list_of_llh[count_tree].append(nbm)
                    count_tree+=1
            toptree_f.close()

    list_of_llh.sort(key=itemgetter(0),reverse=True)

    likelyhood=np.zeros(nb_best_trees)
    nb_mut_first=np.zeros(nb_best_trees)
    nb_mut_rest=np.zeros(nb_best_trees)
    
    for i in range (0,nb_best_trees,1):
        likelyhood[i]=float(list_of_llh[i][0])
        nb_mut_first[i]=int(list_of_llh[i][1])
        nb_mut_rest[i]=subset_nb_mut[i_s]-int(list_of_llh[i][1])
        
    minllh=min(likelyhood)
    maxllh=max(likelyhood)

    if minllh==maxllh:
        weight=np.ones(nb_best_trees)
    else:
        weight=np.zeros(nb_best_trees)
        for i in range (0,nb_best_trees,1):
            weight[i]=(likelyhood[i]-minllh)/(maxllh-minllh)
    
    mean_nb_mut_first=np.average(nb_mut_first,weights=weight)
    mean_nb_mut_rest=np.average(nb_mut_rest,weights=weight)

    subset_nb_mut_first[i_s]=mean_nb_mut_first
    subset_nb_mut_rest[i_s]=mean_nb_mut_rest


############################################################
# Output table
# sample_name  cancer_type     cancer_subtype   nb_mut    nb_seg   nb_clones  tree_score  mean_first   mean_rest

output_fn=path_output+output_filename
output_f=open(output_fn,"w")
output="sample_name\tcancer_type\tcancer_subtype\tnb_reads\tnb_mut\tnb_seg\tnb_clones\ttree_score\tmean_first\tmean_rest\n"


for i_s,sn in enumerate (subset_sample):
    if sn not in subset_blacklist_sample_phylo:
        cancertype=subset_cancer_type[i_s]
        cancersubtype=subset_subtype[i_s]
        nbread=str(subset_nb_reads[i_s])
        nbmut=str(subset_nb_mut[i_s])
        nbseg=str(subset_nb_seg[i_s])
        nbnode=str(subset_nb_nodes[i_s])
        treescore=str(subset_tree_score[i_s])
        nbmutfirst=str(subset_nb_mut_first[i_s])
        nbmutrest=str(subset_nb_mut_rest[i_s])
        
        ss=sn+"\t"+cancertype+"\t"+cancersubtype+"\t"+nbread+"\t"+nbmut+"\t"+nbseg+"\t"+nbnode+"\t"+treescore+"\t"+nbmutfirst+"\t"+nbmutrest+"\n"
        output+=ss
    
output_f.write("%s"%output)
output_f.close()
