#!/usr/bin/Python

# creating a matrix containing FPKM and TPM values for all transcripts (rows) across all Samples (columns)
import pandas as pd
import numpy as np
import glob
import os

n = 10 ** 5 # memory error for a chunksize bigger than this

# using chunks to load only that chunk into the memory for reading instead of loading the whole file at once

i = 1
for chunk in pd.read_csv("all-transcripts/all-genes-tx-abundance.txt", sep = "\t", chunksize = n):	# Size of all-genes-tx-abundance.txt = 8G
	for n,g in chunk.groupby('Sample'):
		g = pd.pivot_table(g,index=["Tx_id"], columns = ["Sample"], values = ["FPKM","TPM"])
		g.fillna(0,inplace = True)

		#final_lst.append(g)	
		# appending dataframe chunks and later concating all dataframes requires more RAM and processes are killed. 
		# Hence, appending dataframe chunks to CSV's

		g.to_csv("all-transcripts/"+str(i)+"all-genes-abundance-matrix.txt", sep = "\t")
		i += 1


# 2448 files, average size of each file betwwen 4-5Mb

files = glob.glob("all-transcripts/*all-genes-abundance-matrix.txt")	

# Loading 2448 * 5Mb = 12240Mb ~ 12.24G into the memory is computer intensive and processes are killed
# when trying to read all sub-files into sub-dataframes and concating them all together.
# Hence, Dividing files list into four equal sized parts

one_fourth = len(files)//4
half = len(files)//2
three_fourth = (len(files) * (3/4))



print("Concating sub-dataframes...")	
# so not all files are loaded into the memory together and taking up all the RAM

df1 = pd.concat([pd.read_csv(f) for f in files[:int(one_fourth)]], axis = 1)
#df1.to_csv("df1-genes-abundance-matrix.txt", sep ="\t")
#del(df1)

df2 = pd.concat([pd.read_csv(f) for f in files[int(one_fourth):int(half)]], axis = 1)
#df2.to_csv("df2-genes-abundance-matrix.txt", sep ="\t")
#del(df2)

df3 = pd.concat([pd.read_csv(f) for f in files[int(half):int(three_fourth)]], axis = 1)
#df3.to_csv("df3-genes-abundance-matrix.txt", sep ="\t")
#del(df3)

df4 = pd.concat([pd.read_csv(f) for f in files[int(three_fourth):]], axis = 1)
#df4.to_csv("df4-genes-abundance-matrix.txt", sep ="\t")
#del(df4)

# removing all previously generated CSV's to free up some memory
print("Removing previously generated CSV's...")
os.system("rm *all-genes-abundance-matrix.txt")



print("Merging df1 and df2...")
final1 = pd.concat([df1,df2])
del(df1,df2)
print("Merging df3 and df4...")
final2 = pd.concat([df3,df4])
del(df3,df4)

print("Merging final1 and final2...")
final_merge = pd.concat([final1,final2])
del(final1,final2)

final_merge.to_csv("Final-all-Tx-abundance-matrix.txt", sep = "\t")	#5.19G

print("Finished!")