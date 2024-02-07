f=open("inputFiles.txt", "r")

lines=f.readlines()

ids=[]
for line in lines:
    ids.append(line.split()[0])

count=0
count2=0
for i in range(len(ids)):
    count+=1
    if count==10:
        count=0
        print("#tbr" + str(count2-8) + "PDBids=(" + ids[count2-9], ids[count2-8], ids[count2-7], ids[count2-6], ids[count2-5], ids[count2-4], ids[count2-3], ids[count2-2], ids[count2-1], ids[count2]+")") 
    count2+=1

f.close()
