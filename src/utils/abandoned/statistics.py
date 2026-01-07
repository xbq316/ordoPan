
import os

filedir="/mnt/f/project/SG-2025.6.25/agora"
dirs=os.listdir(filedir)

a=[]
for i in dirs:
    if i.startswith("imitate"):
        a.append(i)

flag=False
for i in a:
    x=i.split(".")
    filename=f"{filedir}/{i}/1.log"
    for j in open(filename):
        j=j.strip()
        if j.startswith("real"):
            line=j.split("\t")
            mins=line[1].split("m")[0]
            sec=line[1].split("m")[1].split("s")[0]
            time=(float(mins)*60+float(sec))/60
            flag=True
        zz=j.split(" ")
        if len(zz)==3 and flag:
            line=[float(z) for z in zz]
            if time <1:
                time=1
            else:
                time=round(time,2)
            print(f"{x[1]} {x[2]} {x[3]} {line[0]} {line[1]} {time}")
            flag=False
            break
        