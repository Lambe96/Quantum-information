import os

N_min = input('Enter N_min:')
N_min = int(N_min)
N_max = input('Enter N_max:')
N_max = int(N_max)
step = input('Enter step:')
step = int(step)

for i in range(N_min,N_max,step):
    f= open("dim.txt","w")
    f.write('%d,%d,%d,%d' % (i,i,i,i))
    f.close()
    os.system("./ex04")
    print(i)
