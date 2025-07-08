import numpy as np
import matplotlib.pyplot as plt 

nx = [5,10,20,40]
error = [0.0439665423400700,0.0140930018928378 ,0.0037103748109555 ,0.0007727707569152 ]

print((np.log(error[3]) - np.log(error[0]))/(np.log(nx[3]) - np.log(nx[0])))
ny = [5,10,20,40]
error2 = [0.0227961177358106 ,0.0065313696729424, 0.0016595538865112,0.0003418496659838]

plt.title("Convergence of nx and ny")
##plt.loglog(nx,error)
#plt.scatter(nx,error, label = "Data : nx")
#plt.loglog(ny,error2)
#plt.scatter(ny,error2, label = "Data : ny")
#plt.xlabel("Number of space gird points (nx/ny)   [-]")
#plt.ylabel("Relative error   [%]")
#plt.legend()


nt = [300,400,500,600,700,800,900]
error3 = [0.0008625199473337,0.0005547827853593,0.0003699770232044,0.0002467053761400,0.0001586210830352,0.0000925397904067,0.0000411325349956 ]


nt2= [200,400,600,800,1000,2000]
error4 = [0.0000676641066790,0.0000331287124788,0.0000216228862937,0.0000158710901831,0.0000124203709238,0.0000055197379736    ]

plt.title("Convergence of nt")
plt.xlabel("Number of time grid points (nt)   [-]")
plt.ylabel("Relative error   [%]")
plt.loglog(nt2,error4)
plt.scatter(nt2,error4,label = "Data : nt")
print((np.log(error3[5]) - np.log(error3[0]))/(np.log(nt[5]) - np.log(nt[0])))
plt.legend()
plt.savefig("conv_t.pdf")
plt.show()