vector=[1,2,3,4,5,6,7,8,9]
temp= []
temp1=[]
lenvec=len(vector)
for i in range(lenvec):
    temp.append(0)
    temp1.append(0)

for z in range(lenvec):
    temp[z]= vector[z]+temp[z-1]
print ("The final answer is ", temp)

for y in range(lenvec):
    temp1[y]= vector[y]*temp[y]

print ("The final answer is ", temp1)
print ("The final answer is ", vector)  
