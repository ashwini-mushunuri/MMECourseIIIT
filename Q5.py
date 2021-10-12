def plotting(time,p_values):
	plt.plot(time,p_values)
	plt.xlabel("Time")
	plt.ylabel("p(t)=0.75-0.75*((1-(4*alpha/3))^t)         alpha=0.01")
	plt.title("Sites that differ from the original base as a function of time by Jukes-Cantor model")
	plt.show()
	plt.plot(time,p_values)
	plt.xlabel("Time")
	plt.ylabel("p(t)=0.75-0.75*((1-(4*alpha/3))^t)         alpha=0.01")
	plt.title("Sites that differ from the original base as a function of time by Jukes-Cantor model")
	plt.savefig("5.jpg")



def main():
	if os.path.isdir("20172108"):
		os.chdir("20172108")
	else:
		os.mkdir("20172108")
		os.chdir("20172108")
	f=open("p_values.csv","w")
	f.write("p(t)=0.75-0.75*((1-(4*alpha/3))^t)\n")
	f.write("alpha=0.01\n")
	f.write("Time,p(t)\n")
	time,p_values=[],[]
	for values in range(0,501):
		time.append(values)
		p=0.75-0.75*math.pow((1-(4*0.01/3)),values)
		p_values.append(p)
		f.write(str(values))
		f.write(","+str(p)+"\n")
	f.close()
	plotting(time,p_values)
	print "p(t) values as a function of time are stored in p_values.csv file in 20172028_20172072 folder"
	print "Plot is also saved in the file 5.jpg in 201460579 folder"



import math, os
try:
	import matplotlib.pyplot as plt
except:
	"!!!matplotlib package not installed!!!"
else:
	main()
