def kimura3(seq1,seq2):
	f=open("distances.txt","a")
	f.write("\nKimura 3 Parameter Distance Details...\n")
	beta_count,del_count,gamma_count,pairs=0,0,0,[]
	beta_values,gamma_values,del_values=["AG","GA","CT","TC"],["AC","CA","GT","TG"],["AT","TA","GC","CG"]
	for values in zip(seq1,seq2):
		if '-' not in values:
			pairs.append(values)
	for (base1,base2) in pairs:
		if base1+base2 in beta_values:
			beta_count+=1
		elif base1+base2 in gamma_values:
			gamma_count+=1
		elif base1+base2 in del_values:
			del_count+=1
	print "Done..."
	beta=beta_count/float(len(pairs))
	delta=del_count/float(len(pairs))
	gamma=gamma_count/float(len(pairs))
	print "Beta=",beta
	print "Delta=",delta
	print "Gamma=",gamma
	f.write("Beta= "+str(beta)+"\n")
	f.write("Delta= "+str(delta)+"\n")
	f.write("Gamma= "+str(gamma)+"\n")
	try:
		k3_distance=-0.25*(math.log(1-2*beta-2*gamma)+math.log(1-2*beta-2*delta)+math.log(1-2*gamma-2*delta))
	except:
		print "!!! 1 < -(2b+2g) or 1 < -(2b+2d) or 1 < -(2g+2d); -ve logarithm cannot be computed!!!"
	else:
		print "Kimura 3 Distance=",k3_distance
		f.write("Kimura 3 Distance= "+str(k3_distance)+"\n")
	f.close()



def kimura2(seq1,seq2):
	f=open("distances.txt","a")
	f.write("\nKimura 2 Parameter Distance Details...\n")
	ts,tv,pairs,transition,transversion=0,0,[],["AG","GA","CT","TC"],["AC","CA","AT","TA","GC","CG","GT","TG"]
	for values in zip(seq1,seq2):
		if '-' not in values:
			pairs.append(values)
	for (base1,base2) in pairs:
		if base1+base2 in transition:
			ts+=1
		elif base1+base2 in transversion:
			tv+=1
	print "Done..."
	trans_freq=ts/float(len(pairs))
	tranv_freq=tv/float(len(pairs))
	print "Transition Frequency, p=",trans_freq
	print "Transversion Frequency, q=",tranv_freq
	f.write("Transition Frequency, p= "+str(trans_freq)+"\n")
	f.write("Transversion_Frequency, q= "+str(tranv_freq)+"\n")
	try:
		k2_distance=-0.5*math.log(1-2*trans_freq-tranv_freq)-0.25*math.log(1-2*tranv_freq)
	except:
		print "!!! 1 < -(2p+q) or 1 < 2q; -ve logarithm cannot be computed!!!"
	else:
		print "Kimura 2 Distance=",k2_distance
		f.write("Kimura 2 Distance= "+str(k2_distance)+"\n")
	f.close()



def jukes_cantor(seq1,seq2):
	f=open("distances.txt","a")
	f.write("\nJukes Cantor Distance Details...\n")
	b,difference,pairs=0.75,0,[]
	for values in zip(seq1,seq2):
		if '-' not in values:
			pairs.append(values)
	for (base1,base2) in pairs:
		if base1!=base2:
			difference+=1
	p=difference/float(len(pairs))
	try:
		jc_distance=-b*math.log(1-p/b)
	except:
		print "!!! 1 < (p/b); -ve logarithm cannot be computed!!!"
		f.write("Jukes Cantor Distance= NA")
	else:
		print "Done..."
		print "Proportion of sites with different nucleotides, p=",p
		f.write("Proportion of sites with different nucleotides, p= "+str(p)+"\n")
		print "Jukes Cantor Distance=",jc_distance
		f.write("Jukes Cantor Distance= "+str(jc_distance)+"\n")
	f.close()



def frequency(seq1,seq2):
	f=open("distances.txt","w")
	f.write("Frequency Distribution Details...\n")
	print "Done..."
	print "--------------------------------------------------------\n            Seq1                      Seq2\n--------------------------------------------------------"
	f.write("--------------------------------------------------------\n            Seq1                      Seq2\n--------------------------------------------------------\n")
	print "        A     T     G     C        A     T     G    C\n--------------------------------------------------------"
	f.write("        A     T     G     C        A     T     G    C\n--------------------------------------------------------\n")
	print "Count  ",seq1.count("A")," ",seq1.count("T")," ",seq1.count("G")," ",seq1.count("C"),"    ",seq2.count("A")," ",seq2.count("T")," ",seq2.count("G")," ",seq2.count("C")
	f.write("Count  "+str(seq1.count("A"))+"   "+str(seq1.count("T"))+"   "+str(seq1.count("G"))+"   "+str(seq1.count("C"))+"      "+str(seq2.count("A"))+"   "+str(seq2.count("T"))+"   "+str(seq2.count("G"))+"   "+str(seq2.count("C"))+"\n")
	print "Freq  ",round(seq1.count("A")/float(len(seq1)),3),round(seq1.count("T")/float(len(seq1)),3),round(seq1.count("G")/float(len(seq1)),3),round(seq1.count("C")/float(len(seq1)),3),"  ",round(seq2.count("A")/float(len(seq2)),3),round(seq2.count("T")/float(len(seq2)),3),round(seq2.count("G")/float(len(seq2)),3),round(seq2.count("C")/float(len(seq2)),3)
	f.write("Freq  "+str(round(seq1.count("A")/float(len(seq1)),3))+" "+str(round(seq1.count("T")/float(len(seq1)),3))+" "+str(round(seq1.count("G")/float(len(seq1)),3))+" "+str(round(seq1.count("C")/float(len(seq1)),3))+"    "+str(round(seq2.count("A")/float(len(seq2)),3))+" "+str(round(seq2.count("T")/float(len(seq2)),3))+" "+str(round(seq2.count("G")/float(len(seq2)),3))+" "+str(round(seq2.count("C")/float(len(seq2)),3))+"\n")
	print "--------------------------------------------------------"
	f.write("--------------------------------------------------------\n")



def read_file(file1,file2):
	seq1,seq2="",""
	f=open(file1,"r")
	for line in f:
		if line.startswith(">") or re.match("[a-z]",line):
			pass
		else:
			line=re.sub("\n","",line)
			seq1+=line
	f.close()
	f=open(file2,"r")
	for line in f:
		if line.startswith(">") or re.match("[a-z]",line):
			pass
		else:
			line=re.sub("\n","",line)
			seq2+=line
	f.close()
	return seq1,seq2



def main():
	if len(sys.argv)!=3:
		print "!!!Invalid input. Input should contain 2 input file names!!!"
		print "Usage: python 4.py <file_name1> <file_name2>"
	elif os.path.isfile(sys.argv[1]) and os.path.isfile(sys.argv[2]):
		seq1,seq2=read_file(sys.argv[1],sys.argv[2])
		if os.path.isdir("20172108"):
			os.chdir("20172108")
		else:
			os.mkdir("20172108")
			os.chdir("20172108")
		print "Calculating Frequency Distribution..."
		frequency(seq1,seq2)
		print "\nCalculating Jukes Cantor Distance..."
		jukes_cantor(seq1,seq2)
		print "\nCalculating Kimura-2 Distance..."
		kimura2(seq1,seq2)
		print "\nCalculating Kimura-3 Distance..."
		kimura3(seq1,seq2)
		print "\nAll results are stored in distances.txt file in 20172108 folder"
	else:
		print "!!!Input file(s) does not exist in working directory!!!"



import sys, os, re, math
main()
