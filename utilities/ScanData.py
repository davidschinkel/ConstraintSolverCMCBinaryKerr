#!/usr/bin/env python2
# -*- encoding: utf8 -*-

if __name__ == "__main__":

	import argparse
	import glob
	import numpy as np
	parser = argparse.ArgumentParser(description='Gibt die Konfigurationsdaten aus.')
	parser.add_argument('path', help='Dateien welche durchsucht werden sollen. Bei der Verwendung von regex muss das Argument als string übergeben werden.')
	args = parser.parse_args()
	files=glob.glob(args.path)

	data = []
	
	for file in files: 
		stream=open(file, 'r')

		a0, Rp, Ktr, rH1, c1, Sz1, Pz1, rH2, c2, Sz2, Pz2, trapping_surface, ns_domain_0, ns_domain_1, aux , nt = [stream.readline().partition(' ')[0] for i in xrange(16)]

		newrow = [a0, Rp, Ktr, rH1, c1, Sz1, Pz1, rH2, c2, Sz2, Pz2, trapping_surface, ns_domain_0, ns_domain_1, nt]

		condition = (float(rH1)==1.3333 and float(rH2)==0.6666 and float(a0)==4.8 and float(Ktr)==0.05)

		if(condition):
			newrow.insert(0,file)
			data.append(newrow)
			#print file + " \t " + a0 + " \t " + Rp + " \t " + Ktr  + " \t " +  rH1  + " \t " +  c1  + " \t " + Sz1  + " \t " + Pz1  + " \t " + rH2  + " \t " +  c2  + " \t " + Sz2  + " \t " + Pz2  + " \t " + trapping_surface  + " \t " + ns_domain_0  + " \t " + ns_domain_1  + " \t " + nt
			print file + " \t " + a0 + " \t " +  rH1  + " \t " +  c1  + " \t " + Sz1  + " \t " + Pz1  + " \t " + rH2  + " \t " +  c2  + " \t " + Sz2  + " \t " + Pz2  
			#print "gnuplot -e \"a=" + file  + "\" ../../utilities/GnuplotScript.gnu"

	npdata = np.array(data)


	#print npdata[:,1].argsort()
	
	#print np

	#sorting the data
	#print data.sort(key = lambda x: x[1])
	#for i in npdata:
	#	print i