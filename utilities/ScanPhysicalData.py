#!/usr/bin/env python2
# -*- encoding: utf8 -*-

if __name__ == "__main__":

	import argparse
	import glob
	import math as m
	parser = argparse.ArgumentParser(description='Berechnet Größen aus den physikalischen Daten.')
	parser.add_argument('path', help='Dateien welche durchsucht werden sollen. Bei der Verwendung von regex muss das Argument als string übergeben werden.')
	args = parser.parse_args()
	files=glob.glob(args.path)

	for file in files: 
		stream=open(file, 'r')
		[timestamp, MB, M_common, A_common, J_common, M_up, A_up, J_up, M_down, A_down, J_down, pd, pc_common, pc_up, pc_down] = stream.readline().split()
		MB = float(MB)
		M_common = float(M_common)
		A_common = float(A_common)
 		J_common = float(J_common)
 		M_up = float(M_up)
 		A_up = float(A_up)
 		J_up = float(J_up)
 		M_down = float(M_down)
 		A_down = float(A_down)
 		J_down = float(J_down)
#		print [timestamp, MB, M_common, A_common, J_common, M_up, A_up, J_up, M_down, A_down, J_down]
		penrose=16*m.pi*MB*MB - (A_up + A_down)
		dain_bondi = (A_up + A_down)/(8*m.pi*(MB*MB + m.sqrt(MB*MB*MB*MB - (J_down + J_up)*(J_down + J_up))))
		print file + " \tPenrose (>=0): " + str(penrose) + "\t Dain-Bondi (>0 & <=1): " + str(dain_bondi)
