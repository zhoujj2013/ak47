#!/usr/bin/python

import os, sys
import re

def usage():
	print '\nFor creating excel file by input plat text file\n'
	print 'Author: zhoujj2013@gmail.com 4/2/2017\n'
	print 'Usage: python '+sys.argv[0]+' header.txt con.txt result.xlsx\n'
	print 'Example input:'
	print '/x400ifs-accel/zhoujj/Project/15.MendelDeseaseDiagnosis/bin/testExcel/header.txt'
	print '/x400ifs-accel/zhoujj/Project/15.MendelDeseaseDiagnosis/bin/testExcel/con.txt'
	sys.exit(2)

def CreatExcel(meta_arr, header, con, filename):
	#pdir = os.path.split( os.path.realpath(sys.argv[0]))[0]
	#sys.path.append(pdir + "/XlsxWriter-master/build/lib")
	import xlsxwriter
	
	#print filename
	workbook = xlsxwriter.Workbook(filename)
	worksheet = workbook.add_worksheet()
	bold = workbook.add_format({'bold': True})
	numFormat = workbook.add_format()
	numFormat.set_num_format(0x00)
	
	meta_format = workbook.add_format({'bold': True, 'font_size': 12, 'bg_color': 'yellow'})
	
	row = 0
	for mi in range(len(meta_arr)):
		for ci in range(len(meta_arr[mi])):
			worksheet.write(row, ci, meta_arr[mi][ci],meta_format)
		row += 1
	row += 1
	
	for r in con:
		for ci in range(len(r)):
			#print ci
			if r[ci] in header:
				worksheet.write(row, ci , r[ci], bold)
				worksheet.write_comment(row, ci, header[r[ci]], {'author':'zhoujj zhoujiajian[at]link.cuhk.edu.hk', 'width':200, 'height':200})
			elif ci == 0:
				worksheet.write(row, ci , r[ci], bold)
			elif r[ci].replace('.','',1).isdigit() == True:
				#print r[ci]
				worksheet.write(row, ci , float(r[ci]), numFormat)
			else:
				worksheet.write(row, ci , r[ci])

		row += 1
	workbook.close()


if __name__ == "__main__":
	
	if len(sys.argv) < 3:
		usage()
	
	meta_f = open(sys.argv[1], 'rb')
	header_f = open(sys.argv[2], 'rb')
	con_f = open(sys.argv[3], 'rb')
	fname = sys.argv[4]
	
	m = {}
	marr = []
	while True:
		l = meta_f.readline()
		if len(l) == 0:
			break
		l = l.strip("\n")
		lc = l.split("=")
		marr.append(lc)
		m[lc[0]] = lc[1]

	h = {}
	while True:
		l = header_f.readline()
		if len(l) == 0:
			break
		l = l.strip("\n")
		lc = l.split("=")
		h[lc[0]] = lc[1]
	c = []
	while True:
		l = con_f.readline()
		if len(l) == 0:
			break
		l = l.strip("\n")
		lc = l.split("\t")
		c.append(lc)
	#print c
	CreatExcel(marr,h,c,fname)
