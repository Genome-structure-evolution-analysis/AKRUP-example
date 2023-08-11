# -*- coding: utf-8 -*-
# @Author: wjq
# @Date:   2023-07-27 04:50:05
# @Last Modified by:   wjq
# @Last Modified time: 2023-08-11 17:18:54

import os
import glob


files = glob.glob('*.txt')
for file in files:
	name = os.path.basename(file)
	for ty in ['png', 'pdf']:
		sf = open('ancestral_plotfig-auto.conf', 'w')
		text = f"""[ancestralfig]
name = 
anc_reverse = True
frame_flag = True
chr_name = False
ancestor_file = {name}
savefile = {name}.{ty}
"""
		sf.write(text)
		sf.close()
		os.system('AKRUP -akf ancestral_plotfig-auto.conf')


