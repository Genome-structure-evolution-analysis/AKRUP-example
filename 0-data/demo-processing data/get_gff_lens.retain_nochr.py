# -*- coding: utf-8 -*-
# @Author: wjq
# @Date:   2021-10-22 14:56:21
# @Last Modified by:   wjq
# @Last Modified time: 2023-07-18 23:21:43

import pandas as pd
from Bio import SeqIO

def isnum(gff):
    gff_list = gff.to_list()
    for i in range(len(gff_list)):
        if not isinstance(gff_list[i], int) and gff_list[i].isnumeric():
            gff_list[i] = int(gff_list[i])
    args = pd.Series(gff_list)
    return args

def get_gff_lens(gff3_path, name_pos, chr_perfix, spec_name, skiprows, content):
    data = pd.read_csv(gff3_path, sep="\t", header=None,skiprows=skiprows, comment='#')
    data = data[data[2] == 'mRNA'].loc[:, [0, 8, 3, 4, 6]]
    data[8] = data[8].str.split(':|=|;', expand=True)[name_pos]

    data[0] = data[0].str.replace(chr_perfix, '', regex=True).replace('^0', '', regex=True)
    data = data[data[8].str.contains('\.1$')]  #  Removing Variable Shear

    data[3] = data[3].astype('int')
    for name, group in data.groupby([0]):
        group = group.sort_values(by=[3])
        data.loc[group.index, 'order'] = list(range(1, len(group)+1))
        data.loc[group.index, 'newname'] = list(
            # [spec_name+str(name).zfill(2)+'g'+str(i).zfill(5) for i in range(1, len(group)+1)]) 
            [spec_name+str(name).zfill(1)+'g'+str(i).zfill(4) for i in range(1, len(group)+1)])
    data[[4, 'order']] = data[[4, 'order']].astype('int')  # change
    data = data[[0, 3, 4, 6, 8, 'newname', 'order']].sort_values(by=[0, 'order'], key=isnum)
    # data.to_csv(f'{spec_name}.sca.gff', sep="\t", index=False, header=None)
    data = data[~data[0].str.contains(content)]
    data.to_csv(f'{spec_name}.new.gff', sep="\t", index=False, header=None)
    lens = data.groupby(0).max()[[4, 'order']].sort_values(by=[0], key=isnum)
    # data[~data[0].str.contains(content)].to_csv(f'{spec_name}.new.gff', sep="\t", index=False, header=None)
    # lens = data[~data[0].str.contains(content)].groupby(0).max()[[4, 'order']].sort_values(by=[0], key=isnum)
    lens = lens[lens['order'] > 100]
    lens.to_csv(f'{spec_name}.lens.txt', sep="\t", header=None)
    lens['order'].to_csv(f'{spec_name}.lens', sep="\t", header=None)
    print('end!!!')
    return data


def ch_pep_cds(data, fas_f, sf_name, wrap, num=0):
    id_dict = data.set_index([8])['newname'].to_dict()
    writer = SeqIO.FastaIO.FastaWriter(open(sf_name, 'w'), wrap=wrap)
    for seq_record in SeqIO.parse(fas_f, "fasta"):
        if seq_record.id in id_dict:
            num += 1
            seq_record.id = id_dict[seq_record.id]
            seq_record.description = ''
            seq_record.seq = seq_record.seq.upper()
            writer.write_record(seq_record)
    print(f'All seq: {num}')


if __name__ == '__main__':

    wrap = 60  # Retained bases or amino acids in each line
    skiprows = 3  # The number of lines skipped
    spec_name = 'Os'  # Species abbreviation
    gff_name_pos = 3
    chr_perfix = 'Chr'  # Replace delete string
    
    content = 'Un|Sy'

    pep_file = 'Osativa_323_v7.0.protein.fa'
    cds_file = 'Osativa_323_v7.0.cds.fa'
    gff_file = 'Osativa_323_v7.0.gene.gff3'


    new_gff = get_gff_lens(gff_file, gff_name_pos, chr_perfix, spec_name, skiprows, content)
    ch_pep_cds(new_gff, pep_file, f'{spec_name}.pep', wrap)
    ch_pep_cds(new_gff, cds_file, f'{spec_name}.cds', wrap)
