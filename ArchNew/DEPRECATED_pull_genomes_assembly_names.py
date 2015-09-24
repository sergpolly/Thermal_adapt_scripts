from Bio import Entrez as ncbi
from Bio import SeqIO
import pandas as pd
import os
from time import sleep
import subprocess as sub

def get_tax_id(species):
    """to get data from ncbi taxomomy, we need to have the taxid. we can
    get that by passing the species name to esearch, which will return
    the tax id"""
    species = species.replace(' ', "+").strip()
    search = ncbi.esearch(term = species, db = "taxonomy", retmode = "xml")
    record = ncbi.read(search)
    return record['IdList'][0]

def get_tax_data(taxid):
    """once we have the taxid, we can fetch the record"""
    search = ncbi.efetch(id = taxid, db = "taxonomy", retmode = "xml")
    return ncbi.read(search)

# looks like it's the wrong way of doing it ...
# rather use fetched genomes from ftp, but
# mind that there are links instead of genomes there, so still would need to do a lot of wget by the python ...
# maybe there would be no need for biopython stuff though ...


root_path = "."
# archaea/genomes
data_path = os.path.join(root_path,"archaea")
asm_path = os.path.join(root_path,"archaea_links")
result_path = os.path.join(root_path,"archaea/genomes")

# let ncbi know who you are ...
ncbi.email = "sergey.venev@umassmed.edu"

# read table with organisms and their taxonid ...
dat = pd.read_csv(os.path.join(data_path,"termoprop.txt"),sep='\t')


folders = sub.check_output("ls archaea_links",shell=True).rstrip().split('\n')
folders = map(lambda x: x.replace('_',' '),folders)

corresponding_fnames = []

for name,txid in dat[['name','taxonid']].values:
    if name not in folders:
        # get scientific name for the thing ...
        tax_data = get_tax_data(str(txid))[0]
        scientific_name = tax_data['ScientificName']
        if scientific_name not in folders:
            # get taxon parent's name ...
            parent_id = tax_data['ParentTaxId']
            parent_name = get_tax_data(parent_id)[0]['ScientificName']
            if parent_name not in folders:
                print "even parent not in folders!!!"
                print "(%d,%s) -> (%s,%s) not found"%(txid,name,parent_id,parent_name)
            else:
                corresponding_fnames.append(parent_name)
        else:
            corresponding_fnames.append(scientific_name)
    else:
        corresponding_fnames.append(name)


dat['fname'] = map(lambda x: x.replace(' ','_'),corresponding_fnames)



sub.check_output("ls -p %s | grep -v /"%os.path.join(asm_path,yyy),shell=True).strip().split()



links = pd.read_csv("folder_assembly.txt")
links.set_index('folder',inplace=True)
# mind the indexing differences - they are important!!!!!!!!!!!
dat['assembly'] = links.loc[dat.fname].values

dat.to_csv("tax_temp_asm.txt")

# ################################

# import ftplib

# path = 'genomes/all'
# # filename = 'L28POC_B.xpt'

# ftp = ftplib.FTP("ftp.ncbi.nih.gov") 
# ftp.login("anonymous", "sergey.venev@umassmed.edu") 
# ftp.cwd(path)
# # check needed foles - protein fasta, or genome genbank and fetch it!
# ftp.retrbinary("RETR " + filename ,open(filename, 'wb').write)
# ftp.quit()


# for dd in dirs:
#     # sub.check_output("ls",shell=True).strip().split()
#     inside = sub.check_output("ls %s"%dd,shell=True).strip().split()
#     if 'representative' in inside:
#         links = sub.check_output("ls %s/representative"%dd,shell=True).strip().split()
#         if len(links)>1:
#             print "several links in representative of %s, taking first: %s"%(dd,links[0])
#         with open(dd+'/'+links[0],'w') as fp:
#             fp.write("link to ftp.ncbi.nih.gov/genomes/genbank/all/fname ...\n")
#     elif 'latest_assembly_versions' in inside:
#         links = sub.check_output("ls %s/latest_assembly_versions"%dd,shell=True).strip().split()
#         if len(links)>1:
#             print "several links in latest_assembly_versions of %s, taking first: %s"%(dd,links[0])
#         with open(dd+'/'+links[0],'w') as fp:
#             fp.write("link to ftp.ncbi.nih.gov/genomes/genbank/all/fname ...\n")
#     elif 'all_assembly_versions' in inside:
#         links = sub.check_output("ls %s/all_assembly_versions"%dd,shell=True).strip().split()
#         if len(links)>1:
#             print "several links in all_assembly_versions of %s, taking first: %s"%(dd,links[0])
#         with open(dd+'/'+links[0],'w') as fp:
#             fp.write("link to ftp.ncbi.nih.gov/genomes/genbank/all/fname ...\n")




# import subprocess as sub
# dirs = sub.check_output("ls",shell=True).strip().split()
# pairs = []
# for dd in dirs:
#     # sub.check_output("ls",shell=True).strip().split()
#     inside = sub.check_output("ls %s"%dd,shell=True).strip().split()
#     if 'representative' in inside:
#         links = sub.check_output("ls %s/representative"%dd,shell=True).strip().split()
#         if len(links)>1:
#             print "several links in representative of %s, taking first: %s"%(dd,links[0])
#         # with open(dd+'/'+links[0],'w') as fp:
#         #     fp.write("link to ftp.ncbi.nih.gov/genomes/genbank/all/fname ...\n")
#         pairs.append((dd,links[0]))
#     elif 'latest_assembly_versions' in inside:
#         links = sub.check_output("ls %s/latest_assembly_versions"%dd,shell=True).strip().split()
#         if len(links)>1:
#             print "several links in latest_assembly_versions of %s, taking first: %s"%(dd,links[0])
#         # with open(dd+'/'+links[0],'w') as fp:
#         #     fp.write("link to ftp.ncbi.nih.gov/genomes/genbank/all/fname ...\n")
#         pairs.append((dd,links[0]))
#     elif 'all_assembly_versions' in inside:
#         links = sub.check_output("ls %s/all_assembly_versions"%dd,shell=True).strip().split()
#         if len(links)>1:
#             print "several links in all_assembly_versions of %s, taking first: %s"%(dd,links[0])
#         # with open(dd+'/'+links[0],'w') as fp:
#         #     fp.write("link to ftp.ncbi.nih.gov/genomes/genbank/all/fname ...\n")
#         pairs.append((dd,links[0]))
# with open("folder_assembly.txt",'w') as fp:
#     fp.write("folder,assembly\n")
#     for dd, link in pairs:
#         fp.write("%s,%s\n"%(dd,link))




# # import os
# #     import ftplib
# #     from contextlib import closing

# #     with closing(ftplib.FTP()) as ftp:
# #         try:
# #             ftp.connect(host, port, 30*5) #5 mins timeout
# #             ftp.login(login, passwd)
# #             ftp.set_pasv(True)
# #             with open(local_filename, 'w+b') as f:
# #                 res = ftp.retrbinary('RETR %s' % orig_filename, f.write)

# #                 if not res.startswith('226 Transfer complete'):
# #                     print('Downloaded of file {0} is not compile.'.format(orig_filename))
# #                     os.remove(local_filename)
# #                     return None

# #             return local_filename

# #         except:
# #                 print('Error during download from FTP')





































