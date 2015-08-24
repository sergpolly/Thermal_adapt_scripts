import re
import os
import sys
from Bio import Seq
from Bio import SeqIO
from Bio import Entrez as ncbi
import xml.etree.ElementTree as ET
import pandas as pd
import itertools




def get_temperature(text):
    # options are:
    # 1) 35-40
    # 2) 98
    # 3) 100C
    # 4) 34~45
    # 5) 35-40C
    to_str = lambda txt: str(txt) if type(txt)!=unicode else txt
    text = to_str(text)
    patterns = ['(\d+)-(\d+)','(\d+)-(\d+)C','(\d+)~(\d+)','(\d+)~(\d+)C','(\d+)','(\d+.\d+)','(\d+)C',u'(\d+)[\xb0]C','<(\d+)-(\d+)']
    for pattern in patterns:
        res = re.match(pattern,text)
        if res:
            return pd.np.mean([float(temp) for temp in res.groups()])
    return pd.np.nan




def get_env_features(xml_root):
    to_return = []
    xml_entries = list(xml_root.getchildren())
    for entry_root in xml_entries:
        uid = entry_root.attrib['uid']
        xml_org = list(entry_root.iter('OrganismName'))
        # organism name isn't a cumpulsory option for these xml files ...
        # so we have avoid this problem somehow...
        if not xml_org:
            return to_return
        else:
            xml_org = xml_org[0]
        org_name = xml_org.text
        # let's find Environmental description ...
        xml_env = list(entry_root.iter('Environment'))
        if not xml_env:
            print >> sys.stdout, "No 'Environment'-entry for %s in the BioProject %s"%(org_name, uid)
        else:
            print >> sys.stdout, "There is %d 'Environment' entry in %s for %s, as expected!"%(len(xml_env), uid, org_name)
            dict_to_append = [('BioProject',uid),]
            dict_to_append += [('Organism',org_name),]
            dict_to_append += [(node.tag, node.text) for node in xml_env[0].getchildren()]
            to_return.append(dict(dict_to_append))
    return to_return



def xmlSplitter(data,separator=lambda x: x.startswith('<?xml')):
  buff = []
  for line in data:
    if separator(line):
      if buff:
        yield ''.join(buff)
        buff[:] = []
    buff.append(line)
  yield ''.join(buff)



# path = os.path.join(os.path.expanduser('~'),'GENOMES_BACTER/ftp.ncbi.nih.gov/refseq/release/bacteria/genbank')
path = os.path.join(os.path.expanduser('~'),'GENOMES_BACTER_RELEASE69/genbank')

# descriptions_fname = "genbank.inventory.description"
results_xml = "results.xml"


# now iterate over all xml files and whatever is inside those files ...
xml_roots = []
with open(os.path.join(path,results_xml),'r') as fp:
    for batch in xmlSplitter(fp):
        xml_roots.append(ET.fromstring(batch))

# now collecting environmental features about these bacteria ...
results = []
for xml_root in xml_roots:
    results += get_env_features(xml_root)


env_df = pd.DataFrame(results)
env_df = env_df.drop_duplicates(subset=['BioProject','Organism','OptimumTemperature'])


notnnull_temp_entries = env_df.OptimumTemperature.notnull().nonzero()[0].size
entries_parsed = env_df.OptimumTemperature.apply(get_temperature).notnull().nonzero()[0].size
print >> sys.stdout, "There are %d non-empty OptTemp entries in the DataFrame, and %d of them are parseable."%(notnnull_temp_entries,entries_parsed)


env_df.OptimumTemperature = env_df.OptimumTemperature.apply(get_temperature)

# we now for sure that there are 2 entries for '[Clostridium] stercorarium subsp. stercorarium DSM 8532'
# they both have OptTemps that differ slightly 60 and 65 ...
# let's just take one of these:
env_df = env_df[env_df.OptimumTemperature.notnull()]

env_df = env_df.drop_duplicates(subset=['BioProject',],take_last=True)
print >> sys.stdout, "There are %d unique entries with OptTemp after all the adjustments ..."%(env_df.shape[0])


temp_vals_notnull = env_df.OptimumTemperature.notnull().nonzero()[0].size
under50 = env_df.OptimumTemperature.where(env_df.OptimumTemperature<=50.0).notnull().nonzero()[0].size
above50 = env_df.OptimumTemperature.where(env_df.OptimumTemperature>50.0).notnull().nonzero()[0].size
print >> sys.stdout, ''
print >> sys.stdout, "There are %d notnull values of temperature retrieved. %d are above 50C; %d are under 50C(incl.)"%(temp_vals_notnull,above50,under50)
print >> sys.stdout, ''
print >> sys.stdout, 'The temp values are unique in a sense, that each organism is displayed <= once.'

#print which columns do we have:
print >> sys.stdout, "Columns of the env_DataFrame are: ",env_df.columns

env_df[['BioProject', 'Organism', 'OptimumTemperature', 'TemperatureRange', 'OxygenReq', 'Habitat', 'Salinity']].to_csv(os.path.join( path, 'env_catalog_bacter.dat' ), index=False)





# ###########################################################
# # env data extracted using links to BioProject entries ...
# ###########################################################
# print >> sys.stdout, ''
# print >> sys.stdout, 'Environmental data extracted using links to the BioProject uids ...'
# subset = ['BP_Organism','Habitat','OptimumTemperature','OxygenReq','Salinity','TemperatureRange']
# link_dat = pd.read_csv(os.path.join(org_path,'environmental_catalog.dat'))
# link_dat_uniq = link_dat.drop_duplicates(subset=subset)
# link_dat_uniq.OptimumTemperature = link_dat_uniq.OptimumTemperature.apply(get_temperature)
# l_temp_vals_notnull = link_dat_uniq.OptimumTemperature.notnull().nonzero()[0].size
# l_under50 = link_dat_uniq.OptimumTemperature.where(link_dat_uniq.OptimumTemperature<=50.0).notnull().nonzero()[0].size
# l_above50 = link_dat_uniq.OptimumTemperature.where(link_dat_uniq.OptimumTemperature>50.0).notnull().nonzero()[0].size
# print >> sys.stdout, ''
# print >> sys.stdout, "There are %d notnull values of temperature retrieved. %d are above 50C; %d are under 50C(incl.)"%(l_temp_vals_notnull,l_above50,l_under50)
# print >> sys.stdout, ''
# print >> sys.stdout, 'The temp values are unique in a sense, that each organism is displayed <= once.'
# link_dat_uniq[['BP_Organism', 'OptimumTemperature', 'TemperatureRange', 'OxygenReq', 'Habitat', 'Salinity']].rename(columns={'BP_Organism':'Organism'}).to_csv(os.path.join( org_path, 'env_catalog_LINK.dat' ), index=False)















