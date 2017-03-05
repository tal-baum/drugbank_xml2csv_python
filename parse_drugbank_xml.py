from lxml import etree
from StringIO import StringIO
from collections import defaultdict
import csv

# Parse XML
f = open('drugbank.xml','r')
data = f.read()
f.close()

tree = etree.parse(StringIO(data))
context = etree.iterparse(StringIO(data))

root = tree.getroot()
# print len(root), 'drugs'


#######################################################################
# Iterate over drugs
drug2attrib = defaultdict(dict)
# drugbank_id -> {'drugname', 'drug_type', 'groups', 'targets/enzymes/transporters': [_id, _actions]}

target2attrib = defaultdict(dict)
enzyme2attrib = defaultdict(dict)
transporter2attrib = defaultdict(dict)
# drugbank_target_id -> {'gene', 'name', 'organism', 'taxonomy_id', 'uniprot_id', 'genbank_gene_id', 'genbank_protein_id', 'hgnc_id'}

tag_prefix = '{http://www.drugbank.ca}'

for child in root:
    for s in child.findall(tag_prefix+'drugbank-id'):
        if 'primary' in s.attrib:
            drugbank_id = s.text

    drugname = child.findall(tag_prefix+'name')[0].text
    drug2attrib[drugbank_id]['drugname'] = drugname

    drug_type = child.attrib['type']
    drug2attrib[drugbank_id]['drug_type'] = drug_type

    groups = [s.text for s in child.find(tag_prefix+'groups').findall(tag_prefix+'group')]
    drug2attrib[drugbank_id]['groups'] = groups

    # Get targets
    drug2attrib[drugbank_id]['targets'] = []
    for target in child.find(tag_prefix+'targets').findall(tag_prefix+'target'):
        if target.find(tag_prefix+'polypeptide') is None:
            continue
        
        target_id = target.find(tag_prefix+'id').text
        target_gene = target.find(tag_prefix+'polypeptide').find(tag_prefix+'gene-name').text
        target_name = target.find(tag_prefix+'name').text
        target_organism = target.find(tag_prefix+'organism').text
        target_taxonomy_id = target.find(tag_prefix+'polypeptide').find(tag_prefix+'organism').attrib['ncbi-taxonomy-id']
        
        if target_organism is None and target_taxonomy_id == '9606':
            target_organism = 'Human'
        if target_organism == 'Human' and target_taxonomy_id == '':
            target_taxonomy_id = '9606'
        if target_organism == 'Homo sapiens':
            target_organism = 'Human'
        if target_gene is None or target_organism is None or target_taxonomy_id is None:
            continue
        
        target_external_ids = target.find(tag_prefix+'polypeptide').find(tag_prefix+'external-identifiers').findall(tag_prefix+'external-identifier')
        target_uniprot_id = ''
        target_genbank_gene = ''
        target_genbank_protein = ''
        target_hgnc_id = ''

        for external_id in target_external_ids:
            if external_id.find(tag_prefix+'resource').text == 'UniProtKB':
                target_uniprot_id = external_id.find(tag_prefix+'identifier').text
            elif external_id.find(tag_prefix+'resource').text == 'GenBank Gene Database':
                target_genbank_gene = external_id.find(tag_prefix+'identifier').text
            elif external_id.find(tag_prefix+'resource').text == 'GenBank Protein Database':
                target_genbank_protein = external_id.find(tag_prefix+'identifier').text
            elif external_id.find(tag_prefix+'resource').text == 'HUGO Gene Nomenclature Committee (HGNC)':
                target_hgnc_id = external_id.find(tag_prefix+'identifier').text

        target_actions = [s.text.lower() for s in target.find(tag_prefix+'actions').findall(tag_prefix+'action')]
        
        drug2attrib[drugbank_id]['targets'].append((target_id,target_actions))

        if target_id not in target2attrib: #{'gene', 'name', 'organism', 'taxonomy_id', 'uniprot_id', 'genbank_gene_id', 'genbank_protein_id', 'hgnc_id'}
            target2attrib[target_id]['gene'] = target_gene
            target2attrib[target_id]['name'] = target_name
            target2attrib[target_id]['organism'] = target_organism
            target2attrib[target_id]['taxonomy_id'] = target_taxonomy_id
            
            target2attrib[target_id]['uniprot_id'] = target_uniprot_id
            target2attrib[target_id]['genbank_gene_id'] = target_genbank_gene
            target2attrib[target_id]['genbank_protein_id'] = target_genbank_protein
            target2attrib[target_id]['hgnc_id'] = target_hgnc_id

        #print target_id, target_gene, target_name, target_organism, target_taxonomy_id, target_actions


    # Get enzymes
    drug2attrib[drugbank_id]['enzymes'] = []
    for enzyme in child.find(tag_prefix+'enzymes').findall(tag_prefix+'enzyme'):
        if enzyme.find(tag_prefix+'polypeptide') is None:
            continue
        
        enzyme_id = enzyme.find(tag_prefix+'id').text
        enzyme_gene = enzyme.find(tag_prefix+'polypeptide').find(tag_prefix+'gene-name').text
        enzyme_name = enzyme.find(tag_prefix+'name').text
        enzyme_organism = enzyme.find(tag_prefix+'organism').text
        enzyme_taxonomy_id = enzyme.find(tag_prefix+'polypeptide').find(tag_prefix+'organism').attrib['ncbi-taxonomy-id']
        
        if enzyme_organism is None and enzyme_taxonomy_id == '9606':
            enzyme_organism = 'Human'
        if enzyme_organism == 'Human' and enzyme_taxonomy_id == '':
            enzyme_taxonomy_id = '9606'
        if enzyme_organism == 'Homo sapiens':
            enzyme_organism = 'Human'
        if enzyme_gene is None or enzyme_organism is None or enzyme_taxonomy_id is None:
            continue
        
        enzyme_external_ids = enzyme.find(tag_prefix+'polypeptide').find(tag_prefix+'external-identifiers').findall(tag_prefix+'external-identifier')
        enzyme_uniprot_id = ''
        enzyme_genbank_gene = ''
        enzyme_genbank_protein = ''
        enzyme_hgnc_id = ''

        for external_id in enzyme_external_ids:
            if external_id.find(tag_prefix+'resource').text == 'UniProtKB':
                enzyme_uniprot_id = external_id.find(tag_prefix+'identifier').text
            elif external_id.find(tag_prefix+'resource').text == 'GenBank Gene Database':
                enzyme_genbank_gene = external_id.find(tag_prefix+'identifier').text
            elif external_id.find(tag_prefix+'resource').text == 'GenBank Protein Database':
                enzyme_genbank_protein = external_id.find(tag_prefix+'identifier').text
            elif external_id.find(tag_prefix+'resource').text == 'HUGO Gene Nomenclature Committee (HGNC)':
                enzyme_hgnc_id = external_id.find(tag_prefix+'identifier').text

        enzyme_actions = [s.text.lower() for s in enzyme.find(tag_prefix+'actions').findall(tag_prefix+'action')]
        
        drug2attrib[drugbank_id]['enzymes'].append((enzyme_id,enzyme_actions))

        if enzyme_id not in enzyme2attrib: #{'gene', 'name', 'organism', 'taxonomy_id', 'uniprot_id', 'genbank_gene_id', 'genbank_protein_id', 'hgnc_id'}
            enzyme2attrib[enzyme_id]['gene'] = enzyme_gene
            enzyme2attrib[enzyme_id]['name'] = enzyme_name
            enzyme2attrib[enzyme_id]['organism'] = enzyme_organism
            enzyme2attrib[enzyme_id]['taxonomy_id'] = enzyme_taxonomy_id
            
            enzyme2attrib[enzyme_id]['uniprot_id'] = enzyme_uniprot_id
            enzyme2attrib[enzyme_id]['genbank_gene_id'] = enzyme_genbank_gene
            enzyme2attrib[enzyme_id]['genbank_protein_id'] = enzyme_genbank_protein
            enzyme2attrib[enzyme_id]['hgnc_id'] = enzyme_hgnc_id

        #print enzyme_id, enzyme_gene, enzyme_name, enzyme_organism, enzyme_taxonomy_id, enzyme_actions


    # Get transporters
    drug2attrib[drugbank_id]['transporters'] = []
    for transporter in child.find(tag_prefix+'transporters').findall(tag_prefix+'transporter'):
        if transporter.find(tag_prefix+'polypeptide') is None:
            continue
        
        transporter_id = transporter.find(tag_prefix+'id').text
        transporter_gene = transporter.find(tag_prefix+'polypeptide').find(tag_prefix+'gene-name').text
        transporter_name = transporter.find(tag_prefix+'name').text
        transporter_organism = transporter.find(tag_prefix+'organism').text
        transporter_taxonomy_id = transporter.find(tag_prefix+'polypeptide').find(tag_prefix+'organism').attrib['ncbi-taxonomy-id']
        
        if transporter_organism is None and transporter_taxonomy_id == '9606':
            transporter_organism = 'Human'
        if transporter_organism == 'Human' and transporter_taxonomy_id == '':
            transporter_taxonomy_id = '9606'
        if transporter_organism == 'Homo sapiens':
            transporter_organism = 'Human'
        if transporter_gene is None or transporter_organism is None or transporter_taxonomy_id is None:
            continue
        
        transporter_external_ids = transporter.find(tag_prefix+'polypeptide').find(tag_prefix+'external-identifiers').findall(tag_prefix+'external-identifier')
        transporter_uniprot_id = ''
        transporter_genbank_gene = ''
        transporter_genbank_protein = ''
        transporter_hgnc_id = ''

        for external_id in transporter_external_ids:
            if external_id.find(tag_prefix+'resource').text == 'UniProtKB':
                transporter_uniprot_id = external_id.find(tag_prefix+'identifier').text
            elif external_id.find(tag_prefix+'resource').text == 'GenBank Gene Database':
                transporter_genbank_gene = external_id.find(tag_prefix+'identifier').text
            elif external_id.find(tag_prefix+'resource').text == 'GenBank Protein Database':
                transporter_genbank_protein = external_id.find(tag_prefix+'identifier').text
            elif external_id.find(tag_prefix+'resource').text == 'HUGO Gene Nomenclature Committee (HGNC)':
                transporter_hgnc_id = external_id.find(tag_prefix+'identifier').text

        transporter_actions = [s.text.lower() for s in transporter.find(tag_prefix+'actions').findall(tag_prefix+'action')]
        
        drug2attrib[drugbank_id]['transporters'].append((transporter_id,transporter_actions))

        if transporter_id not in transporter2attrib: #{'gene', 'name', 'organism', 'taxonomy_id', 'uniprot_id', 'genbank_gene_id', 'genbank_protein_id', 'hgnc_id'}
            transporter2attrib[transporter_id]['gene'] = transporter_gene
            transporter2attrib[transporter_id]['name'] = transporter_name
            transporter2attrib[transporter_id]['organism'] = transporter_organism
            transporter2attrib[transporter_id]['taxonomy_id'] = transporter_taxonomy_id
            
            transporter2attrib[transporter_id]['uniprot_id'] = transporter_uniprot_id
            transporter2attrib[transporter_id]['genbank_gene_id'] = transporter_genbank_gene
            transporter2attrib[transporter_id]['genbank_protein_id'] = transporter_genbank_protein
            transporter2attrib[transporter_id]['hgnc_id'] = transporter_hgnc_id

        #print transporter_id, transporter_gene, transporter_name, transporter_organism, transporter_taxonomy_id, transporter_actions


    print drugbank_id, drugname, drug_type, groups,
    print 'targets:', len(drug2attrib[drugbank_id]['targets']),
    print 'enzymes:', len(drug2attrib[drugbank_id]['enzymes']),
    print 'transporters:', len(drug2attrib[drugbank_id]['transporters'])

print '\n'


#######################################################################
# List of drugs to save (as long as num_targets = num_enzymes = num_transporters = 0)
drugs = []
for drugbank_id in sorted(drug2attrib.keys()):
    if len(drug2attrib[drugbank_id]['targets']) == 0 and len(drug2attrib[drugbank_id]['enzymes']) == 0 and len(drug2attrib[drugbank_id]['transporters']) == 0:
        continue
    else:
        drugs.append(drugbank_id)
        
print len(drug2attrib), "drugs parsed from XML"
print len(drugs), "drugs with at least 1 target/ enzyme/ transporter"


#######################################################################
# Save drug attributes to CSV {'drugname', 'drug_type', 'groups', 'targets/enzymes/transporters': [_id, _actions]}
longest_drugname = ''
outf = open('drugbank05_drugs.csv', 'w')
writer = csv.writer(outf)
writer.writerow(['drugbank_id', 'drugname', 'drug_type', 'approved', 'experimental', 'illicit', 'investigational', 'nutraceutical', 'withdrawn'])

for drugbank_id in drugs:
    drugname = drug2attrib[drugbank_id]['drugname']
    if isinstance(drugname, unicode):
        if u'\u03b2' in drugname:
            drugname = drugname.replace(u'\u03b2', 'beta')
        if u'\u03b1' in drugname:
            drugname = drugname.replace(u'\u03b1', 'alpha')
        drugname = drugname.encode("utf-8")
    drug_type = drug2attrib[drugbank_id]['drug_type']
    groups = [1 if group in drug2attrib[drugbank_id]['groups'] else 0 for group in ['approved', 'experimental', 'illicit', 'investigational', 'nutraceutical', 'withdrawn'] ]
    
    writer.writerow([drugbank_id, drugname, drug_type]+groups)
    
    if len(drugname) > len(longest_drugname):
        longest_drugname = drugname
    
outf.close()


#######################################################################
# Save all targets, enzymes, transporters to CSV

# drugbank_target_id -> #{'gene', 'name', 'organism', 'taxonomy_id', 'uniprot_id', 'genbank_gene_id', 'genbank_protein_id', 'hgnc_id'}

outf = open('drugbank05_partner_protein.csv', 'w')
outfh = open('drugbank05_partner_protein_human.csv', 'w')
writer = csv.writer(outf)
writerh = csv.writer(outfh)

writer.writerow(['partner_id', 'partner_name', 'gene_name', 'uniprot_id', 'genbank_gene_id', 'genbank_protein_id', 'hgnc_id', 'organism', 'taxonomy_id'])
writerh.writerow(['partner_id', 'partner_name', 'gene_name', 'uniprot_id', 'genbank_gene_id', 'genbank_protein_id', 'hgnc_id', 'organism', 'taxonomy_id'])

partners_written = set()

# Targets
for partner_id in sorted(target2attrib.keys()):
    if partner_id in partners_written:
        # print partner_id, target2attrib[partner_id]['gene'], 'already recorded'
        continue
    
    partner_name = target2attrib[partner_id]['name']
    gene_name = target2attrib[partner_id]['gene']
    organism = target2attrib[partner_id]['organism']
    taxonomy_id = target2attrib[partner_id]['taxonomy_id']
    
    uniprot_id = target2attrib[partner_id]['uniprot_id']
    genbank_gene_id = target2attrib[partner_id]['genbank_gene_id']
    genbank_protein_id = target2attrib[partner_id]['genbank_protein_id']
    hgnc_id = target2attrib[partner_id]['hgnc_id']
    
    partners_written.add(partner_id)
    
    writer.writerow([partner_id, partner_name, gene_name, uniprot_id, genbank_gene_id, genbank_protein_id, hgnc_id, organism, taxonomy_id])
    
    if taxonomy_id == '9606' and organism == 'Human':
        writerh.writerow([partner_id, partner_name, gene_name, uniprot_id, genbank_gene_id, genbank_protein_id, hgnc_id, organism, taxonomy_id])
    
    if taxonomy_id == '9606' and organism.lower() != 'human':
        print partner_id, target2attrib[partner_id]['gene'], organism, taxonomy_id, 'organism mismatch'
        
        
# enzymes
for partner_id in sorted(enzyme2attrib.keys()):
    if partner_id in partners_written:
        # print partner_id, enzyme2attrib[partner_id]['gene'], 'already recorded in targets'
        continue
    
    partner_name = enzyme2attrib[partner_id]['name']
    gene_name = enzyme2attrib[partner_id]['gene']
    organism = enzyme2attrib[partner_id]['organism']
    taxonomy_id = enzyme2attrib[partner_id]['taxonomy_id']
    
    uniprot_id = enzyme2attrib[partner_id]['uniprot_id']
    genbank_gene_id = enzyme2attrib[partner_id]['genbank_gene_id']
    genbank_protein_id = enzyme2attrib[partner_id]['genbank_protein_id']
    hgnc_id = enzyme2attrib[partner_id]['hgnc_id']
    
    partners_written.add(partner_id)
    
    writer.writerow([partner_id, partner_name, gene_name, uniprot_id, genbank_gene_id, genbank_protein_id, hgnc_id, organism, taxonomy_id])
    
    if taxonomy_id == '9606' and organism == 'Human':
        writerh.writerow([partner_id, partner_name, gene_name, uniprot_id, genbank_gene_id, genbank_protein_id, hgnc_id, organism, taxonomy_id])
    
    if taxonomy_id == '9606' and organism.lower() != 'human':
        print partner_id, enzyme2attrib[partner_id]['gene'], organism, taxonomy_id, 'organism mismatch'
        
        
# transporters
for partner_id in sorted(transporter2attrib.keys()):
    if partner_id in partners_written:
        # print partner_id, transporter2attrib[partner_id]['gene'], 'already recorded in targets and/or enzymes'
        continue
    
    partner_name = transporter2attrib[partner_id]['name']
    gene_name = transporter2attrib[partner_id]['gene']
    organism = transporter2attrib[partner_id]['organism']
    taxonomy_id = transporter2attrib[partner_id]['taxonomy_id']
    
    uniprot_id = transporter2attrib[partner_id]['uniprot_id']
    genbank_gene_id = transporter2attrib[partner_id]['genbank_gene_id']
    genbank_protein_id = transporter2attrib[partner_id]['genbank_protein_id']
    hgnc_id = transporter2attrib[partner_id]['hgnc_id']
    
    partners_written.add(partner_id)
    
    writer.writerow([partner_id, partner_name, gene_name, uniprot_id, genbank_gene_id, genbank_protein_id, hgnc_id, organism, taxonomy_id])
    
    if taxonomy_id == '9606' and organism == 'Human':
        writerh.writerow([partner_id, partner_name, gene_name, uniprot_id, genbank_gene_id, genbank_protein_id, hgnc_id, organism, taxonomy_id])
    
    if taxonomy_id == '9606' and organism.lower() != 'human':
        print partner_id, transporter2attrib[partner_id]['gene'], organism, taxonomy_id, 'organism mismatch'
        
outf.close()
outfh.close()


#######################################################################
# Save drug-target, -enzyme, -transporter pairs to CSV

# target [('antagonist', 1374), ('agonist', 857), ('inhibitor', 1818)]
# enzyme [('substrate', 2402), ('inducer', 407), ('inhibitor', 1350)]
# transporter [('substrate', 790), ('inducer', 100), ('inhibitor', 1075)]

# Targets
outf = open('drugbank05_drug2target.csv', 'w')
outfh = open('drugbank05_drug2target_human.csv', 'w')
writer = csv.writer(outf)
writerh = csv.writer(outfh)

target_actions_to_write = ['inhibitor', 'antagonist', 'agonist']
writer.writerow(['drugbank_id', 'partner_id']+target_actions_to_write)
writerh.writerow(['drugbank_id', 'partner_id']+target_actions_to_write)

for drugbank_id in drugs:
    for (target_id, target_actions) in drug2attrib[drugbank_id]['targets']:
        actions = [1 if action in target_actions else 0 for action in target_actions_to_write]
        
        writer.writerow([drugbank_id, target_id]+actions)
        
        if target2attrib[target_id]['organism'] == 'Human' and target2attrib[target_id]['taxonomy_id'] == '9606':
            writerh.writerow([drugbank_id, target_id]+actions)

outf.close()
outfh.close()


# Enzymes
outf = open('drugbank05_drug2enzyme.csv', 'w')
outfh = open('drugbank05_drug2enzyme_human.csv', 'w')
writer = csv.writer(outf)
writerh = csv.writer(outfh)

enzyme_actions_to_write = ['substrate', 'inducer', 'inhibitor']
writer.writerow(['drugbank_id', 'partner_id']+enzyme_actions_to_write)
writerh.writerow(['drugbank_id', 'partner_id']+enzyme_actions_to_write)

for drugbank_id in drugs:
    for (enzyme_id, enzyme_actions) in drug2attrib[drugbank_id]['enzymes']:
        actions = [1 if action in enzyme_actions else 0 for action in enzyme_actions_to_write]
        
        writer.writerow([drugbank_id, enzyme_id]+actions)
        
        if enzyme2attrib[enzyme_id]['organism'] == 'Human' and enzyme2attrib[enzyme_id]['taxonomy_id'] == '9606':
            writerh.writerow([drugbank_id, enzyme_id]+actions)

outf.close()
outfh.close()


# Transporters
outf = open('drugbank05_drug2transporter.csv', 'w')
outfh = open('drugbank05_drug2transporter_human.csv', 'w')
writer = csv.writer(outf)
writerh = csv.writer(outfh)

transporter_actions_to_write = ['substrate', 'inducer', 'inhibitor']
writer.writerow(['drugbank_id', 'partner_id']+transporter_actions_to_write)
writerh.writerow(['drugbank_id', 'partner_id']+transporter_actions_to_write)

for drugbank_id in drugs:
    for (transporter_id, transporter_actions) in drug2attrib[drugbank_id]['transporters']:
        actions = [1 if action in transporter_actions else 0 for action in transporter_actions_to_write]
        
        writer.writerow([drugbank_id, transporter_id]+actions)
        
        if transporter2attrib[transporter_id]['organism'] == 'Human' and transporter2attrib[transporter_id]['taxonomy_id'] == '9606':
            writerh.writerow([drugbank_id, transporter_id]+actions)

outf.close()
outfh.close()

